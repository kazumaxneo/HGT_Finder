"""GTDB taxonomy and DIAMOND inputs. No network calls or automatic downloads."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import itertools
from pathlib import Path
from typing import TextIO

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")
BLAST_COLUMNS = ("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    p = Path(path)
    return gzip.open(p, mode, encoding="utf-8") if p.suffix == ".gz" else p.open(mode.replace("t", ""), encoding="utf-8")


def canonical_accession(value: str) -> str:
    value = value.strip()
    return value[3:] if value.startswith(("RS_", "GB_")) else value


def lineage(value: str) -> dict[str, str]:
    parts = value.strip().split(";")
    if len(parts) != 7 or any(not v.startswith(p) for v, p in zip(parts, PREFIXES)):
        raise ValueError("Taxonomy must contain seven GTDB ranks: d__;p__;c__;o__;f__;g__;s__")
    return {rank: part if len(part) > 3 else "" for rank, part in zip(RANKS, parts)}


def tax_key(tax: dict[str, str], rank: str) -> tuple[str, ...] | None:
    # Full lineage avoids conflating identically named taxa in different parents.
    key = tuple(tax[r] for r in RANKS[:RANKS.index(rank) + 1])
    return key if all(key) else None


def read_taxonomy(path: str | Path) -> dict[str, dict[str, str]]:
    result = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader, None)
        if first is None:
            raise ValueError("Empty taxonomy table")
        if "gtdb_taxonomy" in first or "classification" in first:
            col = first.index("gtdb_taxonomy") if "gtdb_taxonomy" in first else first.index("classification")
            acc = first.index("accession") if "accession" in first else first.index("user_genome") if "user_genome" in first else 0
            rows = reader
        else:
            col, acc = 1, 0
            rows = itertools.chain([first], reader)
        for row in rows:
            if not row or not row[0] or row[0].startswith("#"):
                continue
            if len(row) <= max(col, acc):
                raise ValueError(f"Malformed taxonomy row: {row[0]}")
            if not row[col].startswith("d__"):
                # Header of an exported two-column table.
                if row[acc] in {"genome", "genome_id", "accession", "user_genome"}:
                    continue
                raise ValueError(f"Missing GTDB lineage for {row[acc]}")
            key = canonical_accession(row[acc])
            value = lineage(row[col])
            if key in result and result[key] != value:
                raise ValueError(f"Conflicting taxonomy for {key}")
            result[key] = value
    if not result:
        raise ValueError("No usable taxonomy rows")
    return result


def write_taxonomy(path: Path, taxa: dict[str, dict[str, str]]) -> None:
    with path.open("w") as f:
        for genome, tax in sorted(taxa.items()):
            f.write(genome + "\t" + ";".join(tax[r] or PREFIXES[i] for i, r in enumerate(RANKS)) + "\n")


def tsv_write(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if row.get(k) is None else row.get(k, "") for k in columns})


def read_genes(path: Path) -> list[dict]:
    with open_text(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    required = {"orf_id", "contig", "start", "end"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Genes TSV needs orf_id, contig, start, end columns")
    seen = set()
    for row in rows:
        if row["orf_id"] in seen:
            raise ValueError(f"Duplicate ORF: {row['orf_id']}")
        seen.add(row["orf_id"])
        row["start"], row["end"] = int(row["start"]), int(row["end"])
        if row["start"] < 1 or row["end"] < row["start"]:
            raise ValueError(f"Invalid coordinates for {row['orf_id']}")
        row.setdefault("strand", ".")
        row.setdefault("gene", "")
        row.setdefault("product", "")
    return sorted(rows, key=lambda r: (r["contig"], r["start"], r["end"], r["orf_id"]))


def prepare_genbank(source: Path, out: Path) -> dict:
    from Bio import SeqIO
    records = []
    genes = []
    excluded = []
    seen = set()
    with open_text(source) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for index, feature in enumerate(record.features):
                if feature.type != "CDS":
                    continue
                q = feature.qualifiers
                orf = q.get("locus_tag", q.get("protein_id", [f"{record.id}_{index}"]))[0]
                if orf in seen:
                    raise ValueError(f"Duplicate CDS identifier: {orf}")
                seen.add(orf)
                row = {"orf_id": orf, "contig": record.id, "start": int(feature.location.start) + 1,
                       "end": int(feature.location.end), "strand": feature.location.strand or 0,
                       "gene": q.get("gene", [""])[0], "product": q.get("product", [""])[0]}
                if "translation" not in q:
                    excluded.append({**row, "reason": "no_translation"})
                    continue
                seq = "".join(q["translation"][0].split()).rstrip("*")
                if not seq:
                    excluded.append({**row, "reason": "empty_translation"})
                    continue
                records.append((orf, seq))
                genes.append({**row, "protein_length": len(seq)})
    if not genes:
        raise ValueError("No translated CDS in GenBank")
    out.mkdir(parents=True, exist_ok=True)
    with (out / "query.faa").open("w") as f:
        for orf, seq in records:
            f.write(f">{orf}\n{seq}\n")
    tsv_write(out / "genes.tsv", genes)
    tsv_write(out / "excluded_cds.tsv", excluded,
              ["orf_id", "contig", "start", "end", "strand", "gene", "product", "reason"])
    return {"translated_cds": len(genes), "excluded_cds": len(excluded)}


def read_subject_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or not {"sseqid", "genome_id"}.issubset(reader.fieldnames):
            raise ValueError("Subject map needs sseqid and genome_id columns")
        result = {}
        for r in reader:
            genome = canonical_accession(r["genome_id"])
            if r["sseqid"] in result and result[r["sseqid"]] != genome:
                raise ValueError(f"Conflicting subject map: {r['sseqid']}")
            result[r["sseqid"]] = genome
        return result


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def dump_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
