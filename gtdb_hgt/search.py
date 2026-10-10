from __future__ import annotations

from contextlib import contextmanager
import gzip
import io
import re
import subprocess
import sys
import tarfile
import time
from pathlib import Path

from Bio import SeqIO

from . import __version__
from .io import BLAST_COLUMNS, canonical_accession, dump_json, open_text, read_taxonomy, write_taxonomy


SENSITIVITY_MODES = ("faster", "fast", "default", "mid-sensitive", "sensitive",
                     "more-sensitive", "very-sensitive", "ultra-sensitive")


def execute(cmd: list[str], logfile: Path):
    logfile.parent.mkdir(parents=True, exist_ok=True)
    with logfile.open("a") as f:
        f.write("COMMAND: " + repr(cmd) + "\n")
        f.flush()
        subprocess.run(cmd, check=True, stdout=f, stderr=f)


def protein_sources(path: Path):
    if path.is_dir():
        files = sorted(p for p in path.rglob("*") if p.is_file() and p.name.endswith((".faa", ".faa.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")))
        if not files:
            raise ValueError("No protein FASTA files in directory")
        for p in files:
            with open_text(p) as f:
                yield p.name, f
    elif tarfile.is_tarfile(path):
        # Read members in place. Never extract archive paths to the filesystem.
        with tarfile.open(path, "r|*") as archive:
            for member in archive:
                if not member.isfile() or not member.name.endswith((".faa", ".faa.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")):
                    continue
                raw = archive.extractfile(member)
                if raw is None:
                    continue
                stream = gzip.GzipFile(fileobj=raw) if member.name.endswith(".gz") else raw
                with io.TextIOWrapper(stream) as text:
                    yield Path(member.name).name, text
    else:
        with open_text(path) as f:
            yield path.name, f


def filename_genome(name: str, taxa):
    match = re.search(r"(?:RS_|GB_)?GC[AF]_\d+\.\d+", name)
    if match:
        return canonical_accession(match[0])
    for suffix in (".gz", ".faa", ".fasta", ".fa"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    for suffix in ("_protein", "_proteins", "_translated", "_genomic"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    key = canonical_accession(name)
    return key if key in taxa else None


def build_db(proteins: Path, taxonomy: Path, out: Path, diamond: str, threads: int, label: str):
    started = time.perf_counter()
    if out.exists() and any(out.iterdir()):
        raise ValueError("DB output directory must be empty; existing DBs are never overwritten")
    out.mkdir(parents=True, exist_ok=True)
    taxa = read_taxonomy(taxonomy)
    used, n, letters, skipped = {}, 0, 0, []
    fasta = out / "proteins.faa"
    with fasta.open("w") as dest:
        for name, handle in protein_sources(proteins):
            genome = filename_genome(name, taxa)
            if genome is None or genome not in taxa:
                skipped.append(name)
                continue
            if genome in used:
                raise ValueError(f"Multiple input protein files for {genome}")
            used[genome] = taxa[genome]
            seen = set()
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in seen or "::" in record.id:
                    raise ValueError(f"Duplicate or reserved protein ID in {name}: {record.id}")
                seen.add(record.id)
                dest.write(f">{genome}::{record.id}\n{record.seq}\n")
                n += 1
                letters += len(record.seq)
            if len(used) % 1000 == 0:
                print(f"Prepared {len(used):,} genomes; {n:,} proteins; {letters:,} amino acids; "
                      f"{time.perf_counter() - started:.1f} seconds", file=sys.stderr, flush=True)
    if not n:
        raise ValueError("No proteins matched taxonomy. Supply one FAA per genome with accession in its filename.")
    write_taxonomy(out / "taxonomy.tsv", used)
    execute([diamond, "makedb", "--in", str(fasta), "--db", str(out / "proteins"), "--threads", str(threads)], out / "build.log")
    manifest = {"tool_version": __version__, "label": label, "source_proteins": str(proteins.resolve()),
                "source_taxonomy": str(taxonomy.resolve()), "genomes": len(used), "proteins": n,
                "letters": letters, "build_wall_seconds": time.perf_counter() - started,
                "skipped_files": skipped, "id_encoding": "genome_accession::original_protein_id"}
    dump_json(out / "database.json", manifest)
    return manifest


def search(query: Path, db: Path, out: Path, diamond: str, threads: int,
           evalue=1e-10, max_target_seqs=0, sensitivity="very-sensitive"):
    if sensitivity not in SENSITIVITY_MODES:
        raise ValueError(f"Unknown DIAMOND sensitivity: {sensitivity}")
    if out.exists():
        raise ValueError(f"Search output already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [diamond, "blastp", "--query", str(query), "--db", str(db), "--out", str(out),
           "--threads", str(threads), "--evalue", str(evalue), "--max-target-seqs", str(max_target_seqs),
           "--outfmt", "6", *BLAST_COLUMNS]
    # DIAMOND default sensitivity is selected by omitting all sensitivity flags.
    if sensitivity != "default":
        cmd += ["--" + sensitivity]
    if out.suffix == ".gz":
        cmd += ["--compress", "1"]
    execute(cmd, out.with_suffix(out.suffix + ".log"))
    # A cap of zero only removes the target count limit. Sensitivity, masking,
    # significance and aligner heuristics still determine which hits are found.
    manifest = {"tool_version": __version__, "command": cmd, "max_target_seqs": max_target_seqs,
                "evalue": evalue, "sensitivity": sensitivity, "query": str(query.resolve()),
                "database": str(db.resolve()), "columns": list(BLAST_COLUMNS)}
    dump_json(out.with_suffix(out.suffix + ".json"), manifest)
    return manifest


def search_self(query: Path, out: Path, diamond: str, threads: int, sensitivity: str):
    prefix = out.parent / "self_database"
    execute([diamond, "makedb", "--in", str(query), "--db", str(prefix)], out.parent / "self_search.log")
    return search(query, prefix.with_suffix(".dmnd"), out, diamond, threads,
                  evalue=1e-10, max_target_seqs=0, sensitivity=sensitivity)
