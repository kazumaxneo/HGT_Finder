"""Disk-backed query/species aggregation for unrestricted GTDB hit output."""
import sqlite3

from .io import tax_key


class SpeciesHitStore:
    def __init__(self, path, taxa, close_rank, donor_rank, query_tax):
        self.taxa = taxa
        self.genomes = list(taxa)
        self.genome_index = {g: i for i, g in enumerate(self.genomes)}
        self.species_index = {}
        for tax in taxa.values():
            key = tax_key(tax, "species")
            if key not in self.species_index:
                self.species_index[key] = len(self.species_index)
        self.close_rank, self.donor_rank = close_rank, donor_rank
        self.close_key = tax_key(query_tax, close_rank)
        self.db = sqlite3.connect(path)
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA cache_size=-32000")
        self.db.execute("CREATE TABLE hits (q TEXT, species INTEGER, genome INTEGER, subject TEXT, bits REAL, qcov REAL, scov REAL, pident REAL, evalue REAL, PRIMARY KEY(q,species)) WITHOUT ROWID")
        self.pending = []

    def put(self, q, species, row):
        self.pending.append((q, self.species_index[species], self.genome_index[row["genome_id"]],
                             row["sseqid"], row["bitscore"], row["qcov"], row["scov"], row["pident"], row["evalue"]))
        if len(self.pending) >= 20000:
            self.flush()

    def flush(self):
        self.db.executemany("""INSERT INTO hits VALUES (?,?,?,?,?,?,?,?,?)
            ON CONFLICT(q,species) DO UPDATE SET genome=excluded.genome,
            subject=excluded.subject,bits=excluded.bits,qcov=excluded.qcov,scov=excluded.scov,
            pident=excluded.pident,evalue=excluded.evalue
            WHERE (excluded.bits,excluded.qcov,excluded.scov,excluded.subject) >
                  (hits.bits,hits.qcov,hits.scov,hits.subject)""", self.pending)
        self.db.commit()
        self.pending.clear()

    def get(self, q, default=None):
        result = {}
        for species, genome_index, sid, bits, qcov, scov, pident, evalue in self.db.execute(
                "SELECT species,genome,subject,bits,qcov,scov,pident,evalue FROM hits WHERE q=?", (q,)):
            genome = self.genomes[genome_index]
            tax = self.taxa[genome]
            result[species] = {"genome_id": genome, "tax": tax, "sseqid": sid,
                               "bitscore": bits, "qcov": qcov, "scov": scov, "pident": pident, "evalue": evalue,
                               "group": "close" if tax_key(tax, self.close_rank) == self.close_key else "distal",
                               "donor_key": tax_key(tax, self.donor_rank)}
        return result if result else (default if default is not None else {})

    def close(self):
        self.flush()
        self.db.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        self.db.close()
