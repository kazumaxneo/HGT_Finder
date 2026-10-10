"""WAAFLE-inspired protein alignment quality, separate from HGT evidence.

WAAFLE weights nucleotide identity by modified subject coverage, then averages
per-site scores over a locus (uncovered sites are zero). For a single hit to a
full-length protein ORF we use identity * subject coverage * query coverage.
No contig-edge correction is used; this is an adaptation, not WAAFLE itself.
"""
import math


def alignment_quality(pident, qcov, scov):
    values = (float(pident) / 100, float(qcov), float(scov))
    if any(not math.isfinite(v) or not 0 <= v <= 1 for v in values):
        raise ValueError("Alignment quality needs identity 0..100 and coverage fractions 0..1")
    return math.prod(values)


def best_bit_hit(hits):
    """Keep the existing bit-score ranking; quality must not replace the winner."""
    return max(hits, key=lambda h: (h["bitscore"], h.get("qcov", 0),
                                  h.get("scov", 0), h.get("sseqid", ""))) if hits else None


def hit_quality(hit):
    if hit is None or any(hit.get(k) is None for k in ("pident", "qcov", "scov")):
        return None
    return alignment_quality(hit["pident"], hit["qcov"], hit["scov"])
