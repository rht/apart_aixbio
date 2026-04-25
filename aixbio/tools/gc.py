from __future__ import annotations


def compute_gc(dna: str) -> float:
    if not dna:
        return 0.0
    dna = dna.upper()
    gc_count = dna.count("G") + dna.count("C")
    return gc_count / len(dna)
