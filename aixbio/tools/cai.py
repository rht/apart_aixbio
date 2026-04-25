from __future__ import annotations

import math

from aixbio.tools.codon_tables import get_ecoli_table, split_codons, translate_codon


def compute_cai(dna: str) -> float:
    table = get_ecoli_table()
    codons = split_codons(dna.upper())
    if not codons:
        return 0.0

    log_sum = 0.0
    count = 0
    for codon in codons:
        aa = translate_codon(codon)
        if aa == "*":
            continue
        aa_codons = table.get(aa)
        if aa_codons is None:
            continue
        max_freq = max(aa_codons.values())
        codon_freq = aa_codons.get(codon, 0.0)
        if max_freq > 0 and codon_freq > 0:
            log_sum += math.log(codon_freq / max_freq)
            count += 1

    if count == 0:
        return 0.0
    return math.exp(log_sum / count)
