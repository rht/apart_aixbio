from __future__ import annotations

import RNA


def estimate_five_prime_dg(dna: str, window: int = 50) -> float:
    """Estimate minimum free energy of the 5' mRNA region using ViennaRNA.

    Returns a negative kcal/mol value (more negative = more stable = worse
    for expression).
    """
    seq = dna[:window].upper().replace("T", "U")
    if len(seq) < 10:
        return 0.0
    _structure, mfe = RNA.fold(seq)
    return mfe
