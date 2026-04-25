from __future__ import annotations


def estimate_five_prime_dg(dna: str, window: int = 50) -> float:
    """Estimate minimum free energy of 5' region.

    Rough heuristic: count potential base pairs in the first `window` nt.
    A full implementation would call RNAfold / mFold.
    Returns a negative kcal/mol estimate (more negative = more stable = worse).
    """
    seq = dna[:window].upper().replace("T", "U")
    if len(seq) < 10:
        return 0.0

    pair_energy = {"AU": -0.9, "UA": -0.9, "GC": -2.1, "CG": -2.1, "GU": -0.5, "UG": -0.5}
    dg = 0.0
    n = len(seq)
    for i in range(n // 2):
        j = n - 1 - i
        if i >= j:
            break
        pair = seq[i] + seq[j]
        dg += pair_energy.get(pair, 0.0)
    return dg
