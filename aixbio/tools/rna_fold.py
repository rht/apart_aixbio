from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

# Base-pair energies (kcal/mol estimates for nearest-neighbor stacking)
_PAIR_ENERGY = {
    "AU": -0.9, "UA": -0.9,
    "GC": -2.1, "CG": -2.1,
    "GU": -0.5, "UG": -0.5,
}

_MIN_LOOP = 3  # Minimum hairpin loop size


def _can_pair(a: str, b: str) -> bool:
    return (a + b) in _PAIR_ENERGY


def estimate_five_prime_dg(dna: str, window: int = 50) -> float:
    """Estimate minimum free energy of 5' region using Nussinov-style DP.

    Uses a simplified dynamic programming approach (Nussinov algorithm) to find
    the maximum number of base pairs in the 5' window, then estimates ΔG from
    paired bases.

    .. warning::
        This is a heuristic; for production use, integrate ViennaRNA / RNAfold.
        The energy values are rough approximations and do not account for
        stacking, dangling ends, or loop penalties beyond minimum loop size.

    Returns a negative kcal/mol estimate (more negative = more stable = worse
    for expression).
    """
    logger.debug(
        "RNA fold: using Nussinov DP heuristic (not ViennaRNA). "
        "Results are approximate."
    )
    seq = dna[:window].upper().replace("T", "U")
    n = len(seq)
    if n < 10:
        return 0.0

    # Nussinov DP: dp[i][j] = max base pairs in seq[i..j]
    dp = [[0] * n for _ in range(n)]

    for span in range(_MIN_LOOP + 2, n):
        for i in range(n - span):
            j = i + span
            # Case 1: j is unpaired
            dp[i][j] = dp[i][j - 1]
            # Case 2: j pairs with some k
            for k in range(i, j - _MIN_LOOP):
                if _can_pair(seq[k], seq[j]):
                    score = 1 + (dp[i][k - 1] if k > i else 0) + (dp[k + 1][j - 1] if k + 1 <= j - 1 else 0)
                    dp[i][j] = max(dp[i][j], score)

    # Estimate ΔG from paired bases using average pair energy
    max_pairs = dp[0][n - 1]
    # Use average energy per pair (~-1.5 kcal/mol, weighted toward GC)
    avg_pair_energy = -1.5
    dg = max_pairs * avg_pair_energy

    return dg
