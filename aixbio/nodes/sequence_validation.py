from __future__ import annotations

from aixbio.models.validation import ChainValidation, CheckResult
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.cai import compute_cai
from aixbio.tools.codon_tables import RARE_CODONS_ECOLI, split_codons, translate_dna
from aixbio.tools.gc import compute_gc
from aixbio.tools.restriction_sites import find_restriction_sites
from aixbio.tools.rna_fold import estimate_five_prime_dg


def sequence_validation(state: ChainSubgraphState) -> dict:
    dna_chain = state["optimized_dna"]
    chain = state["chain"]
    cloning_sites = state["cloning_sites"]
    dna = dna_chain.dna_sequence

    checks: list[CheckResult] = []

    gc = compute_gc(dna)
    checks.append(CheckResult(
        name="gc_content",
        passed=0.50 <= gc <= 0.60,
        value=round(gc, 4),
        threshold="0.50-0.60",
    ))

    cai = compute_cai(dna)
    checks.append(CheckResult(
        name="cai_score",
        passed=cai > 0.8,
        value=round(cai, 4),
        threshold="> 0.8",
    ))

    sites = find_restriction_sites(dna, cloning_sites)
    checks.append(CheckResult(
        name="restriction_sites",
        passed=len(sites) == 0,
        value=str(sites) if sites else "none",
        threshold="0 hits",
    ))

    dg = estimate_five_prime_dg(dna)
    checks.append(CheckResult(
        name="rna_secondary_structure",
        passed=dg > -10.0,
        value=round(dg, 2),
        threshold="> -10 kcal/mol",
    ))

    translated = translate_dna(dna)
    expected = chain.aa_sequence
    identity_match = translated.rstrip("*") == expected
    checks.append(CheckResult(
        name="back_translation",
        passed=identity_match,
        value="match" if identity_match else f"mismatch at position {_first_mismatch(translated, expected)}",
        threshold="100% identity",
    ))

    codons = split_codons(dna)
    rare_count = sum(1 for c in codons if c.upper() in RARE_CODONS_ECOLI)
    checks.append(CheckResult(
        name="rare_codons",
        passed=rare_count == 0,
        value=rare_count,
        threshold="0",
    ))

    all_passed = all(c.passed for c in checks)
    failed = tuple(c for c in checks if not c.passed)

    validation = ChainValidation(
        id=chain.id,
        passed=all_passed,
        checks=tuple(checks),
    )
    return {"chain_validation": validation, "failed_checks": failed}


def _first_mismatch(a: str, b: str) -> int:
    for i, (ca, cb) in enumerate(zip(a, b)):
        if ca != cb:
            return i
    return min(len(a), len(b))
