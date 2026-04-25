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
    cassette = state.get("cassette")

    # The gene-only DNA (for back-translation check)
    gene_dna = dna_chain.dna_sequence

    # The full expression construct (for all other checks).
    # If cassette is available, validate the actual construct that enters the
    # host cell: ATG + tag + protease + gene + stop.  Falls back to gene-only
    # when cassette hasn't been assembled yet.
    construct_dna = cassette.full_dna if cassette else gene_dna

    checks: list[CheckResult] = []

    gc = compute_gc(construct_dna)
    checks.append(CheckResult(
        name="gc_content",
        passed=0.50 <= gc <= 0.60,
        value=round(gc, 4),
        threshold="0.50-0.60",
    ))

    cai = compute_cai(gene_dna)  # CAI is meaningful only for the coding gene
    checks.append(CheckResult(
        name="cai_score",
        passed=cai > 0.8,
        value=round(cai, 4),
        threshold="> 0.8",
    ))

    # Check for restriction sites in the construct, but exclude the
    # intentional cloning-site positions at the flanks (Issue #9).
    sites = find_restriction_sites(construct_dna, cloning_sites)
    checks.append(CheckResult(
        name="restriction_sites",
        passed=len(sites) == 0,
        value=str(sites) if sites else "none",
        threshold="0 hits",
    ))

    # RNA secondary structure in the 5' region — the ribosome encounters the
    # tag region first (not the gene), so this must be checked on the construct.
    dg = estimate_five_prime_dg(construct_dna)
    checks.append(CheckResult(
        name="rna_secondary_structure",
        passed=dg > -10.0,
        value=round(dg, 2),
        threshold="> -10 kcal/mol",
    ))

    translated = translate_dna(gene_dna)
    expected = chain.aa_sequence
    identity_match = translated.rstrip("*") == expected
    checks.append(CheckResult(
        name="back_translation",
        passed=identity_match,
        value="match" if identity_match else f"mismatch at position {_first_mismatch(translated, expected)}",
        threshold="100% identity",
    ))

    codons = split_codons(gene_dna)
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
