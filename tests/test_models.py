from aixbio.models import (
    Chain,
    ProteinRecord,
    DNAChain,
    CassetteElement,
    CassetteChain,
    CheckResult,
    ChainValidation,
    ValidationReport,
    RemediationAction,
    PlannedFix,
    RemediationPlan,
    AgentDecision,
)


def test_chain_frozen():
    c = Chain(id="A", aa_sequence="MVLS", length=4)
    assert c.id == "A"
    assert c.length == 4
    try:
        c.id = "B"
        assert False, "Should be frozen"
    except AttributeError:
        pass


def test_protein_record():
    chains = (
        Chain(id="A", aa_sequence="GIVEQCCTSICSLYQLENYCN", length=21),
        Chain(id="B", aa_sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKT", length=30),
    )
    pr = ProteinRecord(uniprot_id="P01308", name="Insulin", chains=chains)
    assert len(pr.chains) == 2
    assert pr.chains[0].id == "A"


def test_validation_report():
    checks = (
        CheckResult(name="gc_content", passed=True, value=0.55, threshold="0.50-0.60"),
        CheckResult(name="cai_score", passed=True, value=0.85, threshold="> 0.8"),
    )
    cv = ChainValidation(id="test", passed=True, checks=checks)
    vr = ValidationReport(chains=(cv,), all_passed=True)
    assert vr.all_passed


def test_remediation_plan():
    fix = PlannedFix(
        check_name="restriction_sites",
        strategy="synonymous_swap",
        target_positions=(5,),
        replacement_codons=("GCG",),
    )
    plan = RemediationPlan(
        actions=(fix,),
        reasoning="Break BamHI site at position 15",
        priority_order=("restriction_sites",),
    )
    assert len(plan.actions) == 1
