from __future__ import annotations

from typing import Literal

from langgraph.types import Send

from aixbio.state.chain_state import ChainSubgraphState
from aixbio.state.pipeline_state import PipelineState


def validation_router(
    state: ChainSubgraphState,
) -> Literal["package_result", "halt_pipeline", "remediation_agent"]:
    validation = state["chain_validation"]
    if validation is None:
        return "halt_pipeline"
    if validation.passed:
        return "package_result"
    for check in validation.checks:
        if check.name == "back_translation" and not check.passed:
            return "halt_pipeline"
    return "remediation_agent"


def revalidation_router(
    state: ChainSubgraphState,
) -> Literal["package_result", "halt_pipeline", "package_result_failed", "remediation_agent"]:
    validation = state["chain_validation"]
    if validation is None:
        return "halt_pipeline"
    if validation.passed:
        return "package_result"
    for check in validation.checks:
        if check.name == "back_translation" and not check.passed:
            return "halt_pipeline"
    if state["remediation_attempt"] >= state["max_remediation_attempts"]:
        return "package_result_failed"
    return "remediation_agent"


def fan_out_to_chains(state: PipelineState) -> list[Send]:
    protein = state["protein_record"]
    if protein is None:
        return []

    sends = []
    for chain in protein.chains:
        chain_state: ChainSubgraphState = {
            "chain": chain,
            "host_organism": state["host_organism"],
            "avoid_sites": state["avoid_sites"],
            "tag_type": state["tag_type"],
            "protease_site": state["protease_site"],
            "vector": state["vector"],
            "cloning_sites": state["cloning_sites"],
            "protein_record": protein,
            "optimized_dna": None,
            "cassette": None,
            "plasmid": None,
            "chain_validation": None,
            "remediation_attempt": 0,
            "max_remediation_attempts": state.get("max_remediation_attempts", 3),
            "failed_checks": (),
            "remediation_plan": None,
            "remediation_history": [],
            "decision_log": [],
            "warnings": [],
        }
        sends.append(Send("chain_processing", chain_state))
    return sends


def structural_router(
    state: PipelineState,
) -> Literal["structural_validation", "__end__"]:
    if state.get("run_structural_validation"):
        return "structural_validation"
    return "__end__"
