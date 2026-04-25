from __future__ import annotations

from aixbio.models.validation import ChainValidation, ValidationReport
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.state.pipeline_state import ChainProcessingResult, PipelineState


def package_result(state: ChainSubgraphState) -> dict:
    return {
        "chain_results": [ChainProcessingResult(
            chain_id=state["chain"].id,
            optimized_dna_sequence=state["optimized_dna"].dna_sequence if state["optimized_dna"] else "",
            cassette_dna=state["cassette"].full_dna if state["cassette"] else "",
            genbank_file=state["plasmid"].genbank_file if state["plasmid"] else "",
            insert_size=state["plasmid"].insert_size if state["plasmid"] else 0,
            validation_passed=state["chain_validation"].passed if state["chain_validation"] else False,
            checks=state["chain_validation"].checks if state["chain_validation"] else (),
            remediation_rounds=state["remediation_attempt"],
            remediation_history=tuple(state.get("remediation_history", [])),
        )],
    }


def package_result_failed(state: ChainSubgraphState) -> dict:
    return package_result(state)


def halt_pipeline(state: ChainSubgraphState) -> dict:
    result = package_result(state)
    result["warnings"] = [
        f"HALT: Back-translation check failed for chain {state['chain'].id}. Pipeline bug detected."
    ]
    return result


def merge_all_chain_results(state: PipelineState) -> dict:
    chain_results = state.get("chain_results", [])
    chain_validations = []
    for cr in chain_results:
        chain_validations.append(ChainValidation(
            id=cr["chain_id"],
            passed=cr["validation_passed"],
            checks=cr["checks"],
        ))

    all_passed = all(cv.passed for cv in chain_validations)

    return {
        "validation_report": ValidationReport(
            chains=tuple(chain_validations),
            all_passed=all_passed,
        ),
        "pipeline_status": "completed" if all_passed else "failed",
    }
