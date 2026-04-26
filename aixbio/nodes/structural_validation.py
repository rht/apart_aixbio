"""DNA-level validation using Evo2 genomic language model.

Scores each chain's optimized DNA sequence via Evo2 log-probability to assess
whether the codon-optimized construct looks like plausible, expressible DNA —
a genuine DNA-level check, unlike protein fold predictions which are identical
for synonymous codon variants.
"""
from __future__ import annotations

import logging

from aixbio.models.structure import Evo2Report
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.evo2 import score_dna

logger = logging.getLogger(__name__)


def structural_validation(state: PipelineState) -> dict:
    chain_results = state.get("chain_results", [])
    if not chain_results:
        return {"structure_report": None}

    results = []
    warnings: list[str] = []
    for cr in chain_results:
        chain_id = cr["chain_id"]
        dna_seq = cr.get("optimized_dna_sequence") or cr.get("cassette_dna", "")
        if not dna_seq:
            msg = f"No DNA sequence available for chain {chain_id}, skipping Evo2"
            logger.warning(msg)
            warnings.append(msg)
            continue
        evo_result = score_dna(chain_id, dna_seq)
        if evo_result.method == "evo2_truncated":
            warnings.append(
                f"Chain {chain_id}: sequence truncated to 4096 bp for Evo2 scoring"
            )
        results.append(evo_result)

    if not results:
        return {"structure_report": None, "warnings": warnings}

    return {"structure_report": Evo2Report(chains=tuple(results)), "warnings": warnings}
