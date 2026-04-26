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
    for cr in chain_results:
        chain_id = cr["chain_id"]
        dna_seq = cr.get("cassette_dna") or cr.get("optimized_dna_sequence", "")
        if not dna_seq:
            logger.warning(f"No DNA sequence available for chain {chain_id}, skipping Evo2")
            continue
        results.append(score_dna(chain_id, dna_seq))

    if not results:
        return {"structure_report": None}

    return {"structure_report": Evo2Report(chains=tuple(results))}
