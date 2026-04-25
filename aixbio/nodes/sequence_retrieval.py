from __future__ import annotations

import logging
from datetime import datetime, timezone

from aixbio.models.audit import AgentDecision
from aixbio.models.protein import Chain, ProteinRecord
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.uniprot import (
    extract_mature_chains,
    extract_protein_name,
    fetch_uniprot_entry_sync,
)

logger = logging.getLogger(__name__)


def sequence_retrieval_agent(state: PipelineState) -> dict:
    compound_id = state["compound_id"]

    entry = fetch_uniprot_entry_sync(compound_id)
    protein_name = extract_protein_name(entry)

    chains_data, reasoning = extract_mature_chains(entry, compound_id)

    chains = tuple(
        Chain(
            id=c["id"],
            aa_sequence=c["aa_sequence"],
            length=c["length"],
        )
        for c in chains_data
    )

    protein_record = ProteinRecord(
        uniprot_id=compound_id,
        name=protein_name,
        chains=chains,
    )

    decision = AgentDecision(
        node="sequence_retrieval_agent",
        reasoning=reasoning,
        action=f"Extracted {len(chains)} chain(s) from {compound_id}",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=f"UniProt {compound_id} ({protein_name})",
        output_summary=", ".join(f"{c.id}: {c.length} aa" for c in chains),
    )

    return {
        "protein_record": protein_record,
        "chain_extraction_reasoning": reasoning,
        "decision_log": [decision],
    }
