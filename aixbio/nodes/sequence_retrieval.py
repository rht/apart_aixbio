from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

from langchain_core.messages import HumanMessage, SystemMessage

from aixbio.config import LLM_MAX_TOKENS, LLM_MODEL, OPENROUTER_API_KEY, OPENROUTER_BASE_URL, ChatOpenRouter
from aixbio.models.audit import AgentDecision
from aixbio.models.protein import Chain, ProteinRecord
from aixbio.prompts.sequence_retrieval import SEQUENCE_RETRIEVAL_SYSTEM
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.uniprot import extract_protein_name, extract_sequence, fetch_uniprot_entry_sync

logger = logging.getLogger(__name__)


def sequence_retrieval_agent(state: PipelineState) -> dict:
    compound_id = state["compound_id"]

    # Use sync wrapper to avoid asyncio.run() nesting issues (Issue #15)
    entry = fetch_uniprot_entry_sync(compound_id)
    full_sequence = extract_sequence(entry)
    protein_name = extract_protein_name(entry)

    features_summary = json.dumps(entry.get("features", []), indent=2, default=str)

    llm = ChatOpenRouter(
        model=LLM_MODEL,
        temperature=0,
        max_tokens=LLM_MAX_TOKENS,
        openai_api_key=OPENROUTER_API_KEY,
        openai_api_base=OPENROUTER_BASE_URL,
    )
    response = llm.invoke([
        SystemMessage(content=SEQUENCE_RETRIEVAL_SYSTEM),
        HumanMessage(content=(
            f"UniProt ID: {compound_id}\n"
            f"Protein name: {protein_name}\n"
            f"Full sequence ({len(full_sequence)} aa): {full_sequence}\n\n"
            f"Feature annotations:\n{features_summary}\n\n"
            "Extract the mature chain(s) for E. coli expression."
        )),
    ])

    response_text = response.content
    if isinstance(response_text, list):
        response_text = response_text[0].get("text", "") if response_text else ""

    try:
        json_start = response_text.index("{")
        json_end = response_text.rindex("}") + 1
        parsed = json.loads(response_text[json_start:json_end])
    except (ValueError, json.JSONDecodeError):
        logger.warning(
            f"FALLBACK ACTIVATED for {compound_id}: LLM response could not be parsed. "
            f"Using full precursor sequence ({len(full_sequence)} aa). This may include "
            f"signal peptides and pro-domains that should NOT be expressed."
        )
        parsed = _fallback_single_chain(compound_id, protein_name, full_sequence)

    chains = tuple(
        Chain(
            id=c["id"],
            aa_sequence=c["aa_sequence"],
            length=c["length"],
        )
        for c in parsed["chains"]
    )

    protein_record = ProteinRecord(
        uniprot_id=compound_id,
        name=protein_name,
        chains=chains,
    )

    reasoning = parsed.get("reasoning", "No reasoning provided")

    decision = AgentDecision(
        node="sequence_retrieval_agent",
        reasoning=reasoning,
        action=f"Extracted {len(chains)} chain(s) from {compound_id}",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=f"UniProt {compound_id} ({protein_name}, {len(full_sequence)} aa)",
        output_summary=", ".join(f"{c.id}: {c.length} aa" for c in chains),
    )

    return {
        "protein_record": protein_record,
        "chain_extraction_reasoning": reasoning,
        "decision_log": [decision],
    }


def _fallback_single_chain(
    compound_id: str, name: str, sequence: str
) -> dict:
    return {
        "chains": [{
            "id": f"{name}_{compound_id}",
            "aa_sequence": sequence,
            "length": len(sequence),
        }],
        "reasoning": (
            "FALLBACK: LLM response parsing failed. Using full precursor sequence. "
            "WARNING: This may include signal peptides and pro-domains. "
            "The resulting protein may not fold correctly in E. coli."
        ),
        "fallback_used": True,
    }
