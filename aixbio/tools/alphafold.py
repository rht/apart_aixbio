from __future__ import annotations

import logging

from aixbio.models.structure import StructureResult

logger = logging.getLogger(__name__)


async def predict_structure(
    chain_id: str,
    aa_sequence: str,
    reference_pdb: str | None = None,
) -> StructureResult:
    """Placeholder for AlphaFold 3 API integration.

    In production, this would call the AlphaFold Server API
    and optionally compute RMSD against a reference PDB.

    WARNING: This function is STUBBED and returns meaningless zero values.
    Do not use pLDDT, RMSD, or perplexity from this output for any decisions.
    """
    logger.warning(
        "AlphaFold integration is STUBBED — returning placeholder zeros for "
        f"chain '{chain_id}'. pLDDT, RMSD, and perplexity values are meaningless. "
        "Integrate the AlphaFold Server API for real structural predictions."
    )
    return StructureResult(
        id=chain_id,
        plddt_mean=0.0,
        rmsd_to_ref=None,
        perplexity=None,
        structure_file="",
    )
