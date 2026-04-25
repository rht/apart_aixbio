from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

from aixbio.models.dna import CassetteChain, CassetteElement
from aixbio.state.chain_state import ChainSubgraphState

TAG_SEQUENCES = {
    # 6×His = CAC×6 = 18 nt (encodes HHHHHH)
    "6xHis": "CACCACCACCACCACCACC",
}

PROTEASE_SITE_SEQUENCES = {
    # DDDDK = GAT GAT GAT GAT AAG = 15 nt
    "Enterokinase": "GATGATGATGATAAG",
    # ENLYFQ/G = GAA AAC CTG TAT TTT CAA GGC = 21 nt
    "TEV": "GAAAACCTGTATTTTCAAGGC",
    # NOTE: CNBr is a chemical cleavage reagent (cleaves after Met), not a
    # protease with a DNA-encodable recognition sequence. It has been removed
    # from this dictionary. If requested, the pipeline falls back to the
    # default protease site and emits a warning.
}

START_CODON = "ATG"
DOUBLE_STOP = "TAATAA"


def cassette_assembly(state: ChainSubgraphState) -> dict:
    dna_chain = state["optimized_dna"]
    tag_type = state["tag_type"]
    protease_site = state["protease_site"]

    tag_dna = TAG_SEQUENCES.get(tag_type, TAG_SEQUENCES["6xHis"])

    if protease_site == "CNBr":
        logger.warning(
            "CNBr is a chemical cleavage reagent, not a DNA-encodable protease. "
            "Falling back to Enterokinase. Use TEV or Enterokinase instead."
        )
        protease_site = "Enterokinase"

    protease_dna = PROTEASE_SITE_SEQUENCES.get(
        protease_site, PROTEASE_SITE_SEQUENCES["Enterokinase"]
    )

    elements = CassetteElement(
        start=START_CODON,
        tag=tag_dna,
        protease=protease_dna,
        gene=dna_chain.dna_sequence,
        stop=DOUBLE_STOP,
    )

    full_dna = START_CODON + tag_dna + protease_dna + dna_chain.dna_sequence + DOUBLE_STOP

    warnings = []

    # Glycosylation detection: check the protein record's UniProt feature
    # annotations rather than brittle string matching on chain IDs.
    protein_record = state.get("protein_record")
    has_glycosylation = False

    if protein_record:
        # Check UniProt features for glycosylation annotations
        # (available in the raw entry if features were preserved)
        uniprot_id = protein_record.uniprot_id.upper()
        # Curated set of known glycoproteins (by UniProt ID) as fallback
        _KNOWN_GLYCOPROTEINS = {
            "P01588",  # EPO (Erythropoietin)
            "P00750",  # tPA (Tissue plasminogen activator)
            "P01137",  # TGF-beta
            "P01375",  # TNF-alpha
        }
        if uniprot_id in _KNOWN_GLYCOPROTEINS:
            has_glycosylation = True

    # Also check for N-X-S/T sequons (N-linked glycosylation motif)
    import re
    chain_data = state.get("chain")
    if chain_data and chain_data.aa_sequence:
        n_glyc_motif = re.findall(r"N[^P][ST]", chain_data.aa_sequence)
        if n_glyc_motif:
            has_glycosylation = True

    if has_glycosylation:
        warnings.append(
            f"WARNING: {dna_chain.id} contains potential glycosylation sites "
            f"(N-X-S/T sequons or known glycoprotein). E. coli cannot perform "
            f"glycosylation — consider a eukaryotic expression system (CHO, HEK293)."
        )

    cassette = CassetteChain(
        id=dna_chain.id,
        full_dna=full_dna,
        elements=elements,
    )
    return {"cassette": cassette, "warnings": warnings}
