from __future__ import annotations

from aixbio.models.dna import CassetteChain, CassetteElement
from aixbio.state.chain_state import ChainSubgraphState

TAG_SEQUENCES = {
    "6xHis": "CACCACCACCACCACCAC",
}

PROTEASE_SITE_SEQUENCES = {
    "Enterokinase": "GATGATGATGATAAAG",
    "TEV": "GAAAACCTGTATTTTCAAGGC",
    "CNBr": "ATG",
}

START_CODON = "ATG"
DOUBLE_STOP = "TAATAA"


def cassette_assembly(state: ChainSubgraphState) -> dict:
    dna_chain = state["optimized_dna"]
    tag_type = state["tag_type"]
    protease_site = state["protease_site"]

    tag_dna = TAG_SEQUENCES.get(tag_type, TAG_SEQUENCES["6xHis"])
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
    glycosylation_compounds = {"EPO", "tPA"}
    chain_id_upper = dna_chain.id.upper()
    for compound in glycosylation_compounds:
        if compound.upper() in chain_id_upper:
            warnings.append(
                f"WARNING: {dna_chain.id} may require glycosylation which E. coli cannot perform."
            )

    cassette = CassetteChain(
        id=dna_chain.id,
        full_dna=full_dna,
        elements=elements,
    )
    return {"cassette": cassette, "warnings": warnings}
