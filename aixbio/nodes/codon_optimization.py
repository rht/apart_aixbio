from __future__ import annotations

from aixbio.models.dna import DNAChain
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.cai import compute_cai
from aixbio.tools.codon_tables import best_ecoli_codon, split_codons, synonymous_alternatives, translate_codon
from aixbio.tools.gc import compute_gc
from aixbio.tools.restriction_sites import ENZYME_SITES, find_restriction_sites


def codon_optimization(state: ChainSubgraphState) -> dict:
    chain = state["chain"]
    avoid = state["avoid_sites"]

    codons = [best_ecoli_codon(aa) for aa in chain.aa_sequence]
    dna = "".join(codons)

    for _pass in range(5):
        hits = find_restriction_sites(dna, avoid)
        if not hits:
            break
        for enzyme, pos in hits:
            site_len = len(ENZYME_SITES[enzyme])
            codon_start = pos // 3
            codon_end = (pos + site_len - 1) // 3 + 1
            for ci in range(codon_start, min(codon_end, len(codons))):
                alts = synonymous_alternatives(codons[ci])
                if alts:
                    codons[ci] = alts[0]
                    break
        dna = "".join(codons)

    cai_score = compute_cai(dna)
    gc_content = compute_gc(dna)

    dna_chain = DNAChain(
        id=chain.id,
        dna_sequence=dna,
        cai_score=cai_score,
        gc_content=gc_content,
    )
    return {"optimized_dna": dna_chain}
