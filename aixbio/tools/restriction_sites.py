from __future__ import annotations

import re

ENZYME_SITES: dict[str, str] = {
    "BamHI": "GGATCC",
    "XhoI": "CTCGAG",
    "EcoRI": "GAATTC",
    "HindIII": "AAGCTT",
    "NdeI": "CATATG",
    "NcoI": "CCATGG",
    "NotI": "GCGGCCGC",
    "SalI": "GTCGAC",
}


def find_restriction_sites(
    dna: str, enzymes: tuple[str, ...] | list[str]
) -> list[tuple[str, int]]:
    dna = dna.upper()
    hits: list[tuple[str, int]] = []
    for enzyme in enzymes:
        site = ENZYME_SITES.get(enzyme)
        if site is None:
            continue
        for m in re.finditer(site, dna):
            hits.append((enzyme, m.start()))
    return hits


def has_restriction_sites(
    dna: str, enzymes: tuple[str, ...] | list[str]
) -> bool:
    return len(find_restriction_sites(dna, enzymes)) > 0
