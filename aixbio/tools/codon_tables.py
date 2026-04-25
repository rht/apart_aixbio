from __future__ import annotations

import python_codon_tables as pct

_ECOLI_TABLE: dict[str, dict[str, float]] | None = None


def get_ecoli_table() -> dict[str, dict[str, float]]:
    global _ECOLI_TABLE
    if _ECOLI_TABLE is None:
        _ECOLI_TABLE = pct.get_codons_table("e_coli")
    return _ECOLI_TABLE


CODON_TO_AA: dict[str, str] = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)

RARE_CODONS_ECOLI = frozenset({"AGG", "AGA", "CGA", "CUA", "AUA", "CCC"})


def translate_codon(codon: str) -> str:
    return CODON_TO_AA[codon.upper()]


def translate_dna(dna: str) -> str:
    codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
    return "".join(translate_codon(c) for c in codons if len(c) == 3)


def best_ecoli_codon(aa: str) -> str:
    table = get_ecoli_table()
    codons = table.get(aa)
    if codons is None:
        raise ValueError(f"Unknown amino acid: {aa}")
    return max(codons, key=lambda c: codons[c])


def split_codons(dna: str) -> list[str]:
    return [dna[i:i+3] for i in range(0, len(dna), 3)]


def synonymous_alternatives(codon: str) -> list[str]:
    aa = translate_codon(codon)
    return [c for c in AA_TO_CODONS[aa] if c != codon.upper()]
