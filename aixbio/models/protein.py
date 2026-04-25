from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Chain:
    id: str
    aa_sequence: str
    length: int


@dataclass(frozen=True)
class ProteinRecord:
    uniprot_id: str
    name: str
    chains: tuple[Chain, ...]
