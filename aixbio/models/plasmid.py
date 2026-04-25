from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PlasmidChain:
    id: str
    genbank_file: str
    vector: str
    insert_size: int


@dataclass(frozen=True)
class PlasmidRecord:
    chains: tuple[PlasmidChain, ...]
