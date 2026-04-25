from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DNAChain:
    id: str
    dna_sequence: str
    cai_score: float
    gc_content: float


@dataclass(frozen=True)
class OptimizedDNA:
    chains: tuple[DNAChain, ...]


@dataclass(frozen=True)
class CassetteElement:
    start: str
    tag: str
    protease: str
    gene: str
    stop: str


@dataclass(frozen=True)
class CassetteChain:
    id: str
    full_dna: str
    elements: CassetteElement


@dataclass(frozen=True)
class CassetteDNA:
    chains: tuple[CassetteChain, ...]
