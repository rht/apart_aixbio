from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Evo2Result:
    id: str
    log_prob: float
    mean_log_prob: float
    sequence_length: int
    method: str = "evo2"


@dataclass(frozen=True)
class Evo2Report:
    chains: tuple[Evo2Result, ...]
