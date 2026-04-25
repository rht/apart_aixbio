from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RemediationAction:
    check_name: str
    fix_type: str
    positions_affected: tuple[int, ...]
    codons_before: tuple[str, ...]
    codons_after: tuple[str, ...]
    reasoning: str


@dataclass(frozen=True)
class PlannedFix:
    check_name: str
    strategy: str
    target_positions: tuple[int, ...]
    replacement_codons: tuple[str, ...]


@dataclass(frozen=True)
class RemediationPlan:
    actions: tuple[PlannedFix, ...]
    reasoning: str
    priority_order: tuple[str, ...]
