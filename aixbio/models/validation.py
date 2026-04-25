from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CheckResult:
    name: str
    passed: bool
    value: float | str
    threshold: str


@dataclass(frozen=True)
class ChainValidation:
    id: str
    passed: bool
    checks: tuple[CheckResult, ...]


@dataclass(frozen=True)
class ValidationReport:
    chains: tuple[ChainValidation, ...]
    all_passed: bool
