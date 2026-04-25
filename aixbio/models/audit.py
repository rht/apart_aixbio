from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AgentDecision:
    node: str
    reasoning: str
    action: str
    timestamp: str
    input_summary: str
    output_summary: str
