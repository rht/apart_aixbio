from __future__ import annotations

from typing import Any

import httpx


UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


async def fetch_uniprot_entry(uniprot_id: str) -> dict[str, Any]:
    url = f"{UNIPROT_BASE}/{uniprot_id}.json"
    async with httpx.AsyncClient(timeout=30) as client:
        resp = await client.get(url)
        resp.raise_for_status()
        return resp.json()


def extract_sequence(entry: dict[str, Any]) -> str:
    return entry["sequence"]["value"]


def extract_features(entry: dict[str, Any], feature_type: str) -> list[dict[str, Any]]:
    return [
        f for f in entry.get("features", [])
        if f.get("type") == feature_type
    ]


def extract_protein_name(entry: dict[str, Any]) -> str:
    protein = entry.get("proteinDescription", {})
    rec_name = protein.get("recommendedName")
    if rec_name:
        return rec_name.get("fullName", {}).get("value", "Unknown")
    sub_names = protein.get("submissionNames", [])
    if sub_names:
        return sub_names[0].get("fullName", {}).get("value", "Unknown")
    return "Unknown"
