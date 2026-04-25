from __future__ import annotations

import logging
from typing import Any

import httpx

logger = logging.getLogger(__name__)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

_MAX_RETRIES = 3
_BACKOFF_BASE = 1.0  # seconds


async def fetch_uniprot_entry(uniprot_id: str) -> dict[str, Any]:
    """Fetch a UniProt entry with retry and exponential backoff.

    Retries up to _MAX_RETRIES times on transient errors (timeouts,
    5xx responses, connection errors).
    """
    import asyncio

    last_error: Exception | None = None
    for attempt in range(_MAX_RETRIES):
        try:
            url = f"{UNIPROT_BASE}/{uniprot_id}.json"
            async with httpx.AsyncClient(timeout=30) as client:
                resp = await client.get(url)
                resp.raise_for_status()
                return resp.json()
        except (httpx.TimeoutException, httpx.ConnectError, httpx.HTTPStatusError) as exc:
            last_error = exc
            if isinstance(exc, httpx.HTTPStatusError) and exc.response.status_code < 500:
                # Client errors (4xx) should not be retried
                raise
            wait = _BACKOFF_BASE * (2 ** attempt)
            logger.warning(
                f"UniProt fetch attempt {attempt + 1}/{_MAX_RETRIES} failed for "
                f"{uniprot_id}: {exc}. Retrying in {wait:.1f}s..."
            )
            await asyncio.sleep(wait)

    raise RuntimeError(
        f"Failed to fetch UniProt entry {uniprot_id} after {_MAX_RETRIES} attempts"
    ) from last_error


def fetch_uniprot_entry_sync(uniprot_id: str) -> dict[str, Any]:
    """Synchronous wrapper for fetch_uniprot_entry.

    Safe to call from sync code regardless of whether an event loop is running.
    Uses httpx.Client directly to avoid asyncio.run() nesting issues.
    """
    last_error: Exception | None = None
    import time
    for attempt in range(_MAX_RETRIES):
        try:
            url = f"{UNIPROT_BASE}/{uniprot_id}.json"
            with httpx.Client(timeout=30) as client:
                resp = client.get(url)
                resp.raise_for_status()
                return resp.json()
        except (httpx.TimeoutException, httpx.ConnectError, httpx.HTTPStatusError) as exc:
            last_error = exc
            if isinstance(exc, httpx.HTTPStatusError) and exc.response.status_code < 500:
                raise
            wait = _BACKOFF_BASE * (2 ** attempt)
            logger.warning(
                f"UniProt fetch attempt {attempt + 1}/{_MAX_RETRIES} failed for "
                f"{uniprot_id}: {exc}. Retrying in {wait:.1f}s..."
            )
            time.sleep(wait)

    raise RuntimeError(
        f"Failed to fetch UniProt entry {uniprot_id} after {_MAX_RETRIES} attempts"
    ) from last_error


def extract_sequence(entry: dict[str, Any]) -> str:
    seq = entry.get("sequence")
    if seq is None or "value" not in seq:
        raise ValueError(
            f"UniProt entry has no sequence data. "
            f"Available keys: {sorted(entry.keys())}"
        )
    return seq["value"]


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


def extract_mature_chains(entry: dict[str, Any], compound_id: str) -> tuple[list[dict], str]:
    """Deterministically extract mature chain(s) from UniProt feature annotations.

    Returns (chains_list, reasoning_string).

    Priority:
    1. If "Chain" features exist, use them (multi-chain or annotated proteins)
    2. Otherwise, strip signal peptide and use the remaining sequence
    """
    import re

    full_sequence = extract_sequence(entry)
    protein_name = extract_protein_name(entry)

    chain_features = extract_features(entry, "Chain")

    if chain_features:
        chains = []
        for feat in chain_features:
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value", 1) - 1
            end = loc.get("end", {}).get("value", len(full_sequence))
            aa_seq = full_sequence[start:end]
            desc = feat.get("description", "")
            clean = re.sub(r"[^a-zA-Z0-9 ]", "", desc).strip().replace(" ", "_")
            chain_id = clean or f"chain_{len(chains) + 1}"
            chains.append({
                "id": chain_id,
                "aa_sequence": aa_seq,
                "length": len(aa_seq),
            })
        reasoning = (
            f"Extracted {len(chains)} chain(s) from UniProt 'Chain' feature annotations "
            f"for {protein_name} ({compound_id})."
        )
        return chains, reasoning

    signal_features = extract_features(entry, "Signal")
    transit_features = extract_features(entry, "Transit peptide")
    mature_start = 0
    for feat in signal_features + transit_features:
        loc = feat.get("location", {})
        end = loc.get("end", {}).get("value", 0)
        mature_start = max(mature_start, end)

    aa_seq = full_sequence[mature_start:]
    clean = re.sub(r"[^a-zA-Z0-9 ]", "", protein_name).strip().replace(" ", "_")
    chain_id = f"{clean}_{compound_id}" if clean else compound_id

    reasoning = (
        f"Single-chain protein {protein_name} ({compound_id}). "
        f"Stripped {mature_start} residue signal/transit peptide. "
        f"Mature sequence: {len(aa_seq)} aa."
    )
    if mature_start == 0:
        reasoning = (
            f"Single-chain protein {protein_name} ({compound_id}). "
            f"No signal peptide annotated. Using full sequence: {len(aa_seq)} aa."
        )

    return [{
        "id": chain_id,
        "aa_sequence": aa_seq,
        "length": len(aa_seq),
    }], reasoning
