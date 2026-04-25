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
