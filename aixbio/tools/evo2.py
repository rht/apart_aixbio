from __future__ import annotations

import logging
import os

from aixbio.models.structure import Evo2Result

logger = logging.getLogger(__name__)

_MODEL_NAME = "evo2-1b-base"
_MAX_SEQUENCE_LENGTH = 4096


def score_dna(chain_id: str, dna_sequence: str) -> Evo2Result:
    """Score a DNA sequence using Evo2 log-probability.

    Returns an Evo2Result with total and mean (per-nucleotide) log-prob.
    Sequences longer than 4096 bp are truncated to stay within the model's
    context window.
    """
    if not (os.environ.get("BIOLMAI_TOKEN") or os.environ.get("BIOLM_TOKEN")):
        raise RuntimeError(
            "Evo2 validation requires BIOLMAI_TOKEN in .env. "
            "Get a token from https://biolm.ai or run `biolmai login`."
        )

    from biolmai import Model

    seq = dna_sequence.upper()
    seq_len = len(seq)

    truncated = False
    if seq_len > _MAX_SEQUENCE_LENGTH:
        logger.warning(
            f"Chain {chain_id}: {seq_len} bp exceeds Evo2 context ({_MAX_SEQUENCE_LENGTH} bp), "
            f"truncating to first {_MAX_SEQUENCE_LENGTH} bp"
        )
        seq = seq[:_MAX_SEQUENCE_LENGTH]
        seq_len = _MAX_SEQUENCE_LENGTH
        truncated = True

    model = Model(_MODEL_NAME, progress=False)
    result = model.predict(items=[{"sequence": seq}])

    record = result[0] if isinstance(result, list) and result else result
    if isinstance(record, dict) and "error" in record:
        logger.error(f"Evo2 API error for chain {chain_id}: {record['error']}")
        return Evo2Result(
            id=chain_id,
            log_prob=None,
            mean_log_prob=None,
            sequence_length=seq_len,
            method="evo2_failed",
        )

    log_prob = record.get("log_prob") if isinstance(record, dict) else None
    if log_prob is None:
        logger.error(f"Evo2 returned no log_prob for chain {chain_id}: {record!r}")
        return Evo2Result(
            id=chain_id,
            log_prob=None,
            mean_log_prob=None,
            sequence_length=seq_len,
            method="evo2_failed",
        )

    mean_log_prob = log_prob / seq_len if seq_len > 0 else 0.0

    logger.info(
        f"Evo2: chain {chain_id} ({seq_len} bp), "
        f"log_prob={log_prob:.2f}, mean={mean_log_prob:.4f}/nt"
    )

    method = "evo2_truncated" if truncated else "evo2"
    return Evo2Result(
        id=chain_id,
        log_prob=log_prob,
        mean_log_prob=round(mean_log_prob, 6),
        sequence_length=seq_len,
        method=method,
    )
