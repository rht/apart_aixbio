"""Tests for aixbio.tools.evo2 — mocked biolmai SDK."""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture(autouse=True)
def _set_token(monkeypatch):
    monkeypatch.setenv("BIOLMAI_TOKEN", "test-token")


def _mock_predict(return_value):
    mock_model = MagicMock()
    mock_model.predict.return_value = return_value
    mock_cls = MagicMock(return_value=mock_model)
    return patch("biolmai.Model", mock_cls)


class TestScoreDnaSuccess:
    def test_dict_response(self):
        with _mock_predict({"log_prob": -120.5}):
            from aixbio.tools.evo2 import score_dna
            result = score_dna("A", "ATGCGTACC")

        assert result.method == "evo2"
        assert result.log_prob == -120.5
        assert result.mean_log_prob == round(-120.5 / 9, 6)
        assert result.sequence_length == 9

    def test_list_response(self):
        with _mock_predict([{"log_prob": -80.0}]):
            from aixbio.tools.evo2 import score_dna
            result = score_dna("B", "ATGCGT")

        assert result.method == "evo2"
        assert result.log_prob == -80.0


class TestScoreDnaErrors:
    def test_dict_error(self):
        with _mock_predict({"error": "model not found"}):
            from aixbio.tools.evo2 import score_dna
            result = score_dna("A", "ATGCGT")

        assert result.method == "evo2_failed"
        assert result.log_prob is None
        assert result.mean_log_prob is None

    def test_list_error(self):
        with _mock_predict([{"error": "rate limited"}]):
            from aixbio.tools.evo2 import score_dna
            result = score_dna("A", "ATGCGT")

        assert result.method == "evo2_failed"
        assert result.log_prob is None

    def test_missing_log_prob_key(self):
        with _mock_predict([{"score": -50.0}]):
            from aixbio.tools.evo2 import score_dna
            result = score_dna("A", "ATGCGT")

        assert result.method == "evo2_failed"
        assert result.log_prob is None


class TestScoreDnaTruncation:
    def test_truncation_sets_method(self):
        long_seq = "ATGC" * 2000  # 8000 bp
        with _mock_predict([{"log_prob": -500.0}]):
            from aixbio.tools.evo2 import score_dna
            result = score_dna("A", long_seq)

        assert result.method == "evo2_truncated"
        assert result.sequence_length == 4096
        assert result.log_prob == -500.0


class TestScoreDnaMissingToken:
    def test_raises_without_token(self, monkeypatch):
        monkeypatch.delenv("BIOLMAI_TOKEN", raising=False)
        monkeypatch.delenv("BIOLM_TOKEN", raising=False)

        from aixbio.tools.evo2 import score_dna
        with pytest.raises(RuntimeError, match="BIOLMAI_TOKEN"):
            score_dna("A", "ATGCGT")
