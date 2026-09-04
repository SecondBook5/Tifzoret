from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "07_networks" / "networks.py"
SPEC = importlib.util.spec_from_file_location("tifzoret_networks", SCRIPT)
NETWORKS = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(NETWORKS)


class FakeResponse:
    def __init__(self, payload: bytes):
        self.payload = payload
        self.headers = {"X-STRING-Version": "fixture-release"}

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def read(self):
        return self.payload


def test_string_provider_caches_payload_and_receipt(tmp_path, monkeypatch):
    payload = b"queryIndex\tqueryItem\tstringId\tpreferredName\n0\tActa2\t10090.ENSMUSP1\tActa2\n"
    calls = []

    def fake_urlopen(request, timeout):
        calls.append((request, timeout))
        return FakeResponse(payload)

    monkeypatch.setattr(NETWORKS.urllib.request, "urlopen", fake_urlopen)
    parameters = {"identifiers": "Acta2", "species": 10090, "limit": 1}
    first = NETWORKS.cached_post("get_string_ids", parameters, tmp_path, False, False)
    second = NETWORKS.cached_post("get_string_ids", parameters, tmp_path, True, False)
    assert first == second == payload.decode()
    assert len(calls) == 1
    receipt = json.loads(next(tmp_path.glob("string_get_string_ids_*.json")).read_text())
    assert receipt["database_release"] == "fixture-release"
    assert receipt["requested_identifier_count"] == 1


def test_large_network_requests_cover_all_chunk_pairs_without_pruning(tmp_path, monkeypatch):
    calls = []

    def fake_cached_post(endpoint, parameters, *_args):
        submitted = parameters["identifiers"].split("\r")
        calls.append(submitted)
        return (
            "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tscore\n"
            f"{submitted[0]}\t{submitted[-1]}\t{submitted[0]}\t{submitted[-1]}\t0.9\n"
        )

    monkeypatch.setattr(NETWORKS, "cached_post", fake_cached_post)
    identifiers = [f"protein_{index}" for index in range(7)]
    rows, api_calls = NETWORKS.fetch_induced_network(
        identifiers, 10090, 700, tmp_path, False, False, batch_size=3
    )

    assert api_calls == 6  # three chunks: all i <= j combinations
    assert all(len(call) <= 6 for call in calls)
    assert set().union(*(set(call) for call in calls)) == set(identifiers)
    assert rows


@pytest.mark.skipif(os.environ.get("TIFZORET_LIVE") != "1", reason="scheduled live-provider test")
def test_live_string_mouse_mapping(tmp_path):
    response = NETWORKS.cached_post(
        "get_string_ids",
        {"identifiers": "Acta2", "species": 10090, "limit": 1, "echo_query": 1},
        tmp_path,
        False,
        True,
    )
    rows = NETWORKS.parse_response(response)
    assert rows
    assert rows[0]["preferredName"].lower() == "acta2"


def test_cached_post_retries_on_transient_failures(tmp_path, monkeypatch):
    """Verify that cached_post retries transient network failures and succeeds."""
    payload = b"queryIndex\tqueryItem\tstringId\tpreferredName\n0\tActa2\t10090.ENSMUSP1\tActa2\n"
    attempts = []

    def fake_urlopen(request, timeout):
        attempts.append(len(attempts) + 1)
        if len(attempts) < 3:
            # First two attempts fail with URLError
            raise NETWORKS.urllib.error.URLError("Temporary network error")
        # Third attempt succeeds
        return FakeResponse(payload)

    monkeypatch.setattr(NETWORKS.urllib.request, "urlopen", fake_urlopen)
    parameters = {"identifiers": "Acta2", "species": 10090, "limit": 1}
    result = NETWORKS.cached_post("get_string_ids", parameters, tmp_path, False, False)
    assert result == payload.decode()
    assert len(attempts) == 3  # Failed twice, succeeded on third


def test_cached_post_fails_with_provider_named_error_after_max_attempts(tmp_path, monkeypatch):
    """Verify that cached_post raises a STRING-named error after exhausting retries."""
    def fake_urlopen(request, timeout):
        raise NETWORKS.urllib.error.URLError("Persistent network error")

    monkeypatch.setattr(NETWORKS.urllib.request, "urlopen", fake_urlopen)
    parameters = {"identifiers": "Acta2", "species": 10090, "limit": 1}

    with pytest.raises(RuntimeError) as exc_info:
        NETWORKS.cached_post("get_string_ids", parameters, tmp_path, False, False)

    error_message = str(exc_info.value)
    assert "STRING" in error_message
    assert "3 attempts" in error_message
    assert "live network call" in error_message
    assert "offline mode" in error_message


def test_cached_post_success_path_unchanged(tmp_path, monkeypatch):
    """Verify that successful calls return identical results (behavior preservation)."""
    payload = b"stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tscore\n10090.P1\t10090.P2\tGene1\tGene2\t0.95\n"
    calls = []

    def fake_urlopen(request, timeout):
        calls.append(timeout)
        return FakeResponse(payload)

    monkeypatch.setattr(NETWORKS.urllib.request, "urlopen", fake_urlopen)
    parameters = {"identifiers": "Gene1\rGene2", "species": 10090}
    result = NETWORKS.cached_post("network", parameters, tmp_path, False, False)

    # Verify the success path is byte-identical
    assert result == payload.decode()
    assert len(calls) == 1
    assert calls[0] == 180  # timeout preserved

    # Verify receipt is written correctly
    receipt = json.loads(next(tmp_path.glob("string_network_*.json")).read_text())
    assert receipt["provider"] == "STRING"
    assert receipt["requested_identifier_count"] == 2
