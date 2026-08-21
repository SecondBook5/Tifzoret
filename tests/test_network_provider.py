from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "networks.py"
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


@pytest.mark.skipif(os.environ.get("BULK_RNAFRAME_LIVE") != "1", reason="scheduled live-provider test")
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
