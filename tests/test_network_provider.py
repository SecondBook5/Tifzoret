from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "bulk_rna_frame" / "workflow" / "scripts" / "networks.py"
SPEC = importlib.util.spec_from_file_location("bulk_rna_frame_networks", SCRIPT)
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
