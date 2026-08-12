from __future__ import annotations
import yaml

REQUIRED_KEYS = {"id","track","area","verdict","severity","summary",
                 "detail","reference","frame","evidence","target","status","commit"}
VERDICTS = {"REGRESSION","FIDELITY-GAP","BEST-PRACTICE","OK"}
SEVERITIES = {"S1","S2","S3","S4","S5"}
TARGETS = {"harness","figures","stats","inline"}
EVIDENCE_MODES = {"run","code-only"}
STATUSES = {"open","fixed-inline"}
TRACKS = {"A","B","C"}

def _check_ref_list(name, value, out, idx):
    if not isinstance(value, list) or not value:
        out.append(f"entry {idx}: {name} must be a non-empty list")
        return
    for j, ref in enumerate(value):
        if not isinstance(ref, dict) or "file" not in ref or "line" not in ref:
            out.append(f"entry {idx}: {name}[{j}] needs file and line")

def validate_fragment(data):
    errors = []
    seen = set()
    if not isinstance(data, list):
        return ["fragment must be a list of entries"]
    for i, e in enumerate(data):
        if not isinstance(e, dict):
            errors.append(f"entry {i}: not a mapping"); continue
        missing = REQUIRED_KEYS - set(e)
        for m in sorted(missing):
            errors.append(f"entry {i}: missing key '{m}'")
        if missing:
            continue
        if e["track"] not in TRACKS: errors.append(f"entry {i}: bad track {e['track']!r}")
        if e["verdict"] not in VERDICTS: errors.append(f"entry {i}: bad verdict {e['verdict']!r}")
        if e["severity"] not in SEVERITIES: errors.append(f"entry {i}: bad severity {e['severity']!r}")
        if e["target"] not in TARGETS: errors.append(f"entry {i}: bad target {e['target']!r}")
        if e["evidence"] not in EVIDENCE_MODES: errors.append(f"entry {i}: bad evidence {e['evidence']!r}")
        if e["status"] not in STATUSES: errors.append(f"entry {i}: bad status {e['status']!r}")
        if e["status"] == "fixed-inline" and not e.get("commit"):
            errors.append(f"entry {i}: status fixed-inline requires a commit sha")
        _check_ref_list("reference", e.get("reference"), errors, i)
        _check_ref_list("frame", e.get("frame"), errors, i)
        if e["id"] in seen:
            errors.append(f"entry {i}: duplicate id {e['id']!r}")
        seen.add(e["id"])
    return errors

def load_and_validate(path):
    data = yaml.safe_load(open(path, encoding="utf-8")) or []
    errs = validate_fragment(data)
    if errs:
        raise ValueError("invalid parity register:\n" + "\n".join(errs))
    return data
