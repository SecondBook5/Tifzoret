from __future__ import annotations
import re
from pathlib import Path
from tools.parity.validate_register import load_and_validate

SEVERITY_RANK = {"S1": 0, "S2": 1, "S3": 2, "S4": 3, "S5": 4}

EXPECTED_STAGES = [
    "qc", "de", "pathways", "ontology", "composition", "regulators",
    "network", "hypothesis", "sva", "wgcna", "mediation", "power", "multilayer"
]

EXPECTED_PANELS = [
    "Set1-A", "Set1-B", "Set1-C", "Set1-D", "Set1-E",
    "Set2-A", "Set2-B", "Set2-C", "Set2-D1", "Set2-D2", "Set2-D3",
    "Set2-E", "Set2-F", "Set2-G", "Set2-H", "Set2-I"
]


def _normalize_for_panel(s: str) -> str:
    """Strip non-alphanumerics for panel matching."""
    return re.sub(r'[^a-z0-9]', '', s.lower())


def _stem(token: str) -> str:
    """Collapse a single trailing plural 's' so 'networks' stems to 'network'."""
    return token[:-1] if len(token) > 1 and token.endswith("s") else token


def _area_tokens(area: str) -> set[str]:
    """
    Split an area label into stemmed alphanumeric tokens for EXACT matching.

    Stage coverage is decided by whole-token equality, not substring
    containment. Substring matching produced silent false passes: 'sva' was
    satisfied by any 'gsva' finding, and the 2-letter 'de' was satisfied by
    'node', 'dendrograms', 'config-defaults', 'seed-determinism', etc. Under
    token matching those never match (gsva != sva; dendrogram != de), so the
    gate flips a genuinely-missing stage to missing. The plural stem keeps the
    real labels matching ('networks' -> network, 'pathways' -> pathway,
    'regulators' -> regulator).
    """
    return {_stem(tok) for tok in re.split(r'[^a-z0-9]+', area.lower()) if tok}


def merge_fragments(fragment_dir: str) -> list[dict]:
    """
    Merge all .yaml fragments in fragment_dir into a single ranked register.

    Returns entries sorted by (severity_rank, track, id).
    Raises ValueError if any fragment is invalid or ids collide across fragments.
    """
    fragment_path = Path(fragment_dir)
    all_entries = []
    seen_ids = {}

    # Load all .yaml files
    for yaml_file in sorted(fragment_path.glob("*.yaml")):
        entries = load_and_validate(str(yaml_file))

        # Check for duplicate IDs across fragments
        for entry in entries:
            entry_id = entry["id"]
            if entry_id in seen_ids:
                raise ValueError(
                    f"duplicate id {entry_id!r} found in {yaml_file.name} "
                    f"and {seen_ids[entry_id]}"
                )
            seen_ids[entry_id] = yaml_file.name
            all_entries.append(entry)

    # Sort by (severity_rank, track, id)
    all_entries.sort(key=lambda e: (
        SEVERITY_RANK[e["severity"]],
        e["track"],
        e["id"]
    ))

    return all_entries


def coverage_report(fragment_dir: str) -> dict:
    """
    Check coverage of stages and panels in fragments.

    Returns dict with:
    - stages_seen: set of covered stage tokens
    - panels_seen: set of covered panel ids
    - stages_missing: list of uncovered stages
    - panels_missing: list of uncovered panels
    """
    entries = merge_fragments(fragment_dir)

    stages_seen = set()
    panels_seen = set()

    for entry in entries:
        area = entry["area"]
        area_normalized = _normalize_for_panel(area)
        tokens = _area_tokens(area)

        # Check stages: a stage is covered when its stemmed name equals a
        # stemmed whole token of the area (exact, not substring — see
        # _area_tokens for why substring matching produced false passes).
        for stage in EXPECTED_STAGES:
            if _stem(stage) in tokens:
                stages_seen.add(stage)

        # Check panels: normalized area starts with normalized panel prefix.
        # Known limitation: the auditor labels the three DE-clustering variants
        # collectively as "Set2-D" (never Set2-D1/D2/D3 individually), so any
        # "set2d" area satisfies all three ids at once. Per-variant coverage is
        # therefore NOT independently detectable here; it is asserted instead by
        # the B2 findings that audit all three options (global / program-grouped
        # / direct-labels). Treat the D-trio as one panel for gate purposes.
        if area_normalized.startswith("set2d"):
            panels_seen.update(["Set2-D1", "Set2-D2", "Set2-D3"])

        for panel in EXPECTED_PANELS:
            panel_normalized = _normalize_for_panel(panel)
            if area_normalized.startswith(panel_normalized):
                panels_seen.add(panel)

    stages_missing = [s for s in EXPECTED_STAGES if s not in stages_seen]
    panels_missing = [p for p in EXPECTED_PANELS if p not in panels_seen]

    return {
        "stages_seen": stages_seen,
        "panels_seen": panels_seen,
        "stages_missing": stages_missing,
        "panels_missing": panels_missing
    }
