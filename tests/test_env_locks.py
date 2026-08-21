"""Every pinned conda dependency in an environment spec must have a matching
line in its frozen lock. This guards the invariant the locks/README states:
"if a YAML pin and its lock disagree, the YAML was changed without refreshing
the lock" — exactly the loose end that opt-in module deps introduced.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml

ENVS = Path(__file__).resolve().parents[1] / "src" / "tifzoret" / "workflow" / "envs"
SPECS = {"core.yaml": "core.lock.txt", "network.yaml": "network.lock.txt", "r.yaml": "r.lock.txt"}


def _conda_deps(spec: Path) -> list[str]:
    """Conda (string) dependencies of an environment YAML, skipping any nested
    pip block, which the locks do not carry."""
    document = yaml.safe_load(spec.read_text())
    return [dep for dep in document.get("dependencies", []) if isinstance(dep, str)]


@pytest.mark.parametrize("spec_name,lock_name", sorted(SPECS.items()))
def test_every_pinned_dependency_is_in_the_lock(spec_name, lock_name):
    spec = ENVS / spec_name
    lock_text = (ENVS / "locks" / lock_name).read_text()
    missing = []
    for dep in _conda_deps(spec):
        name, sep, rest = dep.partition("=")
        if not sep:
            missing.append(f"{dep} (unpinned — every direct dependency must carry a version)")
            continue
        version = rest.partition("=")[0]  # tolerate an optional name=version=build spec
        # Explicit-lock URLs end in `.../<name>-<resolved>-<build>.conda`; a spec
        # may pin a version prefix (e.g. `python=3.12` -> resolved `3.12.11`), so
        # the resolved version follows the pin with either `.` or the build `-`.
        if f"/{name}-{version}-" not in lock_text and f"/{name}-{version}." not in lock_text:
            missing.append(f"{name}={version} (no matching line in {lock_name})")
    assert not missing, f"{spec_name} pins out of sync with {lock_name}:\n  " + "\n  ".join(missing)


def _conda_names(spec: Path) -> set[str]:
    """Conda dependency names (pins/build strings stripped) of an env YAML."""
    return {dep.partition("=")[0] for dep in _conda_deps(spec)}


def test_root_environment_is_a_superset_of_the_per_rule_envs():
    """The root environment.yaml is the single-env `--no-conda` path: it must
    carry every dependency of the three per-rule envs, or a full run under
    `tifzoret run ... --no-conda` fails partway (the class of gap that left
    viper/ComplexHeatmap/scipy out and broke the all-modules path)."""
    root = ENVS.parents[3] / "environment.yaml"
    have = _conda_names(root)
    missing = {}
    for spec_name in SPECS:  # core.yaml, network.yaml, r.yaml
        need = _conda_names(ENVS / spec_name)
        absent = sorted(need - have)
        if absent:
            missing[spec_name] = absent
    assert not missing, (
        "environment.yaml is missing per-rule dependencies (the --no-conda union "
        "env must be a superset):\n  "
        + "\n  ".join(f"{k}: {', '.join(v)}" for k, v in missing.items())
    )
