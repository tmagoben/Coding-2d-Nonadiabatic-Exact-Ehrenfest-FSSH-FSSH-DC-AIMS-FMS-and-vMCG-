"""CI assertions distinguish scientific regression from sealed runtime identity."""
import json
import hashlib
from pathlib import Path
import tarfile

import pytest

from gaussian_dynamics.runtime_compatibility_v233 import (
    assess_runtime_compatibility_v233, release_locked_runtime_profile_v233,
    scientifically_compatible_runtime_profile_v233,
)


def pytest_sessionstart(session):
    """Hydrate missing archived fixtures without replacing local result files."""
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "results/v0302_historical_evidence_recovery.json").read_text())
    with tarfile.open(root / "tests/fixtures/v0300_historical_evidence.tar.xz") as archive:
        assert set(archive.getnames()) == {item["path"] for item in manifest["files"]}
        for item in manifest["files"]:
            relative = Path(item["path"])
            assert len(relative.parts) == 2 and relative.parts[0] == "results" and relative.suffix == ".json"
            member = archive.getmember(item["path"])
            assert member.isfile()
            data = archive.extractfile(member).read()
            assert hashlib.sha256(data).hexdigest() == item["sha256"], item["path"]
            target = root / relative
            if not target.exists():
                target.parent.mkdir(exist_ok=True)
                target.write_bytes(data)


def check_runtime_aware_campaign(result):
    """Require every scientific gate and independently verify the identity gate.

    The campaign's original acceptance result is never modified. A hosted Python
    executable need not match the sealed release binary; its identity gate must
    report rejection, and no other failed gate is allowed.
    """
    v233 = result
    while v233["release"] != "v0.23.3":
        parents = [value for key, value in v233.items() if key.startswith("inherited_v")]
        assert len(parents) == 1, "campaign ancestry must be unambiguous"
        v233 = parents[0]
    runtime = v233["inherited_v0232"]["pyscf_runtime_evidence"]["runtime"]
    locked = assess_runtime_compatibility_v233(runtime, release_locked_runtime_profile_v233())
    scientific = assess_runtime_compatibility_v233(runtime, scientifically_compatible_runtime_profile_v233())
    assert scientific.compatible, scientific.mismatches
    assert v233["runtime_profile_controls"]["release_locked_report"] == locked.as_dict()
    assert v233["runtime_profile_controls"]["scientifically_compatible_report"] == scientific.as_dict()
    acceptance = result["acceptance"]
    checks = acceptance["checks"]
    assert all(type(value) is bool for value in checks.values())
    identity = [key for key in checks if key.endswith("runtime_profile::release_locked_profile_passes_canonical_runtime")]
    assert len(identity) == 1
    assert checks[identity[0]] is locked.compatible
    failed = {key for key, passed in checks.items() if not passed}
    expected = set() if locked.compatible else {identity[0]}
    assert failed == expected, sorted(failed)
    assert acceptance["passed"] is (not failed)
    return {"release": result["release"], "scientific_regression_passed": True,
            "canonical_release_accepted": acceptance["passed"],
            "failed_release_gates": sorted(failed),
            "runtime_identity": locked.as_dict(), "scientific_compatibility": scientific.as_dict()}


@pytest.fixture
def assert_runtime_aware_campaign():
    def check(result):
        report = check_runtime_aware_campaign(result)
        path = Path(__file__).resolve().parents[1] / "results" / f"ci-campaign-{result['release']}.json"
        path.parent.mkdir(exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return check
