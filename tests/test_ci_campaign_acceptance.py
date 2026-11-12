"""A portable CI assertion must never waive a scientific or provenance failure."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from conftest import check_runtime_aware_campaign
from gaussian_dynamics.runtime_compatibility_v233 import (
    assess_runtime_compatibility_v233, release_locked_runtime_profile_v233,
)


@pytest.fixture
def canonical_campaign():
    path = Path(__file__).resolve().parents[1] / "results/v0233_transport_compatibility_campaign.json"
    return json.loads(path.read_text())


def test_archived_canonical_identity_is_still_required_to_pass(canonical_campaign):
    report = check_runtime_aware_campaign(canonical_campaign)
    assert report["canonical_release_accepted"] is True
    assert report["failed_release_gates"] == []


def test_noncanonical_binary_must_report_release_rejection(canonical_campaign):
    result = deepcopy(canonical_campaign)
    runtime = result["inherited_v0232"]["pyscf_runtime_evidence"]["runtime"]
    runtime["python_executable_sha256"] = "0"*64
    locked = assess_runtime_compatibility_v233(runtime, release_locked_runtime_profile_v233())
    result["runtime_profile_controls"]["release_locked_report"] = locked.as_dict()
    gate = "runtime_profile::release_locked_profile_passes_canonical_runtime"
    result["acceptance"]["checks"][gate] = False
    result["acceptance"]["passed"] = False
    report = check_runtime_aware_campaign(result)
    assert report["canonical_release_accepted"] is False
    assert report["failed_release_gates"] == [gate]
    result["acceptance"]["passed"] = True
    with pytest.raises(AssertionError):
        check_runtime_aware_campaign(result)


def test_any_scientific_failure_remains_a_test_failure(canonical_campaign):
    result = deepcopy(canonical_campaign)
    gate = next(key for key in result["acceptance"]["checks"] if not key.startswith("runtime_profile::"))
    result["acceptance"]["checks"][gate] = False
    result["acceptance"]["passed"] = False
    with pytest.raises(AssertionError):
        check_runtime_aware_campaign(result)


def test_forged_identity_report_is_not_trusted(canonical_campaign):
    result = deepcopy(canonical_campaign)
    result["runtime_profile_controls"]["release_locked_report"]["compatible"] = False
    with pytest.raises(AssertionError):
        check_runtime_aware_campaign(result)
