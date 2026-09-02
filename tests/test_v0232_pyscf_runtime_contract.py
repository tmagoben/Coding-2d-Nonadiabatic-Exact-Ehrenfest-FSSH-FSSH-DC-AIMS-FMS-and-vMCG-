from types import SimpleNamespace

import pytest

import gaussian_dynamics.pyscf_runtime_v232 as runtime
from gaussian_dynamics.pyscf_nac_convention_v232 import (
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
    pyscf_state_tuple_for_internal_dij_v232,
    require_exact_pyscf_version_v232,
)


def test_empirical_mapping_is_distinct_and_unambiguous():
    assert PYSCF_NAC_EMPIRICAL_MAPPING_V232 != (
        PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232
    )
    assert pyscf_state_tuple_for_internal_dij_v232(0, 2) == (0, 2)
    assert pyscf_state_tuple_for_internal_dij_v232(2, 0) == (2, 0)
    with pytest.raises(ValueError, match="differ"):
        pyscf_state_tuple_for_internal_dij_v232(1, 1)

    assert require_exact_pyscf_version_v232(
        SimpleNamespace(__version__="2.13.1"),
        distribution_version="2.13.1",
    ) == "2.13.1"
    with pytest.raises(RuntimeError, match="requires exactly"):
        require_exact_pyscf_version_v232(
            SimpleNamespace(__version__="2.13.0"),
            distribution_version="2.13.1",
        )
    with pytest.raises(RuntimeError, match="requires exactly"):
        require_exact_pyscf_version_v232(
            SimpleNamespace(__version__="2.13.1"),
            distribution_version="2.13.0",
        )


def test_runtime_fails_closed_when_pyscf_is_absent(monkeypatch):
    monkeypatch.setattr(runtime.importlib.util, "find_spec", lambda name: None)

    probe = runtime.probe_pyscf_runtime_v232()

    assert not probe.installed
    assert not probe.usable
    with pytest.raises(ImportError, match="fails closed"):
        runtime.require_pyscf_runtime_v232()


def test_memory_probe_policy_rejects_unknown_modes_before_runtime_discovery():
    with pytest.raises(ValueError, match="memory_probe_policy"):
        with runtime.guarded_pyscf_runtime_v232(memory_probe_policy="unknown"):
            pass


def test_runtime_fails_closed_for_wrong_distribution_and_module(monkeypatch):
    fake_pyscf = SimpleNamespace(__version__="2.13.0")
    fake_nac = SimpleNamespace(NonAdiabaticCouplings=object)

    monkeypatch.setattr(runtime.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(
        runtime.metadata,
        "version",
        lambda name: "2.13.0" if name == "pyscf" else "test",
    )
    monkeypatch.setattr(
        runtime.importlib,
        "import_module",
        lambda name: fake_nac if name == "pyscf.nac.sacasscf" else fake_pyscf,
    )

    probe = runtime.probe_pyscf_runtime_v232()

    assert probe.installed
    assert probe.nac_api_available
    assert not probe.exact_version
    assert not probe.usable
    with pytest.raises(RuntimeError, match="version mismatch"):
        runtime.require_pyscf_runtime_v232()


def test_runtime_fails_closed_when_nac_api_is_missing(monkeypatch):
    fake_pyscf = SimpleNamespace(__version__="2.13.1")
    fake_nac = SimpleNamespace()

    monkeypatch.setattr(runtime.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(runtime.metadata, "version", lambda name: "2.13.1")
    monkeypatch.setattr(
        runtime.importlib,
        "import_module",
        lambda name: fake_nac if name == "pyscf.nac.sacasscf" else fake_pyscf,
    )

    probe = runtime.probe_pyscf_runtime_v232()

    assert probe.exact_version
    assert not probe.nac_api_available
    assert not probe.usable
    with pytest.raises(RuntimeError, match="NAC API"):
        runtime.require_pyscf_runtime_v232()
