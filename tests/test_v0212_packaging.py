from pathlib import Path


def test_v0212_package_discovery_is_explicit():
    """Keep non-package release data from confusing setuptools discovery."""
    pyproject = (Path(__file__).parents[1] / "pyproject.toml").read_text()
    expected = '[tool.setuptools]\npackages = ["gaussian_dynamics"]'

    assert expected in pyproject
