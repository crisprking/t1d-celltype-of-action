"""Sanity tests for t1d_coa.config."""

import os
from pathlib import Path
from importlib import reload


def test_project_root_overridable(monkeypatch, tmp_path):
    monkeypatch.setenv("T1D_COA_ROOT", str(tmp_path))
    from t1d_coa import config
    reload(config)
    assert config.PROJECT_ROOT == tmp_path
    assert config.RAW == tmp_path / "data" / "raw"
    assert config.PROCESSED == tmp_path / "data" / "processed"


def test_ensure_dirs_creates_layout(monkeypatch, tmp_path):
    monkeypatch.setenv("T1D_COA_ROOT", str(tmp_path))
    from t1d_coa import config
    reload(config)
    config.ensure_dirs()
    for sub in ("data/raw", "data/interim", "data/processed",
                "results", "data/raw/cellxgene_expr",
                "data/raw/cellxgene_obs"):
        assert (tmp_path / sub).is_dir(), f"missing {sub}"


def test_constants_have_expected_types():
    from t1d_coa import config
    assert isinstance(config.HPAP_API, str)
    assert config.HPAP_API.startswith("https://")
    assert isinstance(config.TAU_THRESHOLD, float)
    assert 0 < config.TAU_THRESHOLD < 1
    assert config.EXPECTED_N_CELLS > 0
