"""
Tests of the YACHT version recorded in database configs: `yacht train` writes `yacht_version`, and a database whose
`min_yacht_version` is newer than the installed YACHT is refused with a clear message (by `yacht run`, the index
tool and the pretrained-database download), while databases without it (YACHT < 2.0, or made by `yacht train`)
are used as before. Also: database files missing at the paths in the config are looked for next to the config.
"""
import argparse
import json
import os
import shutil

import pytest

from yacht import make_training_data_from_sketches, run_YACHT, utils, yacht_index

TESTDATA = os.path.join(os.path.dirname(__file__), "testdata")
REF_ZIP = os.path.join(TESTDATA, "20_genomes_sketches.zip")
SAMPLE_ZIP = os.path.join(TESTDATA, "sample.sig.zip")


@pytest.fixture(scope="module")
def trained_db(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("db")
    args = argparse.Namespace(
        ref_file=REF_ZIP, ksize=31, num_threads=2, ani_thresh=0.95, prefix="db", outdir=str(outdir), force=True
    )
    make_training_data_from_sketches.main(args)
    return os.path.join(str(outdir), "db_config.json")


def _copy_db(trained_db, tmp_path, **config_entries):
    """Copy of the trained database (in tmp_path/db) with extra config entries; returns its config path."""
    src = os.path.dirname(trained_db)
    dst = tmp_path / "db"
    shutil.copytree(src, dst)
    config_path = dst / "db_config.json"
    with open(config_path) as f:
        config = json.load(f)
    for key in ("manifest_file_path", "intermediate_files_dir", "index_file_path"):
        config[key] = str(dst / os.path.basename(config[key]))
    config.update(config_entries)
    with open(config_path, "w") as f:
        json.dump(config, f)
    return str(config_path)


def _run(config, tmp_path):
    return run_YACHT.main(argparse.Namespace(
        json=config, sample_file=SAMPLE_ZIP, significance=0.99, num_threads=None, keep_raw=False, show_all=False,
        min_coverage_list=[1, 0.5], outdir=str(tmp_path),
    ))


@pytest.mark.parametrize("config, version, ok", [
    ({}, "2.0.0", True),
    ({"min_yacht_version": "2.0.0"}, "2.0.0", True),
    ({"min_yacht_version": "2.0.0"}, "2.10.1", True),
    ({"min_yacht_version": "2.0.0"}, "1.3.2", False),
    ({"min_yacht_version": "2.1.0"}, "2.0.0", False),
    ({"min_yacht_version": "10.0"}, "9.9.9", False),
])
def test_database_version_problem(config, version, ok):
    problem = yacht_index.database_version_problem(config, version)
    assert (problem is None) == ok
    if not ok:
        assert config["min_yacht_version"] in problem and version in problem


def test_train_records_the_yacht_version(trained_db):
    with open(trained_db) as f:
        config = json.load(f)
    assert config["yacht_version"] == utils.__version__
    assert "min_yacht_version" not in config  # trained databases keep the sketches that YACHT 1.x needs


def test_run_accepts_a_database_for_this_version(tmp_path, trained_db):
    config = _copy_db(trained_db, tmp_path, min_yacht_version=utils.__version__)
    _run(config, tmp_path)
    assert (tmp_path / "results" / "result_all.txt").exists()


def test_run_refuses_a_database_for_a_newer_version(tmp_path, trained_db):
    config = _copy_db(trained_db, tmp_path, min_yacht_version="99.0.0", yacht_version="99.0.0")
    with pytest.raises(SystemExit) as e:
        _run(config, tmp_path)
    assert e.value.code == 1
    assert not (tmp_path / "results" / "result_all.txt").exists()


def test_index_tool_refuses_a_database_for_a_newer_version(tmp_path, trained_db):
    config = _copy_db(trained_db, tmp_path, min_yacht_version="99.0.0")
    with pytest.raises(ValueError, match="requires YACHT 99.0.0"):
        yacht_index.ensure_index(config)
    with pytest.raises(SystemExit):
        yacht_index.main(["--json", config])
    assert yacht_index.ensure_index_in_folder(os.path.dirname(config)) is None


def test_run_finds_the_files_of_a_manually_extracted_database(tmp_path, trained_db):
    # pretrained zips hold configs with paths like /<folder>/<file>, which `yacht download` rewrites; a database
    # extracted by hand (or moved) is found next to its config file
    config = _copy_db(trained_db, tmp_path)
    with open(config) as f:
        entries = json.load(f)
    for key in ("manifest_file_path", "intermediate_files_dir", "index_file_path"):
        entries[key] = os.path.join("/db_pretrained", os.path.basename(entries[key]))
    with open(config, "w") as f:
        json.dump(entries, f)
    _run(config, tmp_path)
    assert (tmp_path / "results" / "result_all.txt").exists()
    assert yacht_index.ensure_index(config) == str(tmp_path / "db" / "db_index.yidx")
