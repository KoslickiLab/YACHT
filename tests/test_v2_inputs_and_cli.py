"""
Tests of the YACHT 2.0 input checks and CLI changes: sketches that are not valid FracMinHash sketches, samples
sketched with a smaller scaled than the database, `yacht train` inputs that crashed older versions (empty sketches,
duplicate names, several k-mer sizes), the comparison passes, the Excel row limit, the deprecated
`yacht run --num_threads`, `yacht sketch ref` errors and the database download options.
"""
import argparse
import gzip
import io
import json
import os
import subprocess
import zipfile
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import sourmash
from sourmash.minhash import _get_max_hash_for_scaled

from yacht import download_default_ref_db, download_pretrained_ref_db, utils
from yacht import make_training_data_from_sketches, run_YACHT, sketch_ref_genomes

TESTDATA = os.path.join(os.path.dirname(__file__), "testdata")
REF_ZIP = os.path.join(TESTDATA, "20_genomes_sketches.zip")
SAMPLE_ZIP = os.path.join(TESTDATA, "sample.sig.zip")


# ---------------------------------------------------------------------------------------------- sample sketches

def _save_sketch(path, minhash, name="sample"):
    sig = sourmash.SourmashSignature(minhash, name=name)
    with sourmash.save_load.SaveSignaturesToLocation(str(path)) as save:
        save.add(sig)
    return str(path)


def _random_minhash(scaled=1000, n=300, seed=1, abund=True):
    rng = np.random.default_rng(seed)
    max_hash = _get_max_hash_for_scaled(scaled)
    mh = sourmash.MinHash(0, 31, scaled=scaled, track_abundance=abund)
    hashes = sorted(set(int(h) for h in rng.integers(1, max_hash, size=n, dtype=np.uint64)))
    if abund:
        mh.set_abundances({h: int(a) for h, a in zip(hashes, rng.integers(1, 50, size=len(hashes)))})
    else:
        mh.add_many(hashes)
    return mh


def _edit_sketch_json(src, dst, edit):
    """Copy a one-signature sourmash zip, changing the sketch JSON with edit(sketch dict) (the manifest is kept)."""
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(dst, "w") as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename.startswith("signatures/"):
                sigs = json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)
                edit(sigs[0]["signatures"][0])
                data = gzip.compress(json.dumps(sigs).encode())
            zout.writestr(item, data)
    return str(dst)


@pytest.mark.parametrize(
    "edit, reason",
    [
        (lambda s: (s["mins"].insert(1, s["mins"][1]), s["abundances"].insert(1, 3)), "repeated hash values"),
        (lambda s: (s["mins"].append(s["max_hash"] + 1), s["abundances"].append(1)), "above its max_hash"),
        (lambda s: s["abundances"].pop(), "abundances for"),
        (lambda s: s["abundances"].__setitem__(0, 0), "abundance 0"),
    ],
)
def test_invalid_sample_sketches_are_rejected(tmp_path, edit, reason):
    good = _save_sketch(tmp_path / "good.sig.zip", _random_minhash())
    bad = _edit_sketch_json(good, tmp_path / "bad.sig.zip", edit)
    with pytest.raises(utils.InvalidSketchError, match=reason):
        utils.load_sample_sketch(bad, 31)


def test_unsorted_distinct_hashes_are_a_valid_sketch(tmp_path):
    mh = _random_minhash()
    good = _save_sketch(tmp_path / "good.sig.zip", mh)

    def reverse(s):
        s["mins"].reverse()
        s["abundances"].reverse()

    unsorted = _edit_sketch_json(good, tmp_path / "unsorted.sig.zip", reverse)
    _, hashes, info, _ = utils.load_sample_sketch(unsorted, 31)
    assert hashes.tolist() == sorted(mh.hashes)
    assert info[4] == len(mh)


def test_huge_abundances_mean_like_sourmash(tmp_path):
    mh = sourmash.MinHash(0, 31, scaled=1000, track_abundance=True)
    mh.set_abundances({11: 2**63 + 5, 22: 2**64 - 1, 33: 7})
    path = _save_sketch(tmp_path / "huge.sig.zip", mh)
    _, _, info, _ = utils.load_sample_sketch(path, 31)
    assert info[3] == mh.mean_abundance == np.mean(np.array([2**63 + 5, 2**64 - 1, 7], dtype=np.uint64))


def test_finer_sample_is_downsampled_to_the_database_scaled(tmp_path):
    mh = _random_minhash(scaled=100, n=5000, seed=3)
    path = _save_sketch(tmp_path / "s100.sig.zip", mh)
    expected = mh.downsample(scaled=1000)
    _, hashes, info, _ = utils.load_sample_sketch(path, 31, target_scaled=1000)
    assert hashes.tolist() == sorted(expected.hashes)
    assert info[3] == expected.mean_abundance
    assert info[4:] == (len(expected), 1000)
    # the same sketch at scaled 1000 loads identically
    path1000 = _save_sketch(tmp_path / "s1000.sig.zip", expected)
    _, hashes1000, info1000, _ = utils.load_sample_sketch(path1000, 31, target_scaled=1000)
    assert hashes1000.tolist() == hashes.tolist() and info1000[3:] == info[3:]


def test_coarser_sample_is_an_error(tmp_path):
    path = _save_sketch(tmp_path / "s10000.sig.zip", _random_minhash(scaled=10000))
    with pytest.raises(ValueError, match="--scaled 1000"):
        utils.load_sample_sketch(path, 31, target_scaled=1000)


# ------------------------------------------------------------------------------------------------- yacht train

def _train(ref_zip, outdir, prefix="db"):
    args = argparse.Namespace(
        ref_file=str(ref_zip), ksize=31, num_threads=2, ani_thresh=0.95, prefix=prefix, outdir=str(outdir), force=True
    )
    make_training_data_from_sketches.main(args)
    return os.path.join(str(outdir), f"{prefix}_config.json")


@pytest.fixture(scope="module")
def trained_db(tmp_path_factory):
    return _train(REF_ZIP, tmp_path_factory.mktemp("db"))


def _sigs(zip_path):
    return list(sourmash.load_file_as_signatures(zip_path))


def _save_sigs(path, sigs):
    with sourmash.save_load.SaveSignaturesToLocation(str(path)) as save:
        for sig in sigs:
            save.add(sig)
    return str(path)


def test_train_skips_empty_sketches(tmp_path, trained_db):
    empty = sourmash.SourmashSignature(sourmash.MinHash(0, 31, scaled=1000, track_abundance=True), name="empty genome")
    ref = _save_sigs(tmp_path / "ref.sig.zip", _sigs(REF_ZIP) + [empty])
    config = _train(ref, tmp_path)
    with open(trained_db) as f:
        expected = json.load(f)
    with open(config) as f:
        got = json.load(f)
    assert open(got["manifest_file_path"]).read() == open(expected["manifest_file_path"]).read()
    assert open(got["index_file_path"], "rb").read() == open(expected["index_file_path"], "rb").read()


def test_train_rejects_duplicate_names(tmp_path):
    sigs = _sigs(REF_ZIP)
    sigs[1] = sourmash.SourmashSignature(sigs[1].minhash, name=sigs[0].name)
    ref = _save_sigs(tmp_path / "ref.sig.zip", sigs)
    with pytest.raises(ValueError, match="must be unique"):
        _train(ref, tmp_path)


def test_train_rejects_invalid_sketches(tmp_path):
    sig = _sigs(REF_ZIP)[0]
    good = _save_sigs(tmp_path / "one.sig.zip", [sig])

    def add_hash_above_max_hash(s):
        s["mins"].append(s["max_hash"] + 1)
        if s.get("abundances") is not None:
            s["abundances"].append(1)

    bad = _edit_sketch_json(good, tmp_path / "bad.sig.zip", add_hash_above_max_hash)
    with pytest.raises(utils.InvalidSketchError, match="above its max_hash"):
        _train(bad, tmp_path)


def test_train_uses_only_the_requested_ksize(tmp_path, trained_db):
    sigs = _sigs(REF_ZIP)
    other_k = []
    for s in sigs:
        mh = sourmash.MinHash(0, 21, scaled=1000)
        mh.add_many(list(s.minhash.hashes)[:50])
        other_k.append(sourmash.SourmashSignature(mh, name=s.name))
    ref = _save_sigs(tmp_path / "multik.sig.zip", sigs + other_k)
    config = _train(ref, tmp_path)
    with open(trained_db) as f:
        expected = json.load(f)
    with open(config) as f:
        got = json.load(f)
    assert open(got["manifest_file_path"]).read() == open(expected["manifest_file_path"]).read()


def test_comparison_passes_do_not_change_the_selection(tmp_path, trained_db):
    with open(trained_db) as f:
        expected = json.load(f)
    # a matrix budget of ~3 genomes per pass
    with patch.object(utils, "MAX_COMPARISON_MATRIX_BYTES", 4 * 20 * 3):
        config = _train(REF_ZIP, tmp_path)
    with open(config) as f:
        got = json.load(f)
    log_dir = os.path.join(got["intermediate_files_dir"], "comparison_files")
    assert len({name.split("_")[0] for name in os.listdir(log_dir)}) > 1  # several passes were made
    assert open(got["manifest_file_path"]).read() == open(expected["manifest_file_path"]).read()


# --------------------------------------------------------------------------------------------------- yacht run

def _run_args(config, outdir, **kwargs):
    args = dict(
        json=config, sample_file=SAMPLE_ZIP, significance=0.99, num_threads=None, keep_raw=False, show_all=False,
        min_coverage_list=[1, 0.5], outdir=str(outdir),
    )
    args.update(kwargs)
    return argparse.Namespace(**args)


def test_run_reports_the_number_of_sample_hashes(tmp_path, trained_db):
    run_YACHT.main(_run_args(trained_db, tmp_path))
    df = pd.read_csv(tmp_path / "results" / "result_all.txt", sep="\t")
    n_hashes = len(_sigs(SAMPLE_ZIP)[0].minhash)
    assert (df["num_exclusive_kmers_in_sample_sketch"] == n_hashes).all()


def test_run_num_threads_is_deprecated(tmp_path, trained_db):
    parser = argparse.ArgumentParser()
    run_YACHT.add_arguments(parser)
    assert "--num_threads" not in parser.format_help()
    args = parser.parse_args(["--json", trained_db, "--sample_file", SAMPLE_ZIP, "--num_threads", "0", "--outdir", str(tmp_path)])
    run_YACHT.main(args)  # ignored, even if < 1
    assert (tmp_path / "results" / "result.xlsx").exists()


@pytest.mark.parametrize("significance", [0, 1, 1.5])
def test_run_rejects_significance_outside_0_1(tmp_path, trained_db, significance):
    with pytest.raises(SystemExit):
        run_YACHT.main(_run_args(trained_db, tmp_path, significance=significance))


def test_run_excel_row_limit(tmp_path, trained_db):
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "result.xlsx").write_text("from an earlier run")
    # (only the genomes that share hashes with the sample are reported: one here)
    with patch.object(utils, "EXCEL_MAX_ROWS", 1), pytest.raises(SystemExit):
        run_YACHT.main(_run_args(trained_db, tmp_path, show_all=True))
    assert (tmp_path / "results" / "result_all.txt").exists()
    assert not (tmp_path / "results" / "result.xlsx").exists()


def test_run_rejects_other_ksize(tmp_path, trained_db):
    sample = _save_sketch(tmp_path / "k21.sig.zip", sourmash.MinHash(0, 21, scaled=1000, track_abundance=True))
    with pytest.raises(SystemExit):
        run_YACHT.main(_run_args(trained_db, tmp_path, sample_file=sample))


# ------------------------------------------------------------------------------------------ yacht sketch ref

def test_sketch_ref_missing_input_exits_1(tmp_path):
    res = subprocess.run(["yacht", "sketch", "ref", "--infile", str(tmp_path / "missing"), "--outfile", str(tmp_path / "o.zip")])
    assert res.returncode == 1


def test_sketch_ref_does_not_write_into_the_input_folder(tmp_path):
    folder = tmp_path / "genomes"
    folder.mkdir()
    (folder / "g1.fasta").write_text(">s\n" + "ACGTTGCA" * 40 + "\n")
    sketch_ref_genomes.sketch_multiple_files(folder, 31, 1000, tmp_path / "ref.sig.zip")
    assert os.listdir(folder) == ["g1.fasta"]
    assert [s.name for s in _sigs(str(tmp_path / "ref.sig.zip"))] == ["g1"]


# ---------------------------------------------------------------------------------------------- yacht download

def test_default_ref_db_new_versions():
    url = download_default_ref_db.generate_download_url
    base = utils.BASE_URL
    ns = argparse.Namespace
    assert url(ns(database="gtdb", db_version="rs226", gtdb_type="reps", k=31)) == f"{base}gtdb-rs226/gtdb-rs226-reps.k31.sig.zip"
    assert url(ns(database="gtdb", db_version="rs226", gtdb_type=None, k=21)) == f"{base}gtdb-rs226/gtdb-rs226.k21.sig.zip"
    assert url(ns(database="gtdb", db_version="rs220", gtdb_type="reps", k=51)) == f"{base}gtdb-rs220/gtdb-reps-rs220-k51.zip"
    assert url(ns(database="genbank", db_version="genbank-2025.04", ncbi_organism="fungi", k=31)) == (
        f"{base}genbank-2025.04/genbank-20250408-fungi-k31.zip"
    )
    assert url(ns(database="genbank", db_version="genbank-2025.04", ncbi_organism="bacteria", k=31)) is None


ZENODO_RECORDS = [
    {"doi": "10.5281/zenodo.1", "files": [
        {"key": "gtdb-rs214-reps.k31_0.95_pretrained.zip", "size": 2e9, "links": {"self": "https://z/1"}},
    ]},
    {"doi": "10.5281/zenodo.2", "files": [
        {"key": "gtdb-rs214-reps.k31_0.80_pretrained.zip", "size": 4e9, "links": {"self": "https://z/2"}},
    ]},
    {"doi": "10.5281/zenodo.3", "files": [
        {"key": "ictv-vmr-msl40.k31-sc50_0.95_pretrained.zip", "size": 1e8, "links": {"self": "https://z/3"}},
    ]},
    {"doi": "10.5281/zenodo.4", "files": [{"key": "gtdb_ani_thresh_rs214-reps_k31_0.95.zip", "size": 1}]},
]


def _download_args(**kwargs):
    args = dict(database=None, db_version=None, ncbi_organism=None, k=31, ani_thresh=0.95, name=None, list=False, outfolder=".")
    args.update(kwargs)
    return argparse.Namespace(**args)


def test_pretrained_databases_are_listed_from_zenodo(capsys):
    dbs = download_pretrained_ref_db.list_pretrained_databases(ZENODO_RECORDS)
    assert sorted(dbs) == ["gtdb-rs214-reps.k31", "ictv-vmr-msl40.k31-sc50"]
    assert sorted(dbs["gtdb-rs214-reps.k31"]) == [0.8, 0.95]
    with patch.object(download_pretrained_ref_db, "fetch_zenodo_records", return_value=ZENODO_RECORDS):
        download_pretrained_ref_db.main(_download_args(list=True, ani_thresh=None))
    out = capsys.readouterr().out
    assert "gtdb-rs214-reps.k31" in out and "ANI 0.8, 0.95" in out


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        (dict(database="gtdb", db_version="rs214"), "gtdb-rs214-reps.k31_0.95_pretrained.zip"),
        (dict(name="gtdb-rs214-reps.k31", ani_thresh=0.8), "gtdb-rs214-reps.k31_0.80_pretrained.zip"),
        (dict(database="ictv-vmr", db_version="msl40"), "ictv-vmr-msl40.k31-sc50_0.95_pretrained.zip"),
    ],
)
def test_pretrained_database_selection(tmp_path, kwargs, expected):
    with patch.object(download_pretrained_ref_db, "fetch_zenodo_records", return_value=ZENODO_RECORDS), \
            patch.object(download_pretrained_ref_db, "download_file", return_value=True) as download, \
            patch.object(download_pretrained_ref_db, "unzip_file") as unzip:
        download_pretrained_ref_db.main(_download_args(outfolder=str(tmp_path), **kwargs))
    assert download.call_args[0][1] == os.path.join(str(tmp_path), expected)
    unzip.assert_called_with(expected, str(tmp_path))


@pytest.mark.parametrize(
    "kwargs", [dict(name="gtdb-rs214-reps.k31", ani_thresh=0.9995), dict(name="nonexistent"), dict(database="gtdb")]
)
def test_pretrained_database_errors(tmp_path, kwargs):
    with patch.object(download_pretrained_ref_db, "fetch_zenodo_records", return_value=ZENODO_RECORDS), \
            pytest.raises(SystemExit):
        download_pretrained_ref_db.main(_download_args(outfolder=str(tmp_path), **kwargs))
