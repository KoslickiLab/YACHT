#!/usr/bin/env python
"""
Binary index of a trained YACHT reference database (``<prefix>_index.yidx``).

`yacht train` writes the index next to ``<prefix>_config.json`` and ``<prefix>_processed_manifest.tsv``, and
records its location in the config file under ``index_file_path``. `yacht run` memory-maps the index instead
of parsing every genome signature of the database, which is what makes it fast.

Databases trained with YACHT < 2.0 (and pretrained databases downloaded before the index was added to them)
have no index. It can be added in place, without the original reference sketches, with:

    python -m yacht.yacht_index --json <prefix>_config.json
"""
import os
import re
import sys
import csv
import io
import json
import hashlib
import argparse
from typing import List, Optional, Tuple
import numpy as np
import pandas as pd
from loguru import logger

INDEX_SUFFIX = "_index.yidx"
CONFIG_KEY = "index_file_path"
# config entries: the YACHT version that made the database, and the oldest YACHT version that can use it (set e.g. in
# the pretrained databases of YACHT 2.0, which leave out the genome sketches that YACHT 1.x needs)
YACHT_VERSION_KEY = "yacht_version"
MIN_YACHT_VERSION_KEY = "min_yacht_version"
MANIFEST_SUFFIX = "_processed_manifest.tsv"
SIG_SUFFIX = ".sig"


def _core():
    from . import _yacht_core

    return _yacht_core


def manifest_digest(md5sums) -> bytes:
    """
    Fingerprint of a processed manifest: SHA-256 of its md5sum column (in row order).
    :param md5sums: iterable of genome signature md5sums (manifest row order)
    :return: 32 bytes
    """
    return hashlib.sha256("\n".join(str(x) for x in md5sums).encode()).digest()


def default_index_path(manifest_file_path: str) -> str:
    """
    Location of the index that belongs to a processed manifest (``<prefix>_processed_manifest.tsv`` ->
    ``<prefix>_index.yidx``, in the same folder).
    """
    if manifest_file_path.endswith(MANIFEST_SUFFIX):
        return manifest_file_path[: -len(MANIFEST_SUFFIX)] + INDEX_SUFFIX
    return os.path.splitext(manifest_file_path)[0] + INDEX_SUFFIX


def candidate_index_paths(config: dict, json_file_path: Optional[str] = None) -> List[str]:
    """
    Possible locations of the index of a database: the one recorded in the config, the one derived from the
    manifest path, and the same file names next to the config file (for databases that were moved).
    """
    paths = []
    if isinstance(config.get(CONFIG_KEY), str):
        paths.append(config[CONFIG_KEY])
    paths.append(default_index_path(config["manifest_file_path"]))
    if json_file_path:
        config_dir = os.path.dirname(os.path.abspath(json_file_path))
        paths += [os.path.join(config_dir, os.path.basename(p)) for p in list(paths)]
    unique = []
    for p in paths:
        if p not in unique:
            unique.append(p)
    return unique


def genome_sig_paths(manifest: pd.DataFrame, intermediate_files_dir: str) -> List[str]:
    """
    Signature file of every manifest row (``<intermediate_files_dir>/signatures/<md5sum>.sig``, or ``.sig.gz``
    for databases whose signatures were never decompressed).
    """
    sig_dir = os.path.join(intermediate_files_dir, "signatures")
    paths = []
    for md5sum in manifest["md5sum"]:
        path = os.path.join(sig_dir, f"{md5sum}{SIG_SUFFIX}")
        if not os.path.exists(path) and os.path.exists(path + ".gz"):
            path += ".gz"
        paths.append(path)
    return paths


def build_index(
    manifest: pd.DataFrame, intermediate_files_dir: str, ksize: int, scale: int, index_file_path: str
) -> dict:
    """
    Build the index of a trained database and write it (atomically) to index_file_path.
    :param manifest: processed manifest (dataframe with at least the md5sum column), in the order `yacht run` reads it
    :param intermediate_files_dir: the `<prefix>_intermediate_files` folder holding signatures/<md5sum>.sig
    :param ksize: k-mer size
    :param scale: scale factor of the genome sketches
    :param index_file_path: output path
    :return: dictionary of index statistics
    """
    logger.info(f"Building the YACHT index of {len(manifest)} genomes: {index_file_path}")
    stats = _core().build_index(
        genome_sig_paths(manifest, intermediate_files_dir),
        int(ksize),
        int(scale),
        manifest_digest(manifest["md5sum"]),
        index_file_path,
    )
    logger.info(
        f"YACHT index built: {stats['n_hashes']} hashes, {stats['n_keys']} distinct, "
        f"{stats['n_sets']} shared hash sets, {stats['file_len'] / 1e6:.1f} MB."
    )
    return stats


def index_problem(index, manifest: pd.DataFrame, ksize: int, scale: int) -> Optional[str]:
    """
    Check that an opened index belongs to this database. Returns None if it does, otherwise a description of
    the mismatch.
    """
    if index.ksize != ksize or index.scaled != scale:
        return f"the index was built for ksize={index.ksize}, scale={index.scaled} (database: ksize={ksize}, scale={scale})"
    if index.n_rows != len(manifest):
        return f"the index has {index.n_rows} genomes but the manifest has {len(manifest)}"
    if index.digest != manifest_digest(manifest["md5sum"]):
        return "the index does not match the genomes listed in the manifest"
    sizes = index.sketch_sizes()
    if "num_unique_kmers_in_genome_sketch" in manifest.columns and not np.array_equal(
        sizes, manifest["num_unique_kmers_in_genome_sketch"].to_numpy()
    ):
        return "the genome sketch sizes in the index do not match the manifest"
    return None


def find_index(
    config: dict, manifest: pd.DataFrame, json_file_path: Optional[str] = None
) -> Tuple[Optional[object], Optional[str], List[str]]:
    """
    Locate and open a valid index for a database.
    :return: (opened index or None, its path or None, list of problems found with the candidate files)
    """
    problems = []
    for path in candidate_index_paths(config, json_file_path):
        if not os.path.exists(path):
            continue
        try:
            index = _core().YachtIndex(path)
        except ValueError as e:
            problems.append(f"{path}: {e}")
            continue
        problem = index_problem(index, manifest, config["ksize"], config["scale"])
        if problem is None:
            return index, path, problems
        problems.append(f"{path}: {problem}")
    return None, None, problems


def locate_database_files(config: dict, json_file_path: str) -> dict:
    """
    The config with each database file path that does not exist replaced by the same file name next to the config
    file, if that exists: databases that were moved, or extracted from a pretrained zip without `yacht download`
    (whose configs hold paths like ``/<folder>/<file>``).
    """
    config_dir = os.path.dirname(os.path.abspath(json_file_path))
    located = dict(config)
    for key in ("manifest_file_path", "intermediate_files_dir", CONFIG_KEY):
        path = config.get(key)
        if isinstance(path, str) and not os.path.exists(path):
            candidate = os.path.join(config_dir, os.path.basename(path))
            if os.path.exists(candidate):
                located[key] = candidate
    return located


def _version_tuple(version) -> Tuple[int, ...]:
    return tuple(int(part) for part in re.findall(r"\d+", str(version))[:3])


def database_version_problem(config: dict, yacht_version: str) -> Optional[str]:
    """
    Check that this YACHT version can use the database described by config. Databases without
    `min_yacht_version` (made by `yacht train`, or with YACHT < 2.0) can be used by any YACHT 2.x.
    :return: None, or a description of the problem
    """
    required = config.get(MIN_YACHT_VERSION_KEY)
    if not required or _version_tuple(yacht_version) >= _version_tuple(required):
        return None
    made_with = config.get(YACHT_VERSION_KEY)
    return (
        f"This database requires YACHT {required} or later"
        + (f" (it was made with YACHT {made_with})" if made_with else "")
        + f", but this is YACHT {yacht_version}. Please update YACHT, e.g. with "
        "`conda update -c conda-forge -c bioconda yacht`."
    )


def _installed_yacht_version() -> str:
    from .utils import __version__
    return __version__


def open_index_or_exit(config: dict, json_file_path: str, manifest: pd.DataFrame):
    """
    Open the index of the database described by config; print instructions and exit with status 1 if the
    database has no (valid) index, e.g. because it was trained with YACHT < 2.0.
    """
    index, _, problems = find_index(config, manifest, json_file_path)
    if index is not None:
        return index
    expected = candidate_index_paths(config, json_file_path)[0]
    if problems:
        print("The YACHT index of this database is out of date or invalid:\n")
        for p in problems:
            print(f"    {p}")
        print()
    else:
        print(
            f"No YACHT index found for this database (expected {expected}).\n"
            f"This database was likely created with a YACHT version older than 2.0.\n"
        )
    print(
        "Starting with YACHT 2.0, `yacht run` requires a binary index of the reference database. "
        "Create it (in place, from the existing database files) with:\n\n"
        f"    python -m yacht.yacht_index --json {json_file_path}\n\n"
        "Alternatively, re-run `yacht train` on the reference sketches, or re-download pretrained databases with "
        "`yacht download pretrained_ref_db` (which creates the index automatically). Exiting..\n"
    )
    sys.exit(1)


def ensure_index(json_file_path: str, force: bool = False) -> str:
    """
    Make sure the database described by a config file has a valid index: build it if it is missing, stale or
    force=True, and record its location in the config file.
    :return: path of the index
    """
    with open(json_file_path) as fp:
        config = locate_database_files(json.load(fp), json_file_path)
    problem = database_version_problem(config, _installed_yacht_version())
    if problem:
        raise ValueError(problem)
    manifest = pd.read_csv(config["manifest_file_path"], sep="\t", header=0)
    index, path, problems = (None, None, []) if force else find_index(config, manifest, json_file_path)
    if index is not None:
        logger.info(f"Found a valid YACHT index: {path}")
    else:
        for p in problems:
            logger.warning(f"Rebuilding invalid index: {p}")
        path = default_index_path(config["manifest_file_path"])
        build_index(manifest, config["intermediate_files_dir"], config["ksize"], config["scale"], path)
    if config.get(CONFIG_KEY) != path:
        config[CONFIG_KEY] = path
        with open(json_file_path, "w") as fp:
            json.dump(config, fp, indent=4)
    return path


def ensure_index_in_folder(folder: str) -> Optional[str]:
    """
    Make sure the (pretrained) database in a folder has a valid index. Never raises: problems are logged.
    :param folder: folder holding a `*_config.json` file
    :return: path of the index, or None
    """
    try:
        if not os.path.isdir(folder):
            return None
        configs = [f for f in os.listdir(folder) if f.endswith("_config.json")]
        if not configs:
            return None
        with open(os.path.join(folder, configs[0])) as fp:
            problem = database_version_problem(json.load(fp), _installed_yacht_version())
        if problem:
            logger.error(problem)
            return None
        return ensure_index(os.path.join(folder, configs[0]))
    except Exception as e:  # noqa: BLE001 - downloading must not fail because of the index
        logger.error(
            f"Could not create the YACHT index for the database in {folder}: {e}. "
            f"Run `python -m yacht.yacht_index --json <config.json>` to create it."
        )
        return None


def nontrivial_organisms(manifest: pd.DataFrame, names: List[str]) -> pd.Series:
    """
    Boolean mask of the manifest rows `yacht run` analyses: organisms whose name equals the name of a genome
    sharing at least one hash with the sample. Names are round-tripped through CSV and pandas exactly like the
    output of `sourmash scripts multisearch` was in YACHT < 2.0, so that the selection is identical.
    :param manifest: processed manifest
    :param names: signature names of the genomes that share at least one hash with the sample
    """
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["query_name", "match_name"])
    for name in names:
        writer.writerow(["sample", name])
    buffer.seek(0)
    multisearch_result = pd.read_csv(buffer, sep=",", header=0)
    multisearch_result = multisearch_result.drop_duplicates().reset_index(drop=True)
    return manifest["organism_name"].isin(multisearch_result["match_name"].to_list())


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="python -m yacht.yacht_index",
        description="Create (or verify) the binary index that `yacht run` needs for a trained YACHT database. "
        "Use this for databases created with YACHT < 2.0.",
    )
    parser.add_argument("--json", required=True, help="Path to the <prefix>_config.json file of the database.")
    parser.add_argument("--force", action="store_true", help="Rebuild the index even if a valid one exists.")
    args = parser.parse_args(args)
    if not os.path.exists(args.json):
        print(f"Config file {args.json} does not exist.")
        sys.exit(1)
    try:
        path = ensure_index(os.path.abspath(args.json), force=args.force)
    except ValueError as e:
        print(e)
        sys.exit(1)
    print(f"YACHT index ready: {path}")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")
    main()
