import os
import sys
import numpy as np
import sourmash
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool
from loguru import logger
from typing import Optional, List, Set, Dict, Tuple
import collections
from scipy.sparse import lil_matrix
from glob import glob

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)

# Set up contants
COL_NOT_FOUND_ERROR = "Column not found: {}"

# Set up global variables
__version__ = "1.2.3"
GITHUB_API_URL = "https://api.github.com/repos/KoslickiLab/YACHT/contents/demo/{path}"
GITHUB_RAW_URL = "https://raw.githubusercontent.com/KoslickiLab/YACHT/main/demo/{path}"
BASE_URL = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/"
ZENODO_COMMUNITY_URL = "https://zenodo.org/api/records/?communities=yacht&size=100"

def load_signature_with_ksize(filename: str, ksize: int) -> sourmash.SourmashSignature:
    """
    Helper function that loads the signature for a given kmer size from the provided signature file.
    Filename should point to a sourmash signature file. Raises exception if given kmer size is not present in the file.
    This is a wrapper of sourmash.load_file_as_signatures, and accept all types of format: .sig, .sig.zip, .sbt, .lca, and .sqldb.
    However, this function specifically ask for 1 signature so lca format is not appropriate; as of sourmash v4.8, sqldb doesn't accept "abund" parameter for signatures.
    :param filename: string (location of the signature file of any format: .sig, .sig.zip, .sbt, .lca, and .sqldb)
    :param ksize: int (size of kmer)
    :return: sourmash signature
    """
    # Take the first sample signature with the given kmer size
    sketches = list(sourmash.load_file_as_signatures(filename, ksize=ksize))
    if len(sketches) != 1:
        raise ValueError(
            f"Expected exactly one signature with ksize {ksize} in {filename}, found {len(sketches)}"
        )
    if len(sketches[0].minhash.hashes) == 0:
        raise ValueError(
            "Empty sketch in signature. This may be due to too high of a scale factor, please reduce it, eg. --scaled=1, and try again."
        )
    return sketches[0]


def get_num_kmers(
    minhash_mean_abundance: Optional[float],
    minhash_hashes_len: int,
    minhash_scaled: int,
    scale: bool = True,
) -> int:
    """
    Helper function that estimates the total number of kmers in a given sample.
    :param minhash_mean_abundance: float or None (mean abundance of the signature)
    :param minhash_hashes_len: int (number of hashes in the signature)
    :param minhash_scaled: int (scale factor of the signature)
    :param scale: bool (whether to scale the number of kmers by the scale factor)
    :return: int (estimated total number of kmers)
    """
    # Abundances may not have been kept, in which case, just use 1
    if minhash_mean_abundance:
        num_kmers = minhash_mean_abundance * minhash_hashes_len
    else:
        num_kmers = minhash_hashes_len
    if scale:
        num_kmers *= minhash_scaled
    return int(np.round(num_kmers))


def check_file_existence(file_path: str, error_description: str) -> None:
    """
    Helper function that checks if a file exists. If not, raises a ValueError with the given error description.
    :param file_path: string (location of the file)
    :param error_description: string (description of the error)
    :return: None
    """
    if not os.path.exists(file_path):
        raise ValueError(error_description)


def get_info_from_single_sig(
    sig_file: str, ksize: int
) -> Tuple[str, str, float, int, int]:
    """
    Helper function that gets signature information (name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled) from a single sourmash signature file.
    :param sig_file: string (location of the signature file with .sig.gz format)
    :param ksize: int (size of kmer)
    :return: tuple (name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled)
    """
    try:
        sig = load_signature_with_ksize(sig_file, ksize)
        return (
            sig.name,
            sig.md5sum(),
            sig.minhash.mean_abundance,
            len(sig.minhash.hashes),
            sig.minhash.scaled,
        )
    except:
        logger.warning(f"CANNOT extract the relevant info from the signature file: {sig_file}")
        return None

def run_multisearch(
    num_threads: int, ani_thresh: float, ksize: int, scale: int, path_to_temp_dir: str
) -> Dict[str, List[str]]:
    """
    Helper function that runs the sourmash multisearch to find the close related genomes.
    :param num_threads: int (number of threads to use)
    :param ani_thresh: float (threshold for ANI, below which we consider two organisms to be distinct)
    :param ksize: int (size of kmer)
    :param scale: int (scale factor)
    :param path_to_temp_dir: string (path to the folder to store the intermediate files)
    :return: a dataframe with symmetric pairwise multisearch result (query_name, match_name)
    """

    # run the sourmash multisearch
    # save signature files to a text file
    sig_files = pd.DataFrame(
        [
            os.path.join(path_to_temp_dir, "signatures", file)
            for file in os.listdir(os.path.join(path_to_temp_dir, "signatures"))
        ]
    )
    sig_files_path = os.path.join(path_to_temp_dir, "training_sig_files.txt")
    sig_files.to_csv(sig_files_path, header=False, index=False)

    # convert ani threshold to containment threshold
    containment_thresh = ani_thresh**ksize
    cmd = f"sourmash scripts multisearch {sig_files_path} {sig_files_path} -k {ksize} -s {scale} -c {num_threads} -t {containment_thresh} -o {os.path.join(path_to_temp_dir, 'training_multisearch_result.csv')}"
    logger.info(f"Running sourmash multisearch with command: {cmd}")
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise ValueError(f"Error running sourmash multisearch with command: {cmd}")

    # read the multisearch result
    multisearch_result = pd.read_csv(
        os.path.join(path_to_temp_dir, "training_multisearch_result.csv"),
        sep=",",
        header=0,
    )
    multisearch_result = multisearch_result.query(
        "query_name != match_name"
    ).reset_index(drop=True)

    # because the multisearch result is not symmetric, that is
    # we have: A B score but not B A score
    # we need to make it symmetric
    A_TO_B = (
        multisearch_result[["query_name", "match_name"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    B_TO_A = A_TO_B[["match_name", "query_name"]].rename(
        columns={"match_name": "query_name", "query_name": "match_name"}
    )
    multisearch_result = (
        pd.concat([A_TO_B, B_TO_A]).drop_duplicates().reset_index(drop=True)
    )

    # change column type to string
    multisearch_result["query_name"] = multisearch_result["query_name"].astype(str)
    multisearch_result["match_name"] = multisearch_result["match_name"].astype(str)

    return multisearch_result


def remove_corr_organisms_from_ref(
    sig_info_dict: Dict[str, Tuple[str, float, int, int]],
    multisearch_result: pd.DataFrame,
) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
    temp_remove_set: Set[str] = set()  # Annotate as a set of strings
    """
    Helper function that removes the close related organisms from the reference matrix.
    :param sig_info_dict: a dictionary mapping all signature name from reference data to a tuple (md5sum, minhash mean abundance, minhash hashes length, minhash scaled)
    :param multisearch_result: a dataframe with symmetric pairwise multisearch result (query_name, match_name)
    :return 
        remove_corr_df: a dataframe with two columns: removed organism name and its close related organisms
        manifest_df: a dataframe containing the processed reference signature information
    """
    # extract organisms that have close related organisms and their number of unique kmers
    # sort name in order to better check the removed organisms
    corr_organisms = sorted(
        [str(query_name) for query_name in multisearch_result["query_name"].unique()]
    )
    sizes = np.array([sig_info_dict[organism][2] for organism in corr_organisms])
    # sort organisms by size in ascending order, so we keep the largest organism, discard the smallest
    bysize = np.argsort(sizes)
    corr_organisms_bysize = np.array(corr_organisms)[bysize].tolist()

    # use dictionary to store the removed organisms and their close related organisms
    # key: removed organism name
    # value: a set of close related organisms
    mapping = multisearch_result.groupby("query_name")["match_name"].agg(set).to_dict()

    # remove the sorted organisms until all left genomes are distinct (e.g., ANI <= ani_thresh)
    temp_remove_set = set()
    # loop through the organisms size in ascending order
    for organism in tqdm(
        corr_organisms_bysize, desc="Removing close related organisms"
    ):
        ## for a given organism check its close related organisms, see if there are any organisms left after removing those in the remove set
        ## if so, put this organism in the remove set
        left_corr_orgs = mapping[organism].difference(temp_remove_set)
        if len(left_corr_orgs) > 0:
            temp_remove_set.add(organism)

    # generate a dataframe with two columns: removed organism name and its close related organisms
    logger.info(
        "Generating a dataframe with two columns: removed organism name and its close related organisms."
    )
    remove_corr_list = [
        (organism, ",".join(list(mapping[organism])))
        for organism in tqdm(temp_remove_set)
    ]
    remove_corr_df = pd.DataFrame(
        remove_corr_list, columns=["removed_org", "corr_orgs"]
    )

    # remove the close related organisms from the reference genome list
    manifest_df = []
    for sig_name, (
        md5sum,
        minhash_mean_abundance,
        minhash_hashes_len,
        minhash_scaled,
    ) in tqdm(sig_info_dict.items()):
        if sig_name not in temp_remove_set:
            manifest_df.append(
                (
                    sig_name,
                    md5sum,
                    minhash_hashes_len,
                    get_num_kmers(
                        minhash_mean_abundance,
                        minhash_hashes_len,
                        minhash_scaled,
                        False,
                    ),
                    minhash_scaled,
                )
            )
    manifest_df = pd.DataFrame(
        manifest_df,
        columns=[
            "organism_name",
            "md5sum",
            "num_unique_kmers_in_genome_sketch",
            "num_total_kmers_in_genome_sketch",
            "genome_scale_factor",
        ],
    )

    return remove_corr_df, manifest_df

def collect_signature_info(
    num_threads: int, ksize: int, path_to_temp_dir: str
) -> Dict[str, Tuple[str, float, int, int]]:
    """
    Helper function that collects signature information (name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled) from a sourmash signature database.
    :param num_threads: int (number of threads to use)
    :param ksize: int (size of kmer)
    :param path_to_temp_dir: string (path to the folder to store the intermediate files)
    :return: a dictionary mapping signature name to a tuple (md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled)
    """
    ## extract in parallel
    with Pool(num_threads) as p:
        signatures = p.starmap(
            get_info_from_single_sig,
            [
                (os.path.join(path_to_temp_dir, "signatures", file), ksize)
                for file in os.listdir(os.path.join(path_to_temp_dir, "signatures"))
            ],
        )

    return {sig[0]: (sig[1], sig[2], sig[3], sig[4]) for sig in tqdm(signatures) if sig}


def _read_file(file_path):
    temp_df = pd.read_csv(
        file_path,
        sep=",",
        header=None,
        usecols=[0, 1],
        names=["query", "match"]
    )
    return temp_df

def find_and_remove_close_genomes(
    all_genome_name_path: str,
    comparison_path: str,
    sig_info_dict: Dict[str, Tuple[str, float, int, int]],
    num_threads: int = 16,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Helper function that removes closely related organisms from the reference genome list.
    :param all_genome_name_path: Path to a file containing a list of all genome names.
    :param comparison_path: Path to a directory containing pairwise comparison results between genomes.
    :param num_threads: Number of threads to use for multiprocessing when reading the comparison files. Default is 16.
    :param sig_info_dict:
        A dictionary mapping each genome signature name to a tuple containing metadata: 
        (md5sum, minhash mean abundance, minhash hashes length, minhash scaled).
        - md5sum: Checksum for data integrity.
        - minhash mean abundance: The mean abundance for the genome's minhash.
        - minhash hashes length: The length of minhash hashes.
        - minhash scaled: The scaling factor for the minhash.    
    :return 
        remove_corr_df: a dataframe with two columns: removed organism name and its close related organisms
        manifest_df: a dataframe containing the processed reference signature information
    """
    # read the name mapping file
    all_genome_names = pd.read_csv(
        all_genome_name_path,
        sep="\t",
        header=None,
    )
    name_list = all_genome_names[0].to_list()
    
    # read all comparison result files using multiprocessing
    comparison_files = glob(os.path.join(comparison_path, "*.txt"))
    with Pool(num_threads) as p:
        comparison_results = p.map(_read_file, comparison_files)
    comparison_result = pd.concat(comparison_results, ignore_index=True)
    
    # extract unique organisms
    organisms = np.array(sorted(set(comparison_result["query"]).union(set(comparison_result["match"])), key=str))
    organism_index = {organism_id: idx for idx, organism_id in enumerate(organisms)}
    num_organisms = len(organisms)
    
    # create a sparse adjacency matrix to represent the relationships
    adjacency_matrix = lil_matrix((num_organisms, num_organisms), dtype=np.bool_)
    for query_id, match_id in tqdm(zip(comparison_result["query"], comparison_result["match"]), total=len(comparison_result), desc="Creating adjacency matrix"):
        i, j = organism_index[query_id], organism_index[match_id]
        adjacency_matrix[i, j] = True
        adjacency_matrix[j, i] = True
    del comparison_result # free up memory
    
    # extract sizes of organisms
    sizes = np.array([sig_info_dict[name_list[organism_id]][2] for organism_id in tqdm(organisms, desc="Extracting sizes")])
    
    # sort organisms by size in ascending order
    bysize_indices = np.argsort(sizes)
    
    # remove organisms until all left genomes are distinct
    temp_remove_set = set()
    temp_remove_set_name = set()
    removed_mask = np.zeros(num_organisms, dtype=np.bool_)
    
    for idx in tqdm(bysize_indices, desc="Removing closely related organisms"):
        if removed_mask[idx]:
            continue
        # Find related organisms that are not already removed
        related_indices = adjacency_matrix[idx].rows[0]
        related_indices = [i for i in related_indices if not removed_mask[i]]
        if len(related_indices) > 0:
            temp_remove_set.add(organisms[idx])
            temp_remove_set_name.add(name_list[organisms[idx]])
            removed_mask[idx] = True

    # generate a dataframe with two columns: removed organism name and its close related organisms
    logger.info(
        "Generating a dataframe with two columns: removed organism name and its close related organisms."
    )
    remove_corr_list = []
    for organism_id in tqdm(temp_remove_set):
        idx = organism_index[organism_id]
        related_organisms = [name_list[organisms[i]] for i in adjacency_matrix[idx].rows[0]]
        remove_corr_list.append((name_list[organism_id], "#####".join(related_organisms)))
    del temp_remove_set # free up memory
    
    remove_corr_df = pd.DataFrame(remove_corr_list, columns=["removed_org", "corr_orgs"])

    # remove the close related organisms from the reference genome list
    manifest_df = []
    for sig_name, (
        md5sum,
        minhash_mean_abundance,
        minhash_hashes_len,
        minhash_scaled,
    ) in tqdm(sig_info_dict.items(), desc="Removing close related organisms from the reference genome list"):
        if sig_name not in temp_remove_set_name:
            manifest_df.append(
                (
                    sig_name,
                    md5sum,
                    minhash_hashes_len,
                    get_num_kmers(
                        minhash_mean_abundance,
                        minhash_hashes_len,
                        minhash_scaled,
                        False,
                    ),
                    minhash_scaled,
                )
            )
    manifest_df = pd.DataFrame(
        manifest_df,
        columns=[
            "organism_name",
            "md5sum",
            "num_unique_kmers_in_genome_sketch",
            "num_total_kmers_in_genome_sketch",
            "genome_scale_factor",
        ],
    )

    return remove_corr_df, manifest_df


class Prediction:
    def __init__(self):
        self._rank = None
        self._taxid = None
        self._percentage = None
        self._taxpath = None
        self._taxpathsn = None

    @property
    def rank(self):
        return self._rank

    @rank.setter
    def rank(self, value):
        self._rank = value

    @property
    def taxid(self):
        return self._taxid

    @taxid.setter
    def taxid(self, value):
        self._taxid = value

    @property
    def percentage(self):
        return self._percentage

    @percentage.setter
    def percentage(self, value):
        self._percentage = value

    @property
    def taxpath(self):
        return self._taxpath

    @taxpath.setter
    def taxpath(self, value):
        self._taxpath = value

    @property
    def taxpathsn(self):
        return self._taxpathsn

    @taxpathsn.setter
    def taxpathsn(self, value):
        self._taxpathsn = value

    def get_dict(self):
        return self.__dict__

    def get_pretty_dict(self):
        return {
            property.split("_")[1]: value
            for property, value in self.__dict__.items()
            if property.startswith("_")
        }

    def get_metadata(self):
        return {
            "rank": self._rank,
            "taxpath": self._taxpath,
            "taxpathsn": self._taxpathsn,
        }


def get_column_indices(
    column_name_to_index: Dict[str, int]
) -> Tuple[int, int, int, int, Optional[int]]:
    """
    (thanks to https://github.com/CAMI-challenge/OPAL, this function is modified get_column_indices from its load_data.py)
    Helper function that gets the column indices for the following columns: TAXID, RANK, PERCENTAGE, TAXPATH, TAXPATHSN
    :param column_name_to_index: dictionary mapping column name to column index
    :return: indices for TAXID, RANK, PERCENTAGE, TAXPATH, TAXPATHSN
    """
    # Assuming all other indices are mandatory and only index_taxpathsn can be optional
    index_taxpathsn: Optional[int] = None  # Correctly annotated to allow None

    if "TAXID" not in column_name_to_index:
        logger.error(COL_NOT_FOUND_ERROR.format("TAXID"))
        raise RuntimeError
    if "RANK" not in column_name_to_index:
        logger.error(COL_NOT_FOUND_ERROR.format("RANK"))
        raise RuntimeError
    if "PERCENTAGE" not in column_name_to_index:
        logger.error(COL_NOT_FOUND_ERROR.format("PERCENTAGE"))
        raise RuntimeError
    if "TAXPATH" not in column_name_to_index:
        logger.error(COL_NOT_FOUND_ERROR.format("TAXPATH"))
        raise RuntimeError
    index_taxid = column_name_to_index["TAXID"]
    index_rank = column_name_to_index["RANK"]
    index_percentage = column_name_to_index["PERCENTAGE"]
    index_taxpath = column_name_to_index["TAXPATH"]
    if "TAXPATHSN" in column_name_to_index:
        index_taxpathsn = column_name_to_index["TAXPATHSN"]
    else:
        index_taxpathsn = None
    return index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn


def get_cami_profile(
    cami_content: List[str],
) -> List[Tuple[str, Dict[str, str], List[Prediction]]]:
    header: Dict[str, str] = {}  # Dictionary mapping strings to strings
    profile: List[Prediction] = []  # List of Prediction objects
    predictions_dict: Dict[
        str, Prediction
    ] = {}  # Mapping from string to Prediction object
    """
    (thanks to https://github.com/CAMI-challenge/OPAL, this function is modified open_profile_from_tsv from its load_data.py)
    Helper function that opens a CAMI profile file and returns sample profiling information.
    params:cami_content: list of strings (lines of the CAMI profile file)
    return: list of tuples (sample_id, header, profile)
    """
    header = {}
    column_name_to_index = {}
    profile = []
    samples_list = []
    predictions_dict = {}
    reading_data = False
    got_column_indices = False

    for line in cami_content:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip("\n")

        # parse header with column indices
        if line.startswith("@@"):
            for index, column_name in enumerate(line[2:].split("\t")):
                column_name_to_index[column_name] = index
            (
                index_rank,
                index_taxid,
                index_percentage,
                index_taxpath,
                index_taxpathsn,
            ) = get_column_indices(column_name_to_index)
            got_column_indices = True
            reading_data = False
            continue

        # parse header with metadata
        if line.startswith("@"):
            # if last line contained sample data and new header starts, store profile for sample
            if reading_data:
                if "SAMPLEID" in header and "VERSION" in header and "RANKS" in header:
                    if len(profile) > 0:
                        samples_list.append((header["SAMPLEID"], header, profile))
                        profile = []
                        predictions_dict = {}
                else:
                    logger.error(
                        "Header is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS."
                    )
                    raise RuntimeError
                header = {}
            reading_data = False
            got_column_indices = False
            key, value = line[1:].split(":", 1)
            header[key.upper()] = value.strip()
            continue

        if not got_column_indices:
            logger.error(
                "Header line starting with @@ is missing or at wrong position."
            )
            raise RuntimeError

        reading_data = True
        row_data = line.split("\t")

        taxid = row_data[index_taxid]
        # if there is already a prediction for taxon, only sum abundance
        if taxid in predictions_dict:
            prediction = predictions_dict[taxid]
            prediction.percentage += float(row_data[index_percentage])
        else:
            if int(float(row_data[index_percentage])) == 0:
                continue
            prediction = Prediction()
            predictions_dict[taxid] = prediction
            prediction.taxid = row_data[index_taxid]
            prediction.rank = row_data[index_rank]
            prediction.percentage = float(row_data[index_percentage])
            prediction.taxpath = row_data[index_taxpath]
            if isinstance(index_taxpathsn, int):
                prediction.taxpathsn = row_data[index_taxpathsn]
            else:
                prediction.taxpathsn = None
            profile.append(prediction)

    # store profile for last sample
    if "SAMPLEID" in header and "VERSION" in header and "RANKS" in header:
        if reading_data and len(profile) > 0:
            samples_list.append((header["SAMPLEID"], header, profile))
    else:
        logger.error(
            "Header is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS."
        )
        raise RuntimeError

    return samples_list


def create_output_folder(outfolder):
    """
    Helper function that creates the output folder if it does not exist.
    :param outfolder: location of output folder
    :return: None
    """
    if not os.path.exists(outfolder):
        logger.info(f"Creating output folder: {outfolder}")
        os.makedirs(outfolder)


def check_download_args(args, db_type):
    """
    Helper function that checks if the input arguments are valid.
    :param args: input arguments
    :param db_type: type of database options: "pretrained" or "default"
    :return: None
    """
    if args.database not in ["genbank", "gtdb"]:
        logger.error(
            f"Invalid database: {args.database}. Now only support genbank and gtdb."
        )
        sys.exit(1)

    if args.k not in [21, 31, 51]:
        logger.error(f"Invalid k: {args.k}. Now only support 21, 31, and 51.")
        sys.exit(1)

    if args.database == "genbank":
        if args.ncbi_organism is None:
            logger.warning(
                "No NCBI organism specified using parameter --ncbi_organism. Using the default: bacteria"
            )
            args.ncbi_organism = "bacteria"

        if args.ncbi_organism not in [
            "archaea",
            "bacteria",
            "fungi",
            "virus",
            "protozoa",
        ]:
            logger.error(
                f"Invalid NCBI organism: {args.ncbi_organism}. Now only support archaea, bacteria, fungi, virus, and protozoa."
            )
            sys.exit(1)

        if db_type == "pretrained" and args.ncbi_organism == "virus":
            logger.error("We now haven't supported for virus database.")
            sys.exit(1)

def _temp_get_genome_name(sig_file_path, ksize):

    res = get_info_from_single_sig(sig_file_path, ksize)
    if res:
        return res[0]
    else:
        return None

def temp_generate_inputs(
    selected_genomes_file_path: str,
    sig_info_dict: Dict[str, Tuple[str, float, int, int]],
    ksize: int,
    num_threads: int = 16,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Temporary Helper function that generates the required input for `yacht run`.
    :param selected_genomes_file_path: Path to a file containing all the genome file path.
    :param num_threads: Number of threads to use for multiprocessing when reading the comparison files. Default is 16.
    :param sig_info_dict:
        A dictionary mapping each genome signature name to a tuple containing metadata: 
        (md5sum, minhash mean abundance, minhash hashes length, minhash scaled).
        - md5sum: Checksum for data integrity.
        - minhash mean abundance: The mean abundance for the genome's minhash.
        - minhash hashes length: The length of minhash hashes.
        - minhash scaled: The scaling factor for the minhash.    
    :return 
        manifest_df: a dataframe containing the processed reference signature information
    """
    # get info from the signature files of selected genomes
    selected_sig_files = pd.read_csv(selected_genomes_file_path, sep="\t", header=None)
    selected_sig_files = selected_sig_files[0].to_list()
    
    # get the genome name from the signature files using multiprocessing
    with Pool(num_threads) as p:
        result_list = p.starmap(_temp_get_genome_name, [(sig_file_path, ksize) for sig_file_path in selected_sig_files])
    selected_genome_names_set = set([x for x in result_list if x])

    # remove the close related organisms from the reference genome list
    manifest_df = []
    for sig_name, (
        md5sum,
        minhash_mean_abundance,
        minhash_hashes_len,
        minhash_scaled,
    ) in tqdm(sig_info_dict.items(), desc="Removing close related organisms from the reference genome list"):
        if sig_name in selected_genome_names_set:
            manifest_df.append(
                (
                    sig_name,
                    md5sum,
                    minhash_hashes_len,
                    get_num_kmers(
                        minhash_mean_abundance,
                        minhash_hashes_len,
                        minhash_scaled,
                        False,
                    ),
                    minhash_scaled,
                )
            )
    manifest_df = pd.DataFrame(
        manifest_df,
        columns=[
            "organism_name",
            "md5sum",
            "num_unique_kmers_in_genome_sketch",
            "num_total_kmers_in_genome_sketch",
            "genome_scale_factor",
        ],
    )

    return manifest_df