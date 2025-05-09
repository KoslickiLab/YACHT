import os
import sys
import numpy as np
import sourmash
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool
from loguru import logger
from typing import Optional, List, Set, Dict, Tuple
import shutil
import gzip
from glob import glob

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)

# Set up contants
COL_NOT_FOUND_ERROR = "Column not found: {}"
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

# Set up global variables
__version__ = "1.3.0"
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
    Helper function that gets signature information (raw file path, name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled) from a single sourmash signature file.
    :param sig_file: string (location of the signature file with .sig.gz format)
    :param ksize: int (size of kmer)
    :return: tuple (name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled)
    """
    try:
        sig = load_signature_with_ksize(sig_file, ksize)
        return (
            sig_file,
            sig.name,
            sig.md5sum(),
            sig.minhash.mean_abundance,
            len(sig.minhash.hashes),
            sig.minhash.scaled,
        )
    except:
        logger.warning(f"CANNOT extract the relevant info from the signature file: {sig_file}")
        return None

def search_by_indexing(sample_sig_file_path, organism_sig_file_path, path_to_sample_temp_dir,path_to_genome_temp_dir):
    # cmd = f"{FILE_LOCATION}/run_yacht_train_core -t {num_threads} -c {containment_thresh} -p {passes} {sig_files_path} {path_to_temp_dir} {os.path.join(path_to_temp_dir, 'selected_result.tsv')}"
    cmd = f"/usr/bin/time -v {FILE_LOCATION}/run_compute_similarity {sample_sig_file_path} {organism_sig_file_path} {path_to_genome_temp_dir} -t 1 -c 0.9 -p 1 -C -o {os.path.join(path_to_sample_temp_dir, 'no_multisearch_result.csv')}> {path_to_sample_temp_dir}/cpp.log 2>&1"
    # logger.info(f"Running compute similarity algorithm with command: {cmd}")
    logger.info(f"Running compute similarity algorithm with command")
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise ValueError(f"Error running compute similarityn algorithm with command: {cmd}")
    
def run_yacht_train_core(
    num_threads: int, ani_thresh: float, ksize: int, path_to_temp_dir: str, sig_info_dict: Dict[str, Tuple[str, float, int, int, str]], num_genome_threshold: int = 1000000
) -> Dict[str, List[str]]:
    """
    Helper function that runs the cpp script developed by Mahmudur Rahman Hera to find the closely related genomes with ANI > ani_thresh from the reference database,
    then remove them, and generate a dataframe with the selected genomes.
    :param num_threads: int (number of threads to use)
    :param ani_thresh: float (threshold for ANI, below which we consider two organisms to be distinct)
    :param ksize: int (size of kmer)
    :param path_to_temp_dir: string (path to the folder to store the intermediate files)
    :return: a dataframe containing the selected reference signature information
    """

    # run Mahmudur's cpp for genome comparison
    # save signature files to a text file
    sig_files = pd.DataFrame(
        [
            os.path.join(path_to_temp_dir, "signatures", file)
            for file in os.listdir(os.path.join(path_to_temp_dir, "signatures"))
        ]
    )
    sig_files_path = os.path.join(path_to_temp_dir, "training_sig_files.tsv")
    sig_files.to_csv(sig_files_path, header=False, index=False)

    # convert ani threshold to containment threshold
    containment_thresh = ani_thresh**ksize
    total_sig_files = len(sig_files)
    if total_sig_files <= num_genome_threshold:
        passes = 1
    else:
        passes = int(total_sig_files / num_genome_threshold) + 1
    cmd = f"{FILE_LOCATION}/run_yacht_train_core -t {num_threads} -c {containment_thresh} -p {passes} {sig_files_path} {path_to_temp_dir} {os.path.join(path_to_temp_dir, 'selected_result.tsv')}"
    logger.info(f"Running comparison algorithm with command: {cmd}")
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise ValueError(f"Error running comparison algorithm with command: {cmd}")

    # move all split comparison files to a single foldr
    os.makedirs(os.path.join(path_to_temp_dir, "comparison_files"), exist_ok=True)
    for file in glob(os.path.join(path_to_temp_dir, "*.txt")):
        shutil.move(file, os.path.join(path_to_temp_dir, "comparison_files"))

    # get info from the signature files of selected genomes
    selected_sig_files = pd.read_csv(os.path.join(path_to_temp_dir, 'selected_result.tsv'), sep="\t", header=None)
    selected_sig_files = selected_sig_files[0].to_list()
    
    # get the mapping from signature file name to genome name
    mapping = {sig_info_dict[name][-1]:name for name in sig_info_dict}
    selected_genome_names_set = set([mapping[sig_file_path] for sig_file_path in selected_sig_files])

    # remove the close related organisms from the reference genome list
    manifest_df = []
    for sig_name, (
        md5sum,
        minhash_mean_abundance,
        minhash_hashes_len,
        minhash_scaled,
        _
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



def collect_signature_info(
    num_threads: int, ksize: int, path_to_temp_dir: str
) -> Dict[str, Tuple[str, float, int, int]]:
    """
    Helper function that collects signature information (raw file path, name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled) from a sourmash signature database.
    :param num_threads: int (number of threads to use)
    :param ksize: int (size of kmer)
    :param path_to_temp_dir: string (path to the folder to store the intermediate files)
    :return: a dictionary mapping signature name to a tuple (md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled, raw file path)
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

    return {sig[1]: (sig[2], sig[3], sig[4], sig[5], sig[0]) for sig in tqdm(signatures) if sig}


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


def _decompress_and_remove(file_path: str) -> None:
    """
    Decompresses a GZIP-compressed file and removes the original compressed file.
    :param file_path: The path to the .sig.gz file that needs to be decompressed and deleted.
    :return: None
    """
    try:
        output_filename = os.path.splitext(file_path)[0]
        with gzip.open(file_path, 'rb') as f_in:
            with open(output_filename, 'wb') as f_out:
                f_out.write(f_in.read())

        os.remove(file_path)

    except Exception as e:
        logger.info(f"Failed to process {file_path}: {e}")
        
def decompress_all_sig_files(sig_files: List[str], num_threads: int) -> None:
    """
    Decompresses all .sig.gz files in the list using multiple threads.
    :param sig_files: List of .sig.gz files that need to be decompressed.
    :param num_threads: Number of threads to use for decompression.
    :return: None
    """
    with Pool(num_threads) as p:
        p.map(_decompress_and_remove, sig_files)
        
    logger.info("All .sig.gz files have been decompressed.")