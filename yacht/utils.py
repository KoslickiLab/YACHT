import os, sys
import numpy as np
import sourmash
from tqdm import tqdm, trange
import pandas as pd
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
from loguru import logger
from typing import Optional, Union, List, Set, Dict, Tuple
import math

# Configure Loguru logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")

# Set up global variables
__version__ = '1.1.0'
GITHUB_API_URL = "https://api.github.com/repos/KoslickiLab/YACHT/contents/demo/{path}"
GITHUB_RAW_URL = "https://raw.githubusercontent.com/KoslickiLab/YACHT/main/demo/{path}"
BASE_URL = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/"
ZENODO_COMMUNITY_URL = "https://zenodo.org/api/records/?communities=yacht"

def load_signature_with_ksize(filename: str, ksize: int) -> sourmash.SourmashSignature:
    """
    Helper function that loads the signature for a given kmer size from the provided signature file.
    Filename should point to a .sig file. Raises exception if given kmer size is not present in the file.
    :param filename: string (location of the signature file with .sig.gz format)
    :param ksize: int (size of kmer)
    :return: sourmash signature
    """
    # Take the first sample signature with the given kmer size
    sketches = list(sourmash.load_file_as_signatures(filename, ksize=ksize))
    if len(sketches) != 1:
        raise ValueError(f"Expected exactly one signature with ksize {ksize} in {filename}, found {len(sketches)}")
    if len(sketches[0].minhash.hashes) == 0:
        raise ValueError(f"Empty sketch in signature. This may be due to too high of a scale factor, please reduce it, eg. --scaled=1, and try again.")
    if math.isnan(sketches[0].minhash.mean_abundance):
        raise ValueError(f"No mean abundance. This may be due to too high of a scale factor, please reduce it, eg. --scaled=1, and try again.")
    return sketches[0]

def get_num_kmers(minhash_mean_abundance: Optional[float], minhash_hashes_len: int, minhash_scaled: int, scale: bool = True) -> int:
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

def get_info_from_single_sig(sig_file: str, ksize: int) -> Tuple[str, str, float, int, int]:
    """
    Helper function that gets signature information (name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled) from a single sourmash signature file.
    :param sig_file: string (location of the signature file with .sig.gz format)
    :param ksize: int (size of kmer)
    :return: tuple (name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled)
    """
    sig = load_signature_with_ksize(sig_file, ksize)
    return (sig.name, sig.md5sum(), sig.minhash.mean_abundance, len(sig.minhash.hashes), sig.minhash.scaled)

def collect_signature_info(num_threads: int, ksize: int, path_to_temp_dir: str) -> Dict[str, Tuple[str, float, int, int]]:
    """
    Helper function that collects signature information (name, md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled) from a sourmash signature database.
    :param num_threads: int (number of threads to use)
    :param ksize: int (size of kmer)
    :param path_to_temp_dir: string (path to the folder to store the intermediate files)
    :return: a dictionary mapping signature name to a tuple (md5sum, minhash mean abundance, minhash_hashes_len, minhash scaled)
    """
    ## extract in parallel
    with Pool(num_threads) as p:
        signatures = p.starmap(get_info_from_single_sig, [(os.path.join(path_to_temp_dir, 'signatures', file), ksize) for file in os.listdir(os.path.join(path_to_temp_dir, 'signatures'))])
    
    return {sig[0]:(sig[1], sig[2], sig[3], sig[4]) for sig in tqdm(signatures)}

def run_multisearch(num_threads: int, ani_thresh: float, ksize: int, scale: int, path_to_temp_dir: str) -> Dict[str, List[str]]:
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
    sig_files = pd.DataFrame([os.path.join(path_to_temp_dir, 'signatures', file) for file in os.listdir(os.path.join(path_to_temp_dir, 'signatures'))])
    sig_files_path = os.path.join(path_to_temp_dir, 'training_sig_files.txt')
    sig_files.to_csv(sig_files_path, header=False, index=False)
    
    # convert ani threshold to containment threshold
    containment_thresh = (ani_thresh ** ksize)
    cmd = f"sourmash scripts multisearch {sig_files_path} {sig_files_path} -k {ksize} -s {scale} -c {num_threads} -t {containment_thresh} -o {os.path.join(path_to_temp_dir, 'training_multisearch_result.csv')}"
    logger.info(f"Running sourmash multisearch with command: {cmd}")
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise ValueError(f"Error running sourmash multisearch with command: {cmd}")
    
    # read the multisearch result
    multisearch_result = pd.read_csv(os.path.join(path_to_temp_dir, 'training_multisearch_result.csv'), sep=',', header=0)
    multisearch_result = multisearch_result.query('query_name != match_name').reset_index(drop=True)

    # because the multisearch result is not symmetric, that is
    # we have: A B score but not B A score
    # we need to make it symmetric
    A_TO_B = multisearch_result[['query_name','match_name']].drop_duplicates().reset_index(drop=True)
    B_TO_A = A_TO_B[['match_name','query_name']].rename(columns={'match_name':'query_name','query_name':'match_name'})
    multisearch_result = pd.concat([A_TO_B, B_TO_A]).drop_duplicates().reset_index(drop=True)
    
    # change column type to string
    multisearch_result['query_name'] = multisearch_result['query_name'].astype(str)
    multisearch_result['match_name'] = multisearch_result['match_name'].astype(str)
    
    return multisearch_result

def remove_corr_organisms_from_ref(sig_info_dict: Dict[str, Tuple[str, float, int, int]], multisearch_result: pd.DataFrame) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
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
    corr_organisms = sorted([str(query_name) for query_name in multisearch_result['query_name'].unique()])
    sizes = np.array([sig_info_dict[organism][2] for organism in corr_organisms])
    # sort organisms by size in ascending order, so we keep the largest organism, discard the smallest
    bysize = np.argsort(sizes)
    corr_organisms_bysize = np.array(corr_organisms)[bysize].tolist()
    
    # use dictionary to store the removed organisms and their close related organisms
    # key: removed organism name
    # value: a set of close related organisms
    mapping = multisearch_result.groupby('query_name')['match_name'].agg(set).to_dict()
    
    # remove the sorted organisms until all left genomes are distinct (e.g., ANI <= ani_thresh)
    temp_remove_set = set()
    # loop through the organisms size in ascending order
    for organism in tqdm(corr_organisms_bysize, desc='Removing close related organisms'):
        ## for a given organism check its close related organisms, see if there are any organisms left after removing those in the remove set
        ## if so, put this organism in the remove set
        left_corr_orgs = mapping[organism].difference(temp_remove_set)
        if len(left_corr_orgs) > 0:
            temp_remove_set.add(organism)

    # generate a dataframe with two columns: removed organism name and its close related organisms
    logger.info('Generating a dataframe with two columns: removed organism name and its close related organisms.')
    remove_corr_list = [(organism, ','.join(list(mapping[organism]))) for organism in tqdm(temp_remove_set)]
    remove_corr_df = pd.DataFrame(remove_corr_list, columns=['removed_org', 'corr_orgs'])
    
    # remove the close related organisms from the reference genome list
    manifest_df = []
    for sig_name, (md5sum, minhash_mean_abundance, minhash_hashes_len, minhash_scaled) in tqdm(sig_info_dict.items()):
        if sig_name not in temp_remove_set:
            manifest_df.append((sig_name, md5sum, minhash_hashes_len, get_num_kmers(minhash_mean_abundance, minhash_hashes_len, minhash_scaled, False), minhash_scaled))
    manifest_df = pd.DataFrame(manifest_df, columns=['organism_name', 'md5sum', 'num_unique_kmers_in_genome_sketch', 'num_total_kmers_in_genome_sketch', 'genome_scale_factor'])
    
    return remove_corr_df, manifest_df

class Prediction:
    """
    (thanks to https://github.com/CAMI-challenge/OPAL, this class is from its load_data.py)
    """
    
    def __init__(self):
        pass

    @property
    def rank(self):
        return self.__rank

    @property
    def taxid(self):
        return self.__taxid

    @property
    def percentage(self):
        return self.__percentage

    @property
    def taxpath(self):
        return self.__taxpath

    @property
    def taxpathsn(self):
        return self.__taxpathsn

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @taxid.setter
    def taxid(self, taxid):
        self.__taxid = taxid

    @percentage.setter
    def percentage(self, percentage):
        self.__percentage = percentage

    @taxpath.setter
    def taxpath(self, taxpath):
        self.__taxpath = taxpath

    @taxpathsn.setter
    def taxpathsn(self, taxpathsn):
        self.__taxpathsn = taxpathsn

    def get_dict(self):
        return self.__dict__

    def get_pretty_dict(self):
        return {property.split("_")[3]: value for property, value in self.__dict__.items()}

    def get_metadata(self):
        return {'rank': self.__rank, 'taxpath': self.__taxpath, 'taxpathsn': self.__taxpathsn}


def get_column_indices(column_name_to_index: Dict[str, int]) -> Tuple[int, int, int, int, Optional[int]]:
    """
    (thanks to https://github.com/CAMI-challenge/OPAL, this function is modified get_column_indices from its load_data.py)
    Helper function that gets the column indices for the following columns: TAXID, RANK, PERCENTAGE, TAXPATH, TAXPATHSN
    :param column_name_to_index: dictionary mapping column name to column index
    :return: indices for TAXID, RANK, PERCENTAGE, TAXPATH, TAXPATHSN
    """
    
    if "TAXID" not in column_name_to_index:
        logger.error("Column not found: {}".format("TAXID"))
        raise RuntimeError
    if "RANK" not in column_name_to_index:
        logger.error("Column not found: {}".format("RANK"))
        raise RuntimeError
    if "PERCENTAGE" not in column_name_to_index:
        logger.error("Column not found: {}".format("PERCENTAGE"))
        raise RuntimeError
    if "TAXPATH" not in column_name_to_index:
        logger.error("Column not found: {}".format("TAXPATH"))
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

def get_cami_profile(cami_content: List[str]) -> List[Tuple[str, Dict[str, str], List[Prediction]]]:
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
        line = line.rstrip('\n')

        # parse header with column indices
        if line.startswith("@@"):
            for index, column_name in enumerate(line[2:].split('\t')):
                column_name_to_index[column_name] = index
            index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(column_name_to_index)
            got_column_indices = True
            reading_data = False
            continue

        # parse header with metadata
        if line.startswith("@"):
            # if last line contained sample data and new header starts, store profile for sample
            if reading_data:
                if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
                    if len(profile) > 0:
                        samples_list.append((header['SAMPLEID'], header, profile))
                        profile = []
                        predictions_dict = {}
                else:
                    logger.error("Header is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.")
                    raise RuntimeError
                header = {}
            reading_data = False
            got_column_indices = False
            key, value = line[1:].split(':', 1)
            header[key.upper()] = value.strip()
            continue

        if not got_column_indices:
            logger.error("Header line starting with @@ is missing or at wrong position.")
            raise RuntimeError

        reading_data = True
        row_data = line.split('\t')

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
    if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
        if reading_data and len(profile) > 0:
            samples_list.append((header['SAMPLEID'], header, profile))
    else:
        logger.error("Header is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.")
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
        logger.error(f"Invalid database: {args.database}. Now only support genbank and gtdb.")
        os.exit(1)

    if args.k not in [21, 31, 51]:
        logger.error(f"Invalid k: {args.k}. Now only support 21, 31, and 51.")
        os.exit(1)

    if args.database == "genbank":
        if args.ncbi_organism is None:
            logger.warning("No NCBI organism specified using parameter --ncbi_organism. Using the default: bacteria")
            args.ncbi_organism = "bacteria"

        if args.ncbi_organism not in ["archaea", "bacteria", "fungi", "virus", "protozoa"]:
            logger.error(
                f"Invalid NCBI organism: {args.ncbi_organism}. Now only support archaea, bacteria, fungi, virus, and protozoa.")
            os.exit(1)

        if db_type == "pretrained":
            if args.ncbi_organism == "virus":
                logger.error("We now haven't supported for virus database.")
                os.exit(1)