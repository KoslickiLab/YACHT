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
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")
    
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
    :return: a dictionary mapping signature name to a list of its close related genomes (ANI > ani_thresh)
    """
    results = {}
    
    # run the sourmash multisearch
    # save signature files to a text file
    sig_files = pd.DataFrame([os.path.join(path_to_temp_dir, 'signatures', file) for file in os.listdir(os.path.join(path_to_temp_dir, 'signatures'))])
    sig_files_path = os.path.join(path_to_temp_dir, 'training_sig_files.txt')
    sig_files.to_csv(sig_files_path, header=False, index=False)
    
    # convert ani threshold to containment threshold
    containment_thresh = 0.9*(ani_thresh ** ksize)
    cmd = f"sourmash scripts multisearch {sig_files_path} {sig_files_path} -k {ksize} -s {scale} -c {num_threads} -t {containment_thresh} -o {os.path.join(path_to_temp_dir, 'training_multisearch_result.csv')}"
    logger.info(f"Running sourmash multisearch with command: {cmd}")
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise ValueError(f"Error running sourmash multisearch with command: {cmd}")
    
    # read the multisearch result
    multisearch_result = pd.read_csv(os.path.join(path_to_temp_dir, 'training_multisearch_result.csv'), sep=',', header=0)
    multisearch_result = multisearch_result.drop_duplicates().reset_index(drop=True)
    multisearch_result = multisearch_result.query('query_name != match_name').reset_index(drop=True)
    
    for query_name, match_name in tqdm(multisearch_result[['query_name', 'match_name']].to_numpy()):
        if str(query_name) not in results:
            results[str(query_name)] = [str(match_name)]
        else:
            results[str(query_name)].append(str(match_name))
    
    return results

def remove_corr_organisms_from_ref(sig_info_dict: Dict[str, Tuple[str, float, int, int]], sig_same_genoms_dict: Dict[str, List[str]]) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
    """
    Helper function that removes the close related organisms from the reference matrix.
    :param sig_info_dict: a dictionary mapping all signature name from reference data to a tuple (md5sum, minhash mean abundance, minhash hashes length, minhash scaled)
    :param sig_same_genoms_dict: a dictionary mapping signature name to a list of its close related genomes (ANI > ani_thresh)
    :return 
        rep_remove_dict: a dictionary with key as representative signature name and value as a list of signatures to be removed
        manifest_df: a dataframe containing the processed reference signature information
    """
    # for each genome with close related genomes, pick the one with largest number of unique kmers
    rep_remove_dict = {}
    temp_remove_set = set()
    manifest_df = []
    for genome, same_genomes in tqdm(sig_same_genoms_dict.items()):
        # skip if the genome has been removed
        if genome in temp_remove_set:
            continue
        # keep same genome if it is not in the remove set
        same_genomes = list(set(same_genomes).difference(temp_remove_set))
        # get the number of unique kmers for each genome
        unique_kmers = np.array([sig_info_dict[genome][2]] + [sig_info_dict[same_genome][2] for same_genome in same_genomes])
        # get the index of the genome with largest number of unique kmers
        rep_idx = np.argmax(unique_kmers)
        # get the representative genome
        rep_genome = genome if rep_idx == 0 else same_genomes[rep_idx-1]
        # get the list of genomes to be removed
        remove_genomes = same_genomes if rep_idx == 0 else [genome] + same_genomes[:rep_idx-1] + same_genomes[rep_idx:]
        # update remove set
        temp_remove_set.update(remove_genomes)
        if len(remove_genomes) > 0:
            rep_remove_dict[rep_genome] = remove_genomes
    
    # remove the close related organisms from the reference genome list
    for sig_name, (md5sum, minhash_mean_abundance, minhash_hashes_len, minhash_scaled) in tqdm(sig_info_dict.items()):
        if sig_name not in temp_remove_set:
            manifest_df.append((sig_name, md5sum, minhash_hashes_len, get_num_kmers(minhash_mean_abundance, minhash_hashes_len, minhash_scaled, False), minhash_scaled))
    manifest_df = pd.DataFrame(manifest_df, columns=['organism_name', 'md5sum', 'num_unique_kmers_in_genome_sketch', 'num_total_kmers_in_genome_sketch', 'genome_scale_factor'])
    
    return rep_remove_dict, manifest_df
    
# def compute_sample_vector(sample_hashes, hash_to_idx):
#     """
#     Helper function that computes the sample vector for a given sample signature.
#     :param sample_hashes: hashes in the sample signature
#     :param hash_to_idx: dictionary mapping hashes to indices in the training dictionary
#     :return: numpy array (sample vector)
#     """
#     # total number of hashes in the training dictionary
#     hash_to_idx_keys = set(hash_to_idx.keys())
    
#     # total number of hashes in the sample
#     sample_hashes_keys = set(sample_hashes.keys())
    
#     # initialize the sample vector
#     sample_vector = np.zeros(len(hash_to_idx_keys))
    
#     # get the hashes that are in both the sample and the training dictionary
#     sample_intersect_training_hashes = hash_to_idx_keys.intersection(sample_hashes_keys)
    
#     # fill in the sample vector
#     for sh in tqdm(sample_intersect_training_hashes):
#         sample_vector[hash_to_idx[sh]] = sample_hashes[sh]

#     return sample_vector


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
            if float(row_data[index_percentage]) == .0:
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
