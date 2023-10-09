import os, sys
import numpy as np
import sourmash
from tqdm import tqdm, trange
import scipy as sp 
import csv
import sqlite3
import datetime
import zipfile
from loguru import logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")
sig_info_names = ['name', 'minhash_mean_abundance', 'minhash_hashes_len', 'minhash_scaled']
    
def load_signature_with_ksize(filename, ksize):
    """
    Helper function that loads the signature for a given kmer size from the provided signature file.
    Filename should point to a .sig file. Raises exception if given kmer size is not present in the file.
    :param filename: string (location of the signature file)
    :param ksize: kmer size
    :return: sourmash signature
    """
    # Take the first sample signature with the given kmer size
    sketches = list(sourmash.load_file_as_signatures(filename, ksize=ksize))
    if len(sketches) != 1:
        raise ValueError(f"Expected exactly one signature with ksize {ksize} in {filename}, found {len(sketches)}")
    return sketches[0] if sketches else 0


def get_num_kmers(minhash_mean_abundance, minhash_hashes_len, minhash_scaled, scale=True):
    """
    Helper function that estimates the total number of kmers in a given sample.
    :param minhash_mean_abundance: float or None (mean abundance of the signature)
    :return: int (estimated total number of kmers)
    """
    # Abundances may not have been kept, in which case, just use 1
    if minhash_mean_abundance:
        num_kmers = minhash_mean_abundance * minhash_hashes_len
    else:
        num_kmers = minhash_hashes_len
    if scale:
        num_kmers *= minhash_scaled
    return np.round(num_kmers)

def get_num_kmers(sig_info, scale=True):
    """
    Helper function that estimates the total number of kmers in a given sample.
    :param sig_info: tuple (name, minhash mean abundance, minhash hashes length, minhash scaled)
    :return: int (estimated total number of kmers)
    """
    # Abundances may not have been kept, in which case, just use 1
    if sig_info[sig_info_names.index('minhash_mean_abundance')]:
        num_kmers = sig_info[sig_info_names.index('minhash_mean_abundance')] * sig_info[sig_info_names.index('minhash_hashes_len')]
    else:
        num_kmers = sig_info[sig_info_names.index('minhash_hashes_len')]
    if scale:
        num_kmers *= sig_info[sig_info_names.index('minhash_scaled')]
    return np.round(num_kmers)


def check_file_existence(file_path, error_description):
    """
    Helper function that checks if a file exists. If not, raises a ValueError with the given error description.
    :param file_path: string (location of the file)
    :param error_description: string (description of the error)
    :return: None
    """
    if not os.path.exists(file_path):
        raise ValueError(error_description)


def compute_sample_vector(sample_hashes, hash_to_idx):
    """
    Helper function that computes the sample vector for a given sample signature.
    :param sample_hashes: hashes in the sample signature
    :param hash_to_idx: dictionary mapping hashes to indices in the training dictionary
    :return: numpy array (sample vector)
    """
    # total number of hashes in the training dictionary
    hash_to_idx_keys = set(hash_to_idx.keys())
    
    # total number of hashes in the sample
    sample_hashes_keys = set(sample_hashes.keys())
    
    # initialize the sample vector
    sample_vector = np.zeros(len(hash_to_idx_keys))
    
    # get the hashes that are in both the sample and the training dictionary
    sample_intersect_training_hashes = hash_to_idx_keys.intersection(sample_hashes_keys)
    
    # fill in the sample vector
    for sh in tqdm(sample_intersect_training_hashes):
        sample_vector[hash_to_idx[sh]] = sample_hashes[sh]

    return sample_vector

def signatures_to_ref_matrix(signatures, ksize, signature_count, hash_to_idx_db_path):
    """
    Given signature generator, return a sparse matrix with one column per signature and one row per hash
    (union of the hashes)
    :param signatures: sourmash signatures obtained via sourmash.load_file_as_signatures(<file name>)
    :param ksize: kmer size
    :param signature_count: number of signatures in the sourmash signature file
    :param hash_to_idx_db_path: string (location of the hash_to_col_idx.sqlite file)
    :return:
    signature_info_list: list of sourmash signature information tuple (name, minhash mean abundance, minhash hashes length, minhash scaled)
    ref_matrix: sparse matrix with one column per signature and one row per hash (union of the hashes)
    hash_to_idx: dictionary mapping hash to row index in ref_matrix
    is_mismatch: False (if all signatures have the same kmer size) or True (the first signature with a different kmer size)
    """
    row_idx = []
    col_idx = []
    sig_values = []
    signature_info_list = []
    is_mismatch = False

    # Use a dictionary-like database to store hash to index mapping
    dirname, dbname = os.path.split(hash_to_idx_db_path)
    hash_to_idx = DatabaseDict(db_name=dbname, out_dir=dirname)
    next_idx = 0  # Next available index for a new hash

    # Iterate over all signatures
    for col, sig in enumerate(tqdm(signatures, total=signature_count)):
        
        # covert to sourmash list
        signature_info_list.append((sig.name, sig.minhash.mean_abundance, len(sig.minhash.hashes), sig.minhash.scaled))
        
        # check that all signatures have the same ksize as the one provided
        if sig.minhash.ksize != ksize:
            is_mismatch = True
            return signature_info_list, None, None, is_mismatch
        
        sig_hashes = sig.minhash.hashes
        for hash, count in sig_hashes.items():
            # Get the index for this hash， if new hash, and it and increment next_idx
            idx = hash_to_idx.setdefault(hash, next_idx)
            if idx == next_idx:  # New hash was added
                next_idx += 1

            # Append row, col, and value information for creating the sparse matrix
            row_idx.append(idx)
            col_idx.append(col)
            sig_values.append(count)

    # Create the sparse matrix
    ref_matrix = sp.sparse.csc_matrix((sig_values, (row_idx, col_idx)), shape=(next_idx, len(signature_info_list)))

    return signature_info_list, ref_matrix, hash_to_idx, is_mismatch

def count_files_in_zip(zip_path):
    """
    Helper function that counts the number of files in a zip file.
    """
    with zipfile.ZipFile(zip_path, 'r') as z:
        return len(z.namelist())

def get_uncorr_ref(ref_matrix, ksize, ani_thresh):
    """
    Given a reference matrix, return a new reference matrix with only uncorrelated organisms
    :param ref_matrix: sparse matrix with one column per signature and one row per hash (union of the hashes)
    :param ksize: int, size of kmer
    :param ani_thresh: threshold for mutation rate, below which we consider two organisms to be correlated/the same
    :return:
    binary_ref: a new reference matrix with only uncorrelated organisms in binary form (discarding counts)
    uncorr_idx: the indices of the organisms in the reference matrix that are uncorrelated/distinct
    """

    N = ref_matrix.shape[1]  # number of organisms
    immut_prob = ani_thresh ** ksize  # probability of immutation

    # Convert the input matrix to a binary matrix
    # since I don't think we are actually using the counts anywhere, we could probably get a performance
    # boost by just using a binary matrix to begin with
    binary_ref = (ref_matrix > 0).astype(int)
    # number of hashes in each organism
    sizes = np.array(np.sum(binary_ref, axis=0)).reshape(N)

    # sort organisms by size  in ascending order, so we keep the largest organism, discard the smallest
    bysize = np.argsort(sizes)
    # binary_ref sorted by size
    binary_ref_bysize = binary_ref[:, bysize]

    # Compute all pairwise intersections
    intersections = binary_ref_bysize.T.dot(binary_ref_bysize)
    # set diagonal to 0, since we don't want to compare an organism to itself
    intersections.setdiag(0)

    uncorr_idx_bysize = np.arange(N)
    immut_prob_sizes = immut_prob * sizes[bysize]

    for i in trange(N):
        # Remove organisms if they are too similar
        # (Note that we remove organism if there is at least one other organism with intersection > immut_prob_sizes[i])
        if np.max(intersections[i, uncorr_idx_bysize]) > immut_prob_sizes[i]:
            uncorr_idx_bysize = np.setdiff1d(uncorr_idx_bysize, i)

    # Sort the remaining indices, uncorr_idx is now the indices of the organisms in the reference matrix that are uncorrelated
    uncorr_idx = np.sort(bysize[uncorr_idx_bysize])

    return binary_ref[:, uncorr_idx], uncorr_idx


def write_processed_indices(filename, signature_info_list, uncorr_org_idx):
    """
    Write a csv file with the following columns: organism_name, original_index, processed_index,
    num_unique_kmers_in_genome_sketch, num_total_kmers_in_genome_sketch, genome_scale_factor.
    :param filename: output filename
    :param signatures: sourmash signatures
    :param uncorr_org_idx: the indices of the organisms in the reference matrix that are uncorrelated
    (via get_uncorr_ref)
    :return: None
    """
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['organism_name', 'original_index', 'processed_index', 'num_unique_kmers_in_genome_sketch', 'num_total_kmers_in_genome_sketch', 'genome_scale_factor'])
        for i, idx in enumerate(tqdm(uncorr_org_idx)):
            writer.writerow([signature_info_list[idx][sig_info_names.index('name')], idx, i, signature_info_list[idx][sig_info_names.index('minhash_hashes_len')], get_num_kmers(signature_info_list[idx], scale=False), signature_info_list[idx][sig_info_names.index('minhash_scaled')]])

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


def get_column_indices(column_name_to_index):
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

def get_cami_profile(cami_content):
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



class DatabaseDict:
    """
    Creaet a dictionary-like class that stores key-value pairs in a sqlite database in order to save RAM and speed up loading.
    """
    def __init__(self, db_name=None, out_dir=os.getcwd()):
        """
        Initialize the database dictionary.
        param: db_name: string (name of the database)
        param: out_dir: string (location of the database)
        """
        if db_name is None:
            db_name = self.__generate_temp_db_name()
        self.db_path = os.path.join(out_dir, db_name)
        self.conn = sqlite3.connect(self.db_path)
        self.c = self.conn.cursor()
        self.create_table()

    def __generate_temp_db_name(self):
        """
        Generate a temporary database name based on the current time.
        """
        now = datetime.datetime.now()
        return now.strftime("temp_db_%Y%m%d_%H%M%S.db")

    def delete_database(self):
        """
        Delete the database and close the connection.
        """
        try:
            self.conn.close()
            os.remove(self.db_path)
            print(f"Clean the database {self.db_path} done")
        except Exception as e:
            print(f"Error deleting the database: {e}")

    def create_table(self):
        """
        Create a hash_index table in the database if it does not exist.
        """
        self.c.execute("CREATE TABLE IF NOT EXISTS hash_index (hash_key TEXT PRIMARY KEY, index_val INTEGER)")
        self.conn.commit()

    def load_database(self, db_path):
        """
        Load a database into the current connection.
        param: db_path: string (location of the database)
        """
        self.conn.close()
        self.db_path = db_path
        self.conn = sqlite3.connect(self.db_path)
        self.c = self.conn.cursor()
        self.create_table()

    def write_from_dict(self, data_dict):
        """
        Write key-value pairs from an existing dictionary into the database.
        param: data_dict: dictionary (key-value pairs to be written into the database)
        """
        for key, value in tqdm(data_dict.items()):
            self.__setitem__(key, value)

    def setdefault(self, key, default=None):
        """
        Simulate the setdefault method of a dictionary to set a default value for a key if the key does not exist.
        param: key: string (key)
        param: default: any (default value)
        """
        existing_val = self.c.execute("SELECT index_val FROM hash_index WHERE hash_key=?", (key,)).fetchone()
        if existing_val is None:
            self.c.execute("INSERT INTO hash_index (hash_key, index_val) VALUES (?, ?)", (key, default))
            self.conn.commit()
            return default
        else:
            return existing_val[0]

    def get(self, key, default=None):
        """
        Get the value for a key. If the key does not exist, return the default value.
        param: key: string (key)
        param: default: any (default value)
        """
        
        val = self.c.execute("SELECT index_val FROM hash_index WHERE hash_key=?", (key,)).fetchone()
        if val is None:
            return default
        else:
            return val[0]

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.c.execute("REPLACE INTO hash_index (hash_key, index_val) VALUES (?, ?)", (key, value))
        self.conn.commit()

    def close(self):
        self.conn.close()