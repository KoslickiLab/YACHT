import numpy as np
import warnings
from scipy.stats import binom
from scipy.special import betaincinv
import pandas as pd
warnings.filterwarnings("ignore")


def get_nontrivial_idx(A, y):
    """
    This function takes the sparse matrix A and the vector y and returns the indices of the columns of A
    that have a non-zero inner product/overlap with y.
    :param A: sparse matrix
    :param y: vector
    :return: indices of columns of A that have a non-zero inner product with y
    """
    inners = A.T @ y
    nonz_idx = np.nonzero(inners)[0]
    return nonz_idx


def get_exclusive_indicators(A):
    """
    This function takes the sparse matrix A and returns a list of lists,
    where each of the ith list is the set of rows that are non-zero in the ith column.
    :param A: A sparse matrix. Should be binary, but doesn't have to be.
    :return: list(list(int))
    """
    # TODO: currently this isn't finding the unique rows, but rather the rows that are non-zero
    unique_locs = []
    m, N = A.shape
    # sum all the columns up
    col_sums = A.sum(axis=1)
    # look for the rows that have a 1 in them
    #unique_rows = np.nonzero(col_sums > 0)[0]  # FIXME: this is looking for rows that have an entry > 0, not == 1
    unique_rows = np.nonzero(col_sums == 1)[0]
    # turn this into a set
    unique_rows = set(unique_rows)
    # for each column, find the rows that are non-zero
    for i in range(N):
        non_zero_locs = np.nonzero(A[:, i])[0]
        # find the intersection of the two sets
        unique_in_col = list(unique_rows.intersection(non_zero_locs))
        # sort this, if need be
        # unique_in_col = sorted(unique_in_col)
        unique_locs.append(unique_in_col)
    return unique_locs


def get_alt_mut_rate(nu, thresh, ksize, significance=0.99):
    """
    Computes the alternative mutation rate for a given significance level. I.e. how much higher would the mutation rate
    have needed to be in order to have a false positive rate of significance (since we are setting the false negative
    rate to significance by design)?
    :param nu: Number of k-mers exclusive to the organism under consideration
    :param thresh: Number of exclusive k-mers I would need to observe in order to reject the null hypothesis (i.e.
    accept that the organism is present)
    :param ksize: K-mer size
    :param significance: value between 0 and 1 expressing the desired false positive rate (and by design, the false
    negative rate)
    :return: float (alternative mutation rate; how much higher would the mutation rate have needed to be in order to
    make FP and FN rates equal to significance)
    """
    # Replace binary search with the regularized incomplete Gamma function inverse: Solve[significance ==
    #   BetaRegularized[1 - (1 - mutCurr)^k, nu - thresh,
    #    1 + thresh], mutCurr]
    # per mathematica
    mut = 1 - (1 - betaincinv(nu - thresh, 1 + thresh, significance))**(1/ksize)
    if np.isnan(mut):
        return -1
    else:
        return mut


def single_hyp_test(
    y,
    unique_idx,
    ksize,
    significance=0.99,
    ani_thresh=0.95,
    min_coverage=1
):
    """
    Performs a single hypothesis test for the presence of a genome in a metagenome.
    :param y: vector of k-mer counts for the metagenome
    :param unique_idx: indices of the rows/hashes of A that are unique to the genome under consideration
    :param ksize: k-mer size
    :param significance: significance level for the hypothesis test
    :param ani_thresh: threshold for ANI (i.e. how similar do the genomes need to be in order to be considered the same)
    :param min_coverage: minimum coverage of the genome under consideration in the metagenome (float in [0, 1])
    :return: A whole bunch of stuff
    """
    # get the number of unique k-mers
    num_unique_kmers = len(unique_idx)
    # mutation rate
    non_mut_p = (ani_thresh)**ksize
    # assuming coverage of 1, how many unique k-mers would I need to observe in order to reject the null hypothesis?
    acceptance_threshold_wo_coverage = binom.ppf(1-significance, num_unique_kmers, non_mut_p)
    # what is the actual confidence of the test?
    actual_confidence_wo_coverage = 1-binom.cdf(acceptance_threshold_wo_coverage, num_unique_kmers, non_mut_p)
    # number of unique k-mers I would see given a coverage of min_coverage
    num_unique_kmers_coverage = int(num_unique_kmers * min_coverage)
    # how many unique k-mers would I need to observe in order to reject the null hypothesis,
    # assuming coverage of min_cov?
    acceptance_threshold_with_coverage = binom.ppf(1-significance, num_unique_kmers_coverage, non_mut_p)
    # what is the actual confidence of the test, assuming coverage of min_cov?
    actual_confidence_with_coverage = 1-binom.cdf(acceptance_threshold_with_coverage, num_unique_kmers_coverage, non_mut_p)
    # what is the alternative mutation rate? I.e. how much higher would the mutation rate (resp. how low of ANI)
    # have needed to be in order to have a false positive rate of significance
    # (since we are setting the false negative rate to significance by design)?
    alt_confidence_mut_rate = get_alt_mut_rate(num_unique_kmers, acceptance_threshold_wo_coverage, ksize, significance=significance)
    # same as above, but assuming coverage of min_cov
    alt_confidence_mut_rate_with_coverage = get_alt_mut_rate(num_unique_kmers_coverage, acceptance_threshold_with_coverage, ksize, significance=significance)

    # How many unique k-mers do I actually see?
    num_matches = len(np.nonzero(y[unique_idx])[0])
    p_val = binom.cdf(num_matches, num_unique_kmers, non_mut_p)
    # is the genome present? Takes coverage into account
    is_present = (num_matches >= acceptance_threshold_with_coverage)
    return is_present, p_val, num_unique_kmers, num_unique_kmers_coverage, num_matches, acceptance_threshold_wo_coverage, acceptance_threshold_with_coverage, actual_confidence_wo_coverage, actual_confidence_with_coverage, alt_confidence_mut_rate, alt_confidence_mut_rate_with_coverage


def hypothesis_recovery(
    A,
    y,
    ksize,
    significance=0.99,
    ani_thresh=0.95,
    min_coverage=1,
):
    """
    Go through each of the organisms that have non-zero overlap with the sample and perform a hypothesis test for the
    presence of that organism in the sample: have we seen enough k-mers exclusive to that organism to conclude that
    an organism with ANI > ani_thresh (to the one under consideration) is present in the sample?
    :param A: matrix of k-mer counts for the organisms in the training set
    :param y: vector of k-mer counts in the sample
    :param ksize: k-mer size
    :param significance: significance level for the hypothesis test
    :param ani_thresh: threshold for ANI (i.e. how similar do the genomes need to be in order to be considered the same)
    :param min_coverage: minimum coverage of the genome under consideration in the metagenome (float in [0, 1])
    :return: pandas dataframe with the results of the hypothesis tests
    """
    nont_idx = get_nontrivial_idx(A, y)
    N = np.shape(A)[1]
    A_sub = A[:, nont_idx]
    exclusive_indicators = get_exclusive_indicators(A_sub)
    nontriv_flags = np.zeros(N)
    nontriv_flags[nont_idx] = 1
    # Create a pandas dataframe to store the results
    results = pd.DataFrame(
        index=nont_idx,
        columns=[
            'is_present',
            'p_vals',
            'num_unique_kmers',
            'num_unique_kmers_coverage',
            'num_matches',
            'acceptance_threshold_wo_coverage',
            'acceptance_threshold_with_coverage',
            'actual_confidence_wo_coverage',
            'actual_confidence_with_coverage',
            'alt_confidence_mut_rate',
            'alt_confidence_mut_rate_with_coverage',
        ]
    )
    for i in range(len(nont_idx)):
        exclusive_idx = exclusive_indicators[i]
        curr_result = single_hyp_test(
            y,
            exclusive_idx,
            ksize,
            significance=significance,
            ani_thresh=ani_thresh,
            min_coverage=min_coverage,
        )
        curr_idx = nont_idx[i]
        results.loc[curr_idx] = curr_result
    return results, nontriv_flags