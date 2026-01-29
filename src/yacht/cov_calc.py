#!/usr/bin/env python
import sourmash
import math
import numpy as np
import pandas as pd
from loguru import logger
from yacht.utils import ratio_lambda
from yacht.utils import mme_lambda
from yacht.utils import binary_search_lambda
from yacht.utils import mle_zip
from yacht.utils import bootstrap_interval
from yacht.utils import ani_from_lambda
from yacht.utils import _ContainArgs
from yacht.utils import AniResult
from yacht.utils import AdjustStatus, AdjustStatusType
from yacht.utils import SAMPLE_SIZE_CUTOFF, PVALUE_CUTOFF, MEDIAN_ANI_THRESHOLD, MAX_MEDIAN_FOR_MEAN_FINAL_EST, MIN_COUNT_THRESH, ksize
from scipy.stats import poisson, variation
from typing import Optional, Tuple, Dict, Any

no_adj = False #consider updating this in future SUPERYACHT arguments
winner_map = None #skipping this step in this version
kmers_lost_count = None

def cov_calc(sample_sig: sourmash.SourmashSignature, genome_sig: sourmash.SourmashSignature):
    """
    Function that calculates lambda according to Shaw and Yu (2024) from two sourmash.Minshash files (resresenting the sample and the genome sketches). 
    """

    myArgs = _ContainArgs()

    gn_hashes = genome_sig.minhash.hashes
    gn_kmers_keys = genome_sig.minhash.hashes.keys()
    gn_kmers_items = genome_sig.minhash.hashes.items()
    gn_dict = dict(gn_kmers_items)

    sample_hashes_keys = sample_sig.minhash.hashes.keys()
    samp_kmers_items = sample_sig.minhash.hashes.items()
    samp_dict = dict(samp_kmers_items)

    covs = []
    contain_count = 0
    for kmer in gn_hashes:
        if kmer in sample_hashes_keys:
            if samp_dict[kmer] == 0:
                continue
            contain_count += 1
            covs.append(samp_dict[kmer])

    if len(covs)==0:
        return None

    naive_ani = math.pow(contain_count/len(gn_kmers_items),
        1/ksize)

    # Caps naive_ani at 1.0 to prevent biologically impossible ANIs
    if naive_ani > 1.0:
        logger.debug(f"Naive ANI {naive_ani:.6f} exceeds 1.0, capping at 1.0")
        naive_ani = 1.0

    covs.sort()

    if len(covs) == 0: 
        covs.append(0) 
    
    len_ind = len(covs)//2
    median_cov = covs[len(covs)//2]

    pois_obj = poisson(median_cov) #creates a discrete frozen Poisson distribution object
    cov_max = float('inf')
    
    if median_cov < 30: #if median coverage of 30 is not fulfilled
        for i in range(len_ind,len(covs), 1):
            cov = covs[i]
            if pois_obj.cdf(cov) < PVALUE_CUTOFF:
                cov_max = cov
            else:
                break

        # Check if cov_max remains inf (i.e. no valid maximum found)
        if cov_max == float('inf'):
            logger.warning(
                f"Could not determine valid coverage maximum for genome {genome_sig.name}. "
                f"Median coverage: {median_cov}. Returning None."
            )
            return None

    full_covs = [0] * (len(gn_hashes) - contain_count)

    for cov in covs:
        if cov <= cov_max:
            full_covs.append(cov)
    var = variation(full_covs)
    if var is not None:
        logger.debug("VAR {} {}", var, genome_sig.name)

    mean_cov = sum(full_covs)//len(full_covs)
    geq1_mean_cov = sum(full_covs)//len(covs)
    if median_cov > MEDIAN_ANI_THRESHOLD:
        return_lambda = AdjustStatus.high()

    else:
        if (myArgs.ratio == True):
            test_lambda = ratio_lambda(full_covs, MIN_COUNT_THRESH)
        elif (myArgs.mme == True):
            test_lambda = mme_lambda(full_covs)
        elif (myArgs.bin == True):
            test_lambda = binary_search_lambda(full_covs)
        elif (myArgs.mle) == True:
            test_lambda = mle_zip(full_covs, gn_kmers_items)
        else:
            test_lambda = ratio_lambda(full_covs, MIN_COUNT_THRESH)

        if test_lambda is None:
            return_lambda = AdjustStatus.low()
        else:
            return_lambda = AdjustStatus.lambda_value(test_lambda)
            
    match return_lambda.status:

        case AdjustStatusType.LAMBDA:
            # executes if it is the Lambda case
            final_est_cov = return_lambda.value
            opt_lambda = final_est_cov

        case AdjustStatusType.HIGH:
            # executes if it is high coverage case
            if median_cov < MAX_MEDIAN_FOR_MEAN_FINAL_EST:
                final_est_cov = geq1_mean_cov
            else:
                final_est_cov = median_cov
            opt_lambda = final_est_cov

        case AdjustStatusType.LOW:
            # if it is the "low" case
            # final_est_cov logic is handled elsewhere, or use a default
            opt_lambda = None

        # Adding a "wild-card" case, just to be safe
        case _:
            opt_lambda = None

    opt_est_ani = ani_from_lambda(opt_lambda, mean_cov, 31, full_covs)

    if opt_lambda == None or opt_est_ani == None or no_adj == True:
        final_est_ani = naive_ani
    else:
        final_est_ani = opt_est_ani


    low_ani, high_ani, low_lambda, high_lambda= None, None, None, None

#Conditional calculation of confidence intervals

    if myArgs.ci_int==True and opt_lambda is not None:
        bootstrap = bootstrap_interval(full_covs, ksize, myArgs)
        low_ani = bootstrap[0]
        high_ani = bootstrap[1]
        low_lambda = bootstrap[2]
        high_lambda = bootstrap[3]

    if sample_sig.name:
        seq_name = sample_sig.name
    else:
        seq_name = sample_sig.filename

#This is code related to the winner_map situation
    #kmers_lost = kmers_lost_count if winner_map is not None else None 

    ani_result = AniResult(
        naive_ani=naive_ani,
        final_est_ani=final_est_ani,
        final_est_cov=opt_lambda,
        seq_name=seq_name,
        gn_name=genome_sig.filename,
        contig_name=genome_sig.name,
        mean_cov=geq1_mean_cov,
        median_cov=median_cov,
        containment_index=(contain_count, len(gn_hashes)),
        lambda_status=return_lambda,
        ani_ci=(low_ani, high_ani),
        lambda_ci=(low_lambda, high_lambda),
        genome_sketch=genome_sig,
        rel_abund=None,
        seq_abund=None,
        kmers_lost=None,
    )

    results = [ani_result]

    columns_ani = [
        "naive_ani",  # the ani according to naive calculations
        "final_est_ani",  # final estimated ani
        "final_est_cov",  # The final estimated coverage
        "seq_name",  # the name of the sequence
        "gn_name",  # the name of the genome
        "mean_cov",  # the mean coverage observed
        "median_cov", # the median coverage observed
        "containment_index", #the containment index
        "lambda_status", #lambda status
        "ani_ci", #ani confidence interval
        "lambda_ci", #lambda confidence interval
        "genome_sketch", #genome Sourmash signature file
        "rel_abund", #The taxonomic abundance observed (set to None for now)
        "seq_abund", #the number of sequenes that match a given genome (set to None for now)
        "kmers_lost" #the number of kmers that were reassigned during the tax triage step. (set to None for now)
    ]

    cov_calc_df = pd.DataFrame(results, columns=columns_ani)
    
    return cov_calc_df

            
        

        
