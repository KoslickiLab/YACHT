import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn


def load_single_df(folder):
    args = pd.read_csv(folder + '/args.csv')
    results = pd.read_csv(folder + '/results.csv')
    results[args.columns] = pd.concat([args]*len(results.index), ignore_index=True)
    results['s'] = results['s_known']+results['s_unknown']
    results['pct_unknown'] = results['s_unknown']/results['s']
    results['fp_rate'] = results['false_positives']/results['s_unknown']
    results['fn_rate'] = results['false_negatives']/results['s_known']
    return results


def load_all(folder, default_rate = None):
    all_results = []
    result_folders = [f for f in os.listdir(folder) if f[:7] == 'results']
    result_folders.sort()
    for sf in result_folders:
        if os.path.exists(folder + '/' + sf + '/results.csv'):
            all_results.append(load_single_df(folder + '/' + sf))
    all_results_df = pd.concat(all_results, ignore_index=True)
    if default_rate is not None:
        all_results_df['high_fp_mut'] = all_results_df['high_fp_mut'].fillna(default_rate)
        all_results_df['low_fn_mut'] = all_results_df['low_fn_mut'].fillna(default_rate)
    return all_results_df


def gen_plots(data, xcol, xlabel=None, savefolder=None, file_prefix=''):
    plt.figure()
    fp_fn_ax = seaborn.lineplot(data = data, x = xcol, y = 'false_positives')
    seaborn.lineplot(data = data, x = xcol, y = 'false_negatives', ax=fp_fn_ax)
    plt.legend(['False Positives', '__nolegend__', 'False Negatives', '__nolegend__'])
    plt.ylabel('Count (Avg)')
    if xlabel is not None:
        plt.xlabel(xlabel)
    if savefolder:
        plt.savefig(savefolder + '/' + file_prefix + '_fpfn_plot.png')
        
    plt.figure()
    fpfn_rt_ax = seaborn.lineplot(data = data, x = xcol, y = 'fp_rate')
    seaborn.lineplot(data = data, x = xcol, y = 'fn_rate', ax=fpfn_rt_ax)
    plt.legend(['False Positive Rate', '__nolegend__', 'False Negative Rate', '__nolegend__'])
    plt.ylabel('Rate')
    if xlabel is not None:
        plt.xlabel(xlabel)
    if savefolder:
        plt.savefig(savefolder + '/' + file_prefix + '_fpfn_rate_plot.png')
        
    plt.figure()
    mut_rt_ax = seaborn.lineplot(data = data, x = xcol, y = 'high_fp_mut')
    seaborn.lineplot(data = data, x = xcol, y = 'low_fn_mut', ax=mut_rt_ax)
    plt.legend(['High FP Mut Rate', '__nolegend__', 'Low FN Mut Rate', '__nolegend__'])
    plt.ylabel('Mutation Rate')
    if xlabel is not None:
        plt.xlabel(xlabel)
    if savefolder:
        plt.savefig(savefolder + '/' + file_prefix + '_mut_rt_plot.png')
        
    return fp_fn_ax, fpfn_rt_ax, mut_rt_ax


def gen_time_mem_plots(data, xcol, xlabel=None, savefolder=None, file_prefix=''):
    plt.figure()
    fp_fn_ax = seaborn.lineplot(data = data, x = xcol, y = 'recovery_time')
    plt.ylabel('Avg Runtime (s)')
    if xlabel is not None:
        plt.xlabel(xlabel)
    if savefolder:
        plt.savefig(savefolder + '/' + file_prefix + '_runtime.png')
        
    plt.figure()
    fpfn_rt_ax = seaborn.lineplot(data = data, x = xcol, y = 'max_memory_usage')
    plt.ylabel('Peak Memory Usage (MiB)')
    if xlabel is not None:
        plt.xlabel(xlabel)
    if savefolder:
        plt.savefig(savefolder + '/' + file_prefix + '_memory.png')
        
    return fp_fn_ax, fpfn_rt_ax
               
               