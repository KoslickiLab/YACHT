import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import math

def load_data(tab_file, excel_file, tab_col, excel_col, sheet_name=None):
    """Loads and returns organism names from both files, removing non-breaking spaces."""
    tab_df = pd.read_csv(tab_file, sep='\t')

    excel_df = pd.read_excel(excel_file, sheet_name=sheet_name) if sheet_name else pd.read_excel(excel_file)
    excel_df[excel_col] = excel_df[excel_col].astype(str).str.replace('\xa0', ' ', regex=True)

    return dict(zip(tab_df["Paper Identifier"], tab_df[tab_col])), set(excel_df[excel_col])

def find_matches(query_dict, ref_set):
    """Match organisms found in results from synthetic microorganism table"""
    keys_list = list(query_dict.keys())
    for key, query in query_dict.items():
        for organism in ref_set:
            if str(query) in organism:  # Check if query is found in ref
                keys_list = [item.replace(key, organism) for item in keys_list]
    return keys_list

def plot_venn(set1, set2, label1='Singer et al. (2016)', label2='YACHT', output_file='venn_diagram.pdf',title="Comparing YACHT to published synthetic metagenome"):
    """Plots a Venn diagram showing overlap and saves the figure."""
    venn2([set1, set2], set_labels=(label1, label2))
    plt.title(title)
    plt.savefig(output_file)
    plt.clf()

def main(tab_file, excel_file, title, tab_col='GTDB Accession', excel_col='organism_name', sheet_name=None, outfile="venn_diagram.pdf"):
    tab_dict, excel_set = load_data(tab_file, excel_file, tab_col, excel_col, sheet_name)
    tab_set = set(find_matches(tab_dict, excel_set))
    plot_venn(tab_set, excel_set, output_file=outfile,title=title)

# Example usage:
k=31
cov=0.05

ani="0.80"
main('matched_identifiers.tab', f'result_k31_ani0.80.xlsx', sheet_name=f'min_coverage{cov}', outfile=f"venn_k{k}_{ani}_min{cov}.png",title=f"k={k}, ANI={ani}, min_coverage={cov}")

ani="0.95"
main('matched_identifiers.tab', f'result_k31_ani0.95.xlsx', sheet_name=f'min_coverage{cov}', outfile=f"venn_k{k}_{ani}_min{cov}.png",title=f"k={k}, ANI={ani}, min_coverage={cov}")

ani="0.9995"
main('matched_identifiers.tab', f'result_k31_ani0.9995.xlsx', sheet_name=f'min_coverage{cov}', outfile=f"venn_k{k}_{ani}_min{cov}.png",title=f"k={k}, ANI={ani}, min_coverage={cov}")

ani="0.95"
cov=0.0001
main('matched_identifiers.tab', f'result_k31_ani0.95.xlsx', sheet_name=f'min_coverage{cov}', outfile=f"venn_k{k}_{ani}_min{cov}.png",title=f"k={k}, ANI={ani}, min_coverage={cov}")
cov=0.25
main('matched_identifiers.tab', f'result_k31_ani0.95.xlsx', sheet_name=f'min_coverage{cov}', outfile=f"venn_k{k}_{ani}_min{cov}.png",title=f"k={k}, ANI={ani}, min_coverage={cov}")
