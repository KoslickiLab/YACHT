from matplotlib import pyplot as plt
from matplotlib_venn import venn3
import pandas as pd

sheet='min_coverage0.05'

#datasets for different ksizes

file_ani95_k21 = 'result_k21_ani0.95.xlsx'
df_ani95_k21 = pd.read_excel(file_ani95_k21, sheet_name=sheet)
species_ani95_k21_lst = df_ani95_k21['organism_name']
species_ani95_k21_set = {x for x in species_ani95_k21_lst}

file_ani95_k31 = 'result_k31_ani0.95.xlsx'
df_ani95_k31 = pd.read_excel(file_ani95_k31, sheet_name=sheet)
species_ani95_k31_lst = df_ani95_k31['organism_name']
species_ani95_k31_set = {x for x in species_ani95_k31_lst}

file_ani95_k51 = 'result_k51_ani0.95.xlsx'
df_ani95_k51 = pd.read_excel(file_ani95_k51, sheet_name=sheet)
species_ani95_k51_lst = df_ani95_k51['organism_name']
species_ani95_k51_set = {x for x in species_ani95_k51_lst}

venn3([species_ani95_k21_set, species_ani95_k31_set, species_ani95_k51_set], set_labels=('k21', 'k31', 'k51'), alpha=0.5)
plt.title("Different k-sizes at ANI threshold 0.95 at minimal coverage of 0.05")

# Save the figure
plt.savefig('venn_low_abundance_species_ksize.png',dpi=500)


