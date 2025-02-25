from matplotlib import pyplot as plt
from matplotlib_venn import venn3
import pandas as pd

sheet='min_coverage0.05'

#datasets for different ani

file_ani90 = 'result_k31_ani0.90.xlsx'
df_ani90 = pd.read_excel(file_ani90, sheet_name=sheet)
species_ani90_lst = df_ani90['organism_name']
species_ani90_set = {x for x in species_ani90_lst}

file_ani95 = 'result_k31_ani0.95.xlsx'
df_ani95 = pd.read_excel(file_ani95, sheet_name=sheet)
species_ani95_lst = df_ani95['organism_name']
species_ani95_set = {x for x in species_ani95_lst}

file_ani9995 = 'result_k31_ani0.9995.xlsx'
df_ani9995 = pd.read_excel(file_ani9995, sheet_name=sheet)
species_ani9995_lst = df_ani9995['organism_name']
species_ani9995_set = {x for x in species_ani9995_lst}

# Create subplots

venn3([species_ani90_set, species_ani95_set, species_ani9995_set], set_labels=('0.90', '0.95', '0.9995'), alpha=0.5)
plt.title("Different ANI thresholds at k-size=31 using minimal coverage 0.05")

# Save the figure
plt.savefig('venn_low_abundance_species_ani.png',dpi=500)


