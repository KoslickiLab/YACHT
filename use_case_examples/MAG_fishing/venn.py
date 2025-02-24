from matplotlib import pyplot as plt
from matplotlib_venn import venn2  
import pandas as pd

#datasets for different ksizes

sample1 = 'result_SRR32008482.k51.xlsx'
df_sample1 = pd.read_excel(sample1, sheet_name='min_coverage0.5')
sample1_lst = df_sample1['organism_name']
sample1_set = {x for x in sample1_lst}

sample2 = 'result_SRR32008483.k51.xlsx'
df_sample2 = pd.read_excel(sample2, sheet_name='min_coverage0.1')
sample2_lst = df_sample2['organism_name']
sample2_set = {x for x in sample2_lst}

#venn diagram

venn2([sample1_set, sample2_set], alpha=0.5, set_labels=(f'SRR32008482\nmincov=0.5',f'SRR32008483\nmincov=0.1'))

plt.title("k-size=51, ANI=0.95")

# Save the figure
plt.savefig('venn.png',dpi=500)

