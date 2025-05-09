import pandas as pd
import matplotlib.pyplot as plt

# Load both result files
df1 = pd.read_excel("orig.xlsx")
df2 = pd.read_excel("mod.xlsx")

# Truncate organism names at the first space
df1['short_name'] = df1['organism_name'].str.split().str[0]
df2['short_name'] = df2['organism_name'].str.split().str[0]

# Plot both curves
plt.figure(figsize=(12, 6))
plt.plot(df1['short_name'], df1['num_matches'], label='Original', marker='o')
plt.plot(df2['short_name'], df2['num_matches'], label='Optimized', marker='x', linestyle='--')

plt.xticks(rotation=90)
plt.xlabel("Organism")
plt.ylabel("Number of Matched k-mers in Sample")
plt.title("Comparison of Matched k-mers (Original vs. Optimized)")
plt.legend()
plt.tight_layout()
plt.savefig("result_comparison_plot_sample2.png", dpi=300, bbox_inches='tight')
