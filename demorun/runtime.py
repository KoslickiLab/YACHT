import pandas as pd
import matplotlib.pyplot as plt

# Load the runtime data (assuming tab-separated values or CSV with fixed format)
df = pd.read_csv("runtime_data", sep="\s+", engine='python')  # handles space-separated columns

# Convert "MM:SS" time format to total minutes as float
def time_to_minutes(t):
    minutes, seconds = map(int, t.split(":"))
    return minutes + seconds / 60.0

df['original_runtime'] = df['original(minutes)'].apply(time_to_minutes)
df['modified_runtime'] = df['modified(minutes)'].apply(time_to_minutes)

# Plotting
plt.figure(figsize=(8, 5))
x = df['CAMI_sample']
bar_width = 0.35

plt.bar(x - 0.175, df['original_runtime'], width=bar_width, label='Original', color='skyblue')
plt.bar(x + 0.175, df['modified_runtime'], width=bar_width, label='Optimized', color='orange')

plt.xlabel('CAMI Sample')
plt.ylabel('Runtime (minutes)')
plt.title('Runtime Comparison per CAMI Sample')
plt.xticks(df['CAMI_sample'])
plt.legend()
plt.tight_layout()
plt.savefig("runtime_comparison.png", dpi=300)
plt.show()
