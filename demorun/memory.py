import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("memory_data2", sep="\s+", engine="python")


# Plotting
plt.figure(figsize=(8, 5))
plt.plot(df['num_GTDB_files'], df['original(kBytes)'] / 1024, label='Original', marker='o', color='skyblue')
plt.plot(df['num_GTDB_files'], df['modified(kBytes)'] / 1024, label='Modified', marker='x', color='orange')

plt.xlabel("Number of GTDB Genome Files")
plt.ylabel("Memory Usage (MB)")
plt.title("Memory Usage Scaling with Reference Database Size")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("memory_scaling.png", dpi=300)
# plt.show()
