import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

input_file = snakemake.input[0]
output_file = snakemake.output[0]
suffix = snakemake.params.suffix
title = snakemake.params.title

df = pd.read_csv(input_file)
df.set_index("Hugo_Symbol", inplace=True)

selected = df[[col for col in df.columns if col.endswith(suffix)]]

selected = selected.dropna(how="all", axis=0)  # drop genes with all NaNs
selected = selected.dropna(how="all", axis=1)  # drop samples with all NaNs
selected = selected.fillna(0)                  # fill remaining NaNs

sns.clustermap(selected, cmap="vlag", figsize=(14, 10), xticklabels=False)
plt.title(title)

os.makedirs(os.path.dirname(output_file), exist_ok=True)
plt.savefig(output_file)

#NaN errors are the worst, we know now.
#this is called thrice depending on mere_top_gene_multiomics