import pandas as pd

input_file = snakemake.input[0]
output_file = snakemake.output[0]
top_n = int(snakemake.params.top_n)

df = pd.read_csv(input_file, sep="\t", comment="#")

if 'Hugo_Symbol' not in df.columns:
    raise ValueError("Missing 'Hugo_Symbol' column in input.")

gene_counts = df['Hugo_Symbol'].value_counts().head(top_n)
filtered_df = df[df['Hugo_Symbol'].isin(gene_counts.index)]

filtered_df.to_csv(output_file, sep="\t", index=False)

#not quite sure why we put a summary step but it made getting
#the output muchhh easier