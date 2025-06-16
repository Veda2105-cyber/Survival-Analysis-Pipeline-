import pandas as pd

input_file = snakemake.input[0]
output_file = snakemake.output[0]
top_n = int(snakemake.params.n)

df = pd.read_csv(input_file, sep="\t", comment="#")

top_genes = (
    df["Hugo_Symbol"]
    .value_counts()
    .head(top_n)
    .reset_index()
)
top_genes.columns = ["Hugo_Symbol", "Mutation_Count"]

top_genes.to_csv(output_file, sep="\t", index=False)

#removes mutation data that are least mutated, and uses top mutated genes
#no visbile output. nexttt.
