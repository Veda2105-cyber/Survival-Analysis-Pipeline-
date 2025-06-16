import pandas as pd

top_genes_file = snakemake.input["top_genes"]
expression_file = snakemake.input["expression"]
methylation_file = snakemake.input["methylation"]
cna_file = snakemake.input["cna"]
output_file = snakemake.output[0]

top_genes = pd.read_csv(top_genes_file, sep="\t")["Hugo_Symbol"].tolist()

expr = pd.read_csv(expression_file, sep="\t")
meth = pd.read_csv(methylation_file, sep="\t")
cna = pd.read_csv(cna_file, sep="\t")

expr_filtered = expr[expr["Hugo_Symbol"].isin(top_genes)]
meth_filtered = meth[meth["Hugo_Symbol"].isin(top_genes)]
cna_filtered = cna[cna["Hugo_Symbol"].isin(top_genes)]

merged = expr_filtered.merge(meth_filtered, on="Hugo_Symbol", suffixes=("_expr", "_meth"))
merged = merged.merge(cna_filtered, on="Hugo_Symbol")

merged.to_csv(output_file, index=False)

#what is HUGO?
#HUGO: The Human Genome Organisation HUGO Gene Nomenclature Committee (HGNC) assigns names to genes
#this file is for merging exp, methylation, and cna data
#depends on get_top_mutated_genes.py and normalize_expression.py