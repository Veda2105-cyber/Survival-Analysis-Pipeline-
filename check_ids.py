import pandas as pd

df1 = pd.read_csv("results/merged_top_gene_multiomics.csv")
print("\n omics matrix samples:")
print(df1.columns[:10].tolist())

df2 = pd.read_csv("results/labeled_survival_data.csv")
print("\n labelled survival samples:")
print(df2["PATIENT_ID"].head())

#comparing omics data w survival for consistency