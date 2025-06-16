import pandas as pd
import os

omics_file = snakemake.input["omics"]
labels_file = snakemake.input["labels"]
output_file = snakemake.output[0]

omics_df = pd.read_csv(omics_file)
labels_df = pd.read_csv(labels_file)

omics_df = omics_df.set_index("Hugo_Symbol").T
omics_df.index = omics_df.index.str.replace("_expr", "", regex=False)
omics_df.index = omics_df.index.str.replace("_x", "", regex=False)
omics_df.index = omics_df.index.str.replace("_meth", "", regex=False)
omics_df.index.name = "PATIENT_ID"
omics_df = omics_df.reset_index()

merged = pd.merge(labels_df, omics_df, on="PATIENT_ID", how="inner")

os.makedirs(os.path.dirname(output_file), exist_ok=True)
merged.to_csv(output_file, index=False)

#survival labels meet multi-omics data
#heavily reliant on the prev 2 scripts to work properly