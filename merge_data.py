import pandas as pd

with open(snakemake.input.clinical_sample, 'r') as f:
    header_line = f.readline().lstrip("#").strip().split("\t")

clinical_df = pd.read_csv(snakemake.input.clinical_sample, sep="\t", comment='#', names=header_line, skiprows=1)

with open(snakemake.input.mutations, 'r') as f:
    first_line = f.readline()
    if not first_line.startswith("#Sequenced_Samples:"):
        raise ValueError("wrong format in mut file")
    sample_ids = first_line.strip().split(":")[1].strip().split()

if "Sample Identifier" not in clinical_df.columns:
    raise ValueError("'Sample Identifier' column not found")

filtered_clinical = clinical_df[clinical_df["Sample Identifier"].isin(sample_ids)]

filtered_clinical.to_csv(snakemake.output[0], index=False)

#we homogenize the mutation and clincal data 
#so that we are working with something standard