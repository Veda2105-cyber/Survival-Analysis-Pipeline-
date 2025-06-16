import pandas as pd
import os


input_files = snakemake.input
output_file = snakemake.output[0]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, 'w') as f:
    for file in input_files:
        try:
            df = pd.read_csv(file, sep="\t", low_memory=False)
            f.write(f"{file}: {df.shape[0]} rows, {df.shape[1]} cols\n")
        except Exception as e:
            f.write(f"{file}: ERROR - {str(e)}\n")

# summerises dimensions of inputs, text output for ref. this is just for
# the sake of knowing what data we are working with atp
