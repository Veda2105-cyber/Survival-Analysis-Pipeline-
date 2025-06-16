import pandas as pd
import numpy as np

input_file = snakemake.input[0]
output_file = snakemake.output[0]

df = pd.read_csv(input_file, sep='\t', index_col=0)

df_log = np.log2(df + 1)

df_log.to_csv(output_file, sep='\t')

#log2 transorm, we had to read quite a bit on it
#when it was suggested by AI