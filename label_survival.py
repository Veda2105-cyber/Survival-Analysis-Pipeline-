import pandas as pd
import os

input_file = snakemake.input[0]
output_file = snakemake.output[0]

df = pd.read_csv(input_file, sep="\t|,", engine="python", skiprows=4)

df.columns = [col.strip() for col in df.columns]

df["label_event"] = df["OS_STATUS"].apply(lambda x: 1 if "DECEASED" in str(x).upper() else 0)
df["label_time"] = df["OS_MONTHS"]

out_df = df[["PATIENT_ID", "label_time", "label_event"]]

os.makedirs(os.path.dirname(output_file), exist_ok=True)
out_df.to_csv(output_file, index=False)

#we told you that we couldnt deliver a perfect survival analysis
#but we did.
#theres loading, cleaning, and standardizing of the data