import pandas as pd
from lifelines import CoxPHFitter

df = pd.read_csv(snakemake.input[0])

df = df.dropna(subset=["label_time", "label_event"])

id_cols = ["PATIENT_ID", "label_time", "label_event"]
features = [col for col in df.columns if col not in id_cols]

results = []

for gene in features:
    temp_df = df[[gene, "label_time", "label_event"]].dropna()
    if temp_df[gene].nunique() > 1:
        try:
            cph = CoxPHFitter()
            cph.fit(temp_df, duration_col="label_time", event_col="label_event")
            p = cph.summary.loc[gene, "p"]
            results.append((gene, p))
        except Exception as e:
            continue

top_genes = pd.DataFrame(results, columns=["gene", "p_value"])
top_genes = top_genes.sort_values(by="p_value").head(50)
top_genes.to_csv(snakemake.output[0], index=False)

#cox regression on all genes and select the top 50 genes with lowest p-value
#Cox regression models time-event data. 
#aka proportional hazards regression
#p-value is the stat significance for each gene