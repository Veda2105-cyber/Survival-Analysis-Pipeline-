import pandas as pd
from lifelines import CoxPHFitter

df = pd.read_csv("results/final_model_data.csv").dropna()

if 'PATIENT_ID' in df.columns:
    df = df.drop(columns=['PATIENT_ID'])

cph = CoxPHFitter()
cph.fit(df, duration_col="label_time", event_col="label_event")

cph.print_summary()

cph.summary.to_csv("results/cox_model_summary.csv")

significant = cph.summary[cph.summary["p"] < 0.05][["p"]]
significant.to_csv("results/top_cox_features.csv")

#finallyyy we train the model
#depends on merge w labels file