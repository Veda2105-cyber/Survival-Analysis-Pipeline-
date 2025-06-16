import pandas as pd

df = pd.read_csv("results/final_model_data.csv")

print("dimensions of final matrix:", df.shape)
print("\n columns:", df.columns[:10].tolist(), "...")

print("\n event lable counts:\n", df["label_event"].value_counts())

print("\n null values/col (10):")
print(df.isnull().sum().sort_values(ascending=False).head(10))

print("\n feature variance (10):")
print(df.drop(columns=["PATIENT_ID", "label_time", "label_event"]).var().sort_values(ascending=False).head(10))

#final summary data matrix :)