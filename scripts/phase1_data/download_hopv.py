import deepchem as dc
import pandas as pd

print("Downloading HOPV dataset...")
tasks, datasets, transformers = dc.molnet.load_hopv(
    featurizer='Raw',
    splitter=None
)

train = datasets[0]
df = train.to_dataframe()
print("Columns:", df.columns.tolist())
print("Shape:", df.shape)

df.to_csv("hopv_raw.csv", index=False)
print("Saved to hopv_raw.csv")
print(df.head(3))