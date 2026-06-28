import pandas as pd

# Try reading moldata
df = pd.read_csv("moldata.csv")
print("Shape:", df.shape)
print("Columns:", df.columns.tolist())
print(df.head(3))