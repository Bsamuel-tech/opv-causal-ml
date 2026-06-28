import pandas as pd

df = pd.read_csv("hopv_raw.csv")

# y1-y8 are: HOMO, LUMO, electrochemical_gap, optical_gap, PCE, Voc, Jsc, fill_factor
df = df.rename(columns={
    'ids': 'SMILES',
    'y1': 'HOMO',
    'y2': 'LUMO',
    'y3': 'electrochemical_gap',
    'y4': 'optical_gap',
    'y5': 'PCE',
    'y6': 'Voc',
    'y7': 'Jsc',
    'y8': 'fill_factor'
})

# Keep only science columns
df = df[['SMILES', 'HOMO', 'LUMO', 'electrochemical_gap',
         'optical_gap', 'PCE', 'Voc', 'Jsc', 'fill_factor']]

# Add source label
df['source_db'] = 'HOPV15'

# Drop rows with all NaN targets
df = df.dropna(subset=['HOMO', 'LUMO', 'PCE'], how='all')

print("Shape after cleaning:", df.shape)
print("\nColumn names:", df.columns.tolist())
print("\nSample:")
print(df[['SMILES', 'HOMO', 'LUMO', 'PCE']].head(5))
print("\nMissing values:")
print(df.isnull().sum())

df.to_csv("hopv_clean.csv", index=False)
print("\nSaved to hopv_clean.csv")