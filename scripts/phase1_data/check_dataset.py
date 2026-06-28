import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

df = pd.read_csv("master_dataset.csv")

print("=== Dataset Summary ===")
print(df[['HOMO','LUMO','PCE','Voc','Jsc']].describe().round(3))

# Check HOMO range per source
print("\n=== HOMO range by source ===")
print(df.groupby('source_db')['HOMO'].agg(['min','max','mean']).round(3))

# Check PCE range
print("\n=== PCE range by source ===")
print(df.groupby('source_db')['PCE'].agg(['min','max','mean']).round(3))

# Save a distribution plot
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
df['HOMO'].hist(ax=axes[0], bins=50, color='steelblue')
axes[0].set_title('HOMO distribution')
df['LUMO'].hist(ax=axes[1], bins=50, color='coral')
axes[1].set_title('LUMO distribution')
df['PCE'].hist(ax=axes[2], bins=50, color='green')
axes[2].set_title('PCE distribution')
plt.tight_layout()
plt.savefig("dataset_distributions.png", dpi=150)
print("\nSaved dataset_distributions.png")