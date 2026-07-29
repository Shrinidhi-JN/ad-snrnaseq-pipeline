import pandas as pd

cov = pd.read_csv("data/raw/GSE138852_covariates.csv", index_col=0)
print("Shape:", cov.shape)
print("\nColumns:", list(cov.columns))

print("\n--- condition counts ---")
print(cov["oupSample.batchCond"].value_counts())

print("\n--- cell type counts (author-assigned) ---")
print(cov["oupSample.cellType"].value_counts())

# Barcode format looks like <16bp_barcode>_<library_id>. Extract the library/pool ID.
cov["library_id"] = cov.index.str.split("_", n=1).str[1]
print("\n--- unique library/pool IDs and their condition ---")
print(cov.groupby("library_id")["oupSample.batchCond"].first())
print("\nNuclei per library:")
print(cov["library_id"].value_counts())
