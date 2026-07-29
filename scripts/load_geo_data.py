import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

print("Loading counts matrix (genes x cells) - this may take a minute or two...")
counts = pd.read_csv("data/raw/GSE138852_counts.csv", index_col=0)
print("Raw shape (genes x cells):", counts.shape)
counts = counts.astype(np.int32)  # downcast after loading, not during parsing

counts = counts.T  # transpose -> cells x genes (the AnnData convention)
gene_names = counts.columns
cell_barcodes = counts.index

X = sp.csr_matrix(counts.values)
del counts

cov = pd.read_csv("data/raw/GSE138852_covariates.csv", index_col=0)
cov = cov.loc[cell_barcodes]

cov["library_id"] = cov.index.str.split("_", n=1).str[1]
diagnosis_map = {"AD": "AD", "ct": "Control"}

obs = pd.DataFrame({
    "donor_id": cov["library_id"].values,
    "diagnosis": pd.Categorical(
        cov["oupSample.batchCond"].map(diagnosis_map).values, categories=["Control", "AD"]
    ),
    "author_cell_type": cov["oupSample.cellType"].values,
}, index=cell_barcodes)

var = pd.DataFrame(index=gene_names)

adata = ad.AnnData(X=X, obs=obs, var=var)
print("\nFinal AnnData object:")
print(adata)
print("\nDiagnosis counts:", adata.obs["diagnosis"].value_counts().to_dict())
print("Donor/library counts:", adata.obs["donor_id"].value_counts().to_dict())

adata.write_h5ad("data/processed/GSE138852_raw.h5ad")
print("\nSaved -> data/processed/GSE138852_raw.h5ad")
