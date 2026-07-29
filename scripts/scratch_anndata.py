import anndata as ad
import numpy as np
import pandas as pd

X = np.array([
    [5, 0, 3, 1],
    [0, 0, 2, 0],
    [1, 4, 0, 0],
    [0, 5, 1, 0],
    [3, 0, 0, 2],
    [0, 1, 0, 3],
])

obs = pd.DataFrame({
    "donor_id": ["d1", "d1", "d2", "d2", "d3", "d3"],
    "diagnosis": ["AD", "AD", "Control", "Control", "AD", "Control"],
})

var = pd.DataFrame({
    "gene_role": ["marker", "marker", "filler", "filler"],
}, index=["GFAP", "AQP4", "GENE1", "GENE2"])

adata = ad.AnnData(X=X, obs=obs, var=var)

print(adata)
print("\nShape (cells x genes):", adata.shape)
print("\n--- obs (per-cell metadata) ---")
print(adata.obs)
print("\n--- var (per-gene metadata) ---")
print(adata.var)
print("\n--- Expression of GFAP across all cells ---")
print(adata[:, "GFAP"].X)
print("\n--- All cells from donor d1 ---")
print(adata[adata.obs["donor_id"] == "d1"].obs)
