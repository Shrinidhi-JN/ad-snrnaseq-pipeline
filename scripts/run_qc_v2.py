import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("data/processed/GSE138852_raw.h5ad")
mito_genes = [g for g in adata.var_names if g.startswith("MT-")]
adata.var["mt"] = adata.var_names.isin(mito_genes)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

fixed_mito_outlier = adata.obs["pct_counts_mt"] > 5.0
author_flagged = adata.obs["author_cell_type"].isin(["doublet", "unID"])

print("Cells over 5% mito (fixed threshold):", fixed_mito_outlier.sum())
print("\n--- Cross-tab: fixed 5% mito filter vs author doublet/unID ---")
print(pd.crosstab(fixed_mito_outlier, author_flagged, rownames=["over_5pct_mt"], colnames=["author_doublet_or_unID"]))

final_exclude = fixed_mito_outlier | author_flagged
print("\nTotal cells excluded (mito>5% OR author doublet/unID):", final_exclude.sum())
print("Cells remaining for analysis:", (~final_exclude).sum())
