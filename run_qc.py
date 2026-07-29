import numpy as np
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("data/processed/GSE138852_raw.h5ad")

mito_genes = [g for g in adata.var_names if g.startswith("MT-")]
adata.var["mt"] = adata.var_names.isin(mito_genes)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

def mad_outlier(values, n_mads, one_sided_upper=False):
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    upper = med + n_mads * mad
    lower = -np.inf if one_sided_upper else med - n_mads * mad
    return (values < lower) | (values > upper)

# Counts/genes: outlier in EITHER direction (too low = empty/degraded nucleus,
# too high = possible doublet), evaluated on log scale since count data is skewed.
counts_outlier = mad_outlier(np.log1p(adata.obs["total_counts"]), n_mads=5)
genes_outlier = mad_outlier(np.log1p(adata.obs["n_genes_by_counts"]), n_mads=5)
# Mito %: only HIGH is bad (low mito is perfectly fine/expected for nuclei).
mito_outlier = mad_outlier(adata.obs["pct_counts_mt"], n_mads=3, one_sided_upper=True)

adata.obs["qc_outlier"] = counts_outlier | genes_outlier | mito_outlier

print("Flagged by total_counts:", counts_outlier.sum())
print("Flagged by n_genes:", genes_outlier.sum())
print("Flagged by pct_mt:", mito_outlier.sum())
print("Combined QC-metric outliers:", adata.obs["qc_outlier"].sum(), "/", adata.n_obs)

# Validation check: do our statistically-flagged outliers overlap with the
# authors' own doublet/unID calls? They SHOULD overlap somewhat if our
# metric-based approach is catching real problem cells rather than noise.
author_flagged = adata.obs["author_cell_type"].isin(["doublet", "unID"])
print("\n--- Cross-tab: our QC outlier flag vs authors' doublet/unID label ---")
print(pd.crosstab(adata.obs["qc_outlier"], author_flagged, rownames=["our_outlier"], colnames=["author_doublet_or_unID"]))
