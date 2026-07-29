import scanpy as sc

adata = sc.read_h5ad("data/processed/GSE138852_raw.h5ad")
print(adata)

mito_genes = [g for g in adata.var_names if g.startswith("MT-")]
print(f"\nMitochondrial genes found in this matrix: {len(mito_genes)}")
print(mito_genes)

adata.var["mt"] = adata.var_names.isin(mito_genes)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

print("\n--- QC metric distributions across all 13,214 nuclei ---")
print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe())

print("\n--- Author-assigned cell type breakdown ---")
print(adata.obs["author_cell_type"].value_counts())
