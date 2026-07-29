import scanpy as sc

adata = sc.read_h5ad("data/processed/GSE138852_raw.h5ad")

mito_genes = [g for g in adata.var_names if g.startswith("MT-")]
adata.var["mt"] = adata.var_names.isin(mito_genes)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

exclude = (adata.obs["pct_counts_mt"] > 5.0) | adata.obs["author_cell_type"].isin(["doublet", "unID"])
print(f"Excluding {exclude.sum()} / {adata.n_obs} cells")
adata = adata[~exclude].copy()

sc.pp.filter_genes(adata, min_cells=3)  # drop genes expressed in almost no cells
print(f"Shape after QC filtering: {adata.shape}")

adata.layers["counts"] = adata.X.copy()  # keep raw counts for later pseudobulk DE

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # keep full log-normalized data for marker scoring later

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
print(f"Highly variable genes selected: {adata.var['highly_variable'].sum()}")

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=40)
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

adata.write_h5ad("data/processed/GSE138852_qc_normalized.h5ad")
print("Saved -> data/processed/GSE138852_qc_normalized.h5ad")
