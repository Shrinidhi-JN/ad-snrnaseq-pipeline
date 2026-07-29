import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("data/processed/GSE138852_qc_normalized.h5ad")

sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2)
print("Number of Leiden clusters found:", adata.obs["leiden"].nunique())

marker_sets = {
    "Oligodendrocyte": ["MOBP", "MBP", "PLP1"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2"],
    "OPC": ["PDGFRA", "CSPG4", "VCAN"],
    "Microglia": ["CSF1R", "C3", "CX3CR1"],
    "Neuron": ["SLC17A7", "GAD1", "GAD2"],
    "Endothelial": ["CLDN5", "FLT1"],
}

for cell_type, genes in marker_sets.items():
    present = [g for g in genes if g in adata.raw.var_names]
    print(f"{cell_type}: using markers {present} (missing: {set(genes) - set(present)})")
    sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{cell_type}", use_raw=True)

score_cols = [f"score_{ct}" for ct in marker_sets]
cluster_scores = adata.obs.groupby("leiden")[score_cols].mean()
cluster_scores.columns = list(marker_sets.keys())
best_type = cluster_scores.idxmax(axis=1)
print("\n--- Per-cluster marker scores ---")
print(cluster_scores.round(3))
print("\n--- Assigned cell type per cluster ---")
print(best_type)

adata.obs["cell_type"] = adata.obs["leiden"].map(best_type.to_dict())

print("\n--- Validation: our cell_type vs authors' author_cell_type ---")
print(pd.crosstab(adata.obs["cell_type"], adata.obs["author_cell_type"]))

adata.write_h5ad("data/processed/GSE138852_annotated.h5ad")
print("\nSaved -> data/processed/GSE138852_annotated.h5ad")
