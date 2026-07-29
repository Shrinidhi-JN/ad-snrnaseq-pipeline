import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

adata = sc.read_h5ad("data/processed/GSE138852_annotated.h5ad")

# Known AD-associated genes from the literature, used only as a sanity check
# on whatever DE results we get - not fed into the model in any way.
KNOWN_AD_GENES = {
    "Microglia": {"TREM2": "up", "APOE": "up", "C1QB": "up", "CD74": "up"},
    "Astrocyte": {"GFAP": "up", "APOE": "up", "CLU": "up", "VIM": "up"},
    "Neuron": {"SNAP25": "down", "SYT1": "down", "RAB3A": "down"},
}

cell_types = [ct for ct in adata.obs["cell_type"].unique() if ct != "Endothelial"]
all_results = {}

for cell_type in cell_types:
    sub = adata[adata.obs["cell_type"] == cell_type]

    pseudobulk, meta_rows = [], []
    for donor in sub.obs["donor_id"].unique():
        mask = (sub.obs["donor_id"] == donor).values
        n_cells = int(mask.sum())
        if n_cells < 10:
            continue
        pseudobulk.append(np.asarray(sub.layers["counts"][mask].sum(axis=0)).ravel())
        meta_rows.append({"donor_id": donor, "n_cells": n_cells,
                           "diagnosis": sub.obs["diagnosis"][mask].iloc[0]})

    counts_df = pd.DataFrame(pseudobulk, index=[m["donor_id"] for m in meta_rows], columns=sub.var_names).astype(int)
    meta_df = pd.DataFrame(meta_rows).set_index("donor_id")

    group_sizes = meta_df["diagnosis"].value_counts()
    print(f"\n=== {cell_type}: {group_sizes.to_dict()} libraries, {counts_df.shape[1]} genes ===")
    if group_sizes.min() < 2:
        print("  Skipping - fewer than 2 libraries in one group.")
        continue

    keep_genes = counts_df.sum(axis=0) >= 10
    counts_df = counts_df.loc[:, keep_genes]

    dds = DeseqDataSet(counts=counts_df, metadata=meta_df, design="~diagnosis", quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["diagnosis", "AD", "Control"], quiet=True)
    stats.summary()
    res = stats.results_df.sort_values("padj")
    all_results[cell_type] = res

    n_sig = (res["padj"] < 0.05).sum()
    print(f"  {n_sig} genes significant at padj < 0.05")
    print(res.head(10)[["log2FoldChange", "pvalue", "padj"]])

    if cell_type in KNOWN_AD_GENES:
        print(f"  --- Known AD gene check for {cell_type} ---")
        for gene, direction in KNOWN_AD_GENES[cell_type].items():
            if gene not in res.index:
                print(f"    {gene}: not detected")
                continue
            row = res.loc[gene]
            correct_dir = (row["log2FoldChange"] > 0) == (direction == "up")
            print(f"    {gene}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.3f}, "
                  f"expected={direction}, direction_correct={correct_dir}")

for cell_type, res in all_results.items():
    res.to_csv(f"results_de_{cell_type.replace(' ', '_')}.csv")

print("\nSaved per-cell-type DE result CSVs.")
