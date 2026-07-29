# Cell-Type-Resolved Analysis of Alzheimer's Disease Brain Tissue (snRNA-seq)

A single-nucleus RNA-seq analysis pipeline applied to real published Alzheimer's
disease brain tissue, built around a specific methodological question: when the
biological unit of replication is a *donor* (or, in this dataset, a *pooled pair
of donors*), what does it take to draw statistically defensible conclusions —
rather than the inflated, easily-overturned "significance" you get from treating
every individual cell as an independent replicate?

## Data

[GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852) —
Chew, Grubman et al., *"A single-cell atlas of the human cortex reveals drivers
of transcriptional changes in Alzheimer's disease in specific cell
subpopulations"* (entorhinal cortex, 6 AD and 6 control donors, snRNA-seq,
10x Genomics).

**Important caveat about replication structure**, discovered while building
this pipeline (not stated up front in the processed file): the 12 donors were
sequenced as 6 pooled libraries of 2 individuals each, with no genotype-based
demultiplexing available in the processed data. This means the true unit of
independent replication available to us is the **library** (n=6: 3 AD, 3
control), not the individual (n=12). Every statistical step in this pipeline —
differential expression and classification alike — is built around that
constraint rather than around the more optimistic (and incorrect) assumption
of 12 independent replicates.

## Pipeline

1. **QC** (`scripts/run_qc_normalize.py`) — the deposited data is already a
   pre-filtered "high quality nuclei" set (13,214 nuclei; minimum genes-per-nucleus
   is 274, never near zero), so this step layers additional, more targeted QC on
   top: a 5% mitochondrial-read ceiling (appropriate for *nuclear*, not whole-cell,
   RNA — nuclei shouldn't contain much mitochondrial RNA at all) combined with the
   original authors' own doublet/unassigned-nucleus calls. Removes 1,398 nuclei
   (10.6%), leaving 11,816 for analysis.

2. **Clustering & annotation** (`scripts/run_clustering.py`) — Leiden clustering
   on the top 2,000 highly variable genes (40 PCs), followed by marker-gene-set
   scoring per cluster to assign cell types, with no dependency on any
   internet-hosted reference model. Validated directly against the original
   authors' manual cell-type calls: **98.3% overall agreement**. The one
   consistent weak point: endothelial cells (<1% of all nuclei) were mostly
   absorbed into the oligodendrocyte cluster rather than forming their own group
   — a known limitation of unsupervised clustering on rare populations, and
   excluded from the differential expression step below as a result (too few
   cells per library to build a reliable pseudobulk profile).

3. **Pseudobulk differential expression** (`scripts/run_pseudobulk_de.py`) —
   per cell type, raw counts are summed per library (not per cell) before
   running DESeq2, specifically to avoid the well-documented false-discovery
   inflation that comes from treating thousands of cells from the same 2 donors
   as independent replicates (Squair et al. 2021, *Nat. Commun.*). Selected
   findings from the real 3-vs-3 library comparison:
   - **Astrocytes** (202 significant genes, padj<0.05): GFAP up in AD
     (log2FC=2.07, padj≈0) — the classic marker of reactive astrogliosis,
     recovered here with no prior information fed into the model.
   - **Microglia** (27 significant genes): SPP1 significantly elevated — a
     well-known disease-associated-microglia marker from the mouse
     neurodegeneration literature, showing up unprompted in real human tissue.
   - **FKBP5** upregulated independently across three separate cell types
     (oligodendrocyte, OPC, microglia) — FKBP5/glucocorticoid signaling has
     documented links to tau pathology in AD.
   - **Neurons** showed a much weaker signal (4 significant genes) — likely
     an underpowered comparison given the combined excitatory+inhibitory
     category and relatively few neuronal nuclei captured (634 total). All
     three tested synaptic genes (SNAP25, SYT1, RAB3A) trended in the
     expected AD-down direction without reaching significance — consistent
     with low power rather than absence of effect.
   - Not everything matched prior expectation: APOE and CD74 in
     astrocytes/microglia came out statistically significant but in the
     *opposite* direction from the literature genes I checked against. Reported
     as-is rather than adjusted to match expectation — AD transcriptional
     effects on APOE are known to vary by brain region and cell subpopulation
     across studies.

4. **Donor-level disease classifier** (`scripts/run_classifier.py`) — features
   built at the library level (cell-type composition + pseudobulk expression)
   to avoid cell-level data leakage, evaluated with leave-one-library-out CV.
   Achieved a perfect ROC-AUC of 1.0 — **but** a 500-permutation label-shuffle
   test showed random labels achieve this same perfect separation about 8.6%
   of the time given the same pipeline and n=6, meaning this result is **not
   distinguishable from chance at conventional significance levels**. Reported
   deliberately as a negative/inconclusive result: with 6 samples and
   thousands of candidate features, a flexible classifier can perfectly
   separate randomly shuffled groups too often to trust a raw accuracy number
   without this kind of check.

## Reproducing this analysis

```bash
conda create -n adpipeline python=3.11 -y
conda activate adpipeline
pip install -r requirements.txt

mkdir -p data/raw data/processed
cd data/raw
curl -L -o GSE138852_counts.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138852&format=file&file=GSE138852%5Fcounts%2Ecsv%2Egz"
curl -L -o GSE138852_covariates.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138852&format=file&file=GSE138852%5Fcovariates%2Ecsv%2Egz"
gunzip -k GSE138852_counts.csv.gz GSE138852_covariates.csv.gz
cd ../..

python scripts/load_geo_data.py
python scripts/run_qc_normalize.py
python scripts/run_clustering.py
python scripts/run_pseudobulk_de.py
python scripts/run_classifier.py
```

## What this project demonstrates

Cell-type annotation without a black-box reference model; awareness of and
correction for pseudoreplication in single-cell differential expression;
donor-level (not cell-level) ML evaluation with leakage-free feature
selection; permutation-based significance testing appropriate for small
sample sizes; and — throughout — reporting findings that contradicted
expectations or came out statistically inconclusive, rather than only the
ones that looked good.

## License

MIT
