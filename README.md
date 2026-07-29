# Cell-Type-Resolved Analysis of Alzheimer's Disease Brain Tissue (snRNA-seq)

A single-nucleus RNA-seq pipeline built on real, published Alzheimer's disease brain
tissue. The question driving most of the design choices here: when your actual unit
of biological replication is a donor (or in this dataset, a pooled pair of donors),
what do you have to do differently to get statistically defensible results, instead
of the inflated significance you get from treating every individual cell as its own
independent replicate.

## Data

[GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852), from
Chew, Grubman et al., "A single-cell atlas of the human cortex reveals drivers of
transcriptional changes in Alzheimer's disease in specific cell subpopulations"
(entorhinal cortex, 6 AD and 6 control donors, snRNA-seq, 10x Genomics).

**A caveat about replication structure**, found partway through building this
pipeline and not obvious from the processed file alone: the 12 donors were
sequenced as 6 pooled libraries of 2 individuals each, with no genotype-based
demultiplexing available in the processed data. That means the real unit of
independent replication available here is the library (n=6: 3 AD, 3 control),
not the individual (n=12). Every statistical step below, differential expression
and classification alike, is built around that constraint rather than the more
convenient but wrong assumption of 12 independent replicates.

## Pipeline

1. **QC** (`scripts/run_qc_normalize.py`). The deposited data is already a
   pre-filtered "high quality nuclei" set (13,214 nuclei, minimum genes per
   nucleus is 274, never close to zero), so this step adds a second, more
   targeted layer of QC on top: a 5% mitochondrial-read ceiling (appropriate
   for nuclear RNA specifically, since nuclei shouldn't contain much
   mitochondrial RNA at all) combined with the original authors' own doublet
   and unassigned-nucleus calls. Removes 1,398 nuclei (10.6%), leaving 11,816
   for analysis.

2. **Clustering and annotation** (`scripts/run_clustering.py`). Leiden
   clustering on the top 2,000 highly variable genes (40 PCs), then marker
   gene set scoring per cluster to assign cell types, with no dependency on
   any internet-hosted reference model. Checked directly against the original
   authors' manual cell-type calls: 98.3% overall agreement. The one
   consistent weak point is endothelial cells (under 1% of all nuclei), which
   mostly got absorbed into the oligodendrocyte cluster instead of forming
   their own group. That's a known limitation of unsupervised clustering on
   rare populations, and it's why endothelial cells are excluded from the
   differential expression step below (too few cells per library to build a
   reliable pseudobulk profile).

3. **Pseudobulk differential expression** (`scripts/run_pseudobulk_de.py`).
   For each cell type, raw counts are summed per library, not per cell,
   before running DESeq2. This specifically avoids the false-discovery
   inflation that comes from treating thousands of cells from the same 2
   donors as independent replicates (Squair et al. 2021, *Nat. Commun.*). A
   few findings from the real 3-vs-3 library comparison:

   - Astrocytes (202 significant genes at padj<0.05): GFAP up in AD
     (log2FC=2.07, padj close to 0). That's the classic marker of reactive
     astrogliosis, and it came out of the model with no prior information
     fed in.
   - Microglia (27 significant genes): SPP1 significantly elevated, a
     well-known disease-associated-microglia marker from the mouse
     neurodegeneration literature, showing up on its own in real human
     tissue.
   - FKBP5 came up significant in three separate cell types independently
     (oligodendrocyte, OPC, microglia). FKBP5 and glucocorticoid signaling
     have documented links to tau pathology in AD.
   - Neurons showed a much weaker signal (4 significant genes), probably
     because excitatory and inhibitory neurons were combined into one
     category and relatively few neuronal nuclei were captured overall (634
     total). All three synaptic genes checked (SNAP25, SYT1, RAB3A) trended
     in the expected AD-down direction without reaching significance, which
     looks more like low power than an absent effect.
   - Not everything lined up with prior expectation. APOE and CD74 in
     astrocytes and microglia came out significant but in the opposite
     direction from what the literature suggested going in. That's reported
     as-is instead of adjusted to match expectation. AD's effect on APOE
     transcription is known to vary by brain region and cell subpopulation
     across different studies.

4. **Donor-level disease classifier** (`scripts/run_classifier.py`). Features
   built at the library level (cell-type composition plus pseudobulk
   expression) to avoid cell-level data leakage, evaluated with
   leave-one-library-out cross-validation. This got a perfect ROC-AUC of 1.0.
   But a 500-permutation label-shuffle test showed random labels hit that
   same perfect separation about 8.6% of the time with the same pipeline and
   n=6. So this result isn't distinguishable from chance at conventional
   significance levels. It's reported here as a negative or inconclusive
   result on purpose: with 6 samples and thousands of candidate features, a
   flexible enough classifier can separate randomly shuffled groups often
   enough that you can't trust a raw accuracy number without a check like
   this.

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

Cell-type annotation without relying on a black-box reference model.
Recognizing and correcting for pseudoreplication in single-cell differential
expression. Donor-level (not cell-level) ML evaluation with leakage-free
feature selection. Permutation-based significance testing suited to small
sample sizes. And throughout, reporting the findings that contradicted
expectations or came out statistically inconclusive, not just the ones that
looked good.

## License

MIT
