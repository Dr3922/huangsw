# Data sources and public-content boundary

## Public source datasets

- GSE32894, GSE13507, and GSE31684 for bulk expression validation.
- GSE169379 for epithelial single-nucleus analysis.
- GSE171351 for linked Visium spatial analysis.
- cBioPortal/TCGA-BLCA study releases for release-level FGFR3 mutation and CNA summaries.

Accession identifiers and the exact release-level metadata used in the revision are
recorded in the derived ledgers and the manuscript's Supplementary Data 1.

## Included in this release candidate

- nonidentifying sample and donor ledgers;
- derived expression scores and model summaries;
- spatial fold manifests, marker lists, and spot-level score tables;
- figure source data and analysis manifests;
- URLs and expected hashes for large public inputs.

## Excluded

- raw h5ad objects and raw GEO matrices;
- controlled-access or private clinical data;
- manuscript and response documents;
- internal audit, strategy, and writing files;
- large high-resolution histology arrays unless explicitly approved for the final release.

The public repository should cite the original GEO and cBioPortal records and should not
imply ownership or redistribution rights over those source datasets.

