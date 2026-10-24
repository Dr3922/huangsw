# FGFR3 urothelial cancer — analysis code and derived data

Analysis scripts and nonidentifying derived data underlying the study **“FGFR3-high urothelial programs are associated with attenuated epithelial antigen-presentation transcription and heterogeneous spatial immune context.”**

## Contents

- `scripts/`: cohort selection, module scoring, regression, random-effects synthesis, donor-cluster bootstrap, strict spatial leave-one-patient-out (LOPO) projection, and figure-generation scripts
- `data/derived/bulk/`: cohort-internal module scores and bulk model outputs
- `data/derived/cohort_ledgers/`: source-study identifiers and eligibility information
- `data/derived/snrna/`: epithelial pseudobulk, state-specific results, and bootstrap draws (gzip-compressed CSV)
- `data/derived/spatial/`: fold/section mapping, projection scores, and spatial sensitivity results
- `data/derived/figure_data/`: the numeric data underlying published figure panels
- `data/derived/supplementary_tables/`: machine-readable supplementary model and ledger tables
- `docs/`: input-data source descriptions, analysis scope, and reproducibility instructions

The original GEO/Visium `h5ad` data objects, cBioPortal cache, and image arrays are not redistributed here. Retrieve original public datasets from the named accessions in the manuscript.

## Minimal verification

Using Python 3.13 or a compatible recent Python environment:

```bash
python -m pip install -r requirements.txt
python scripts/smoke_reproduce.py --handoff data/derived
```

This verifies sample counts, selected effect estimates, meta-analysis calculations, donor-cluster bootstrap and LOPO spatial score summaries from the distributed *derived data*. An independent end-to-end raw-matrix analysis requires source downloads and optional packages in `requirements-full.txt` (see `docs/REPRODUCTION.md`).

## Scientific scope

Four linked Visium sections offer descriptive spatial transcriptomic information. The composite projection score shows substantial component-construction sensitivity, documented in the accompanying sensitivity data. The epithelial analyses describe HLA-I/APM-related *RNA* signals. Neither assay measures protein-level antigen presentation or dendritic-cell abundance.

## Citation and licensing

Citation metadata are in `CITATION.cff`. This repository does not include a software license because no license has been selected by the authors. No reuse license is asserted here.
