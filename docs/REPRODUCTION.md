# Reproduction guide

## Derived-data smoke test

From the repository root:

```powershell
python -m pip install -r requirements.txt
python scripts/smoke_reproduce.py --handoff data/derived
```

The test checks the principal reported sample counts, bulk effect estimates, DerSimonian-
Laird synthesis, donor-cluster bootstrap intervals, BH correction, TCGA bridge denominator,
and four strict spatial correlations.

To regenerate figure outputs:

```powershell
$env:FGFR3_OUT = (Resolve-Path "data/derived").Path
python scripts/smoke_reproduce.py --handoff data/derived --render
```

Figure 4 additionally needs either the optional spatial asset arrays under
`data/derived/figure_data/spatial_assets/` or `FGFR3_SPATIAL_H5AD`.

## Full analysis inputs

The full scripts read public raw/processed inputs through environment variables. The
following are intentionally not bundled:

- GSE169379_MIBC_snSeq.h5ad
- GSE171351_combined_visium.h5ad
- processed GEO matrices and platform annotations
- cBioPortal API caches
- historical CNV-window score files

Use `scripts/download_public_inputs.py` for the two large GEO objects, then set
`FGFR3_SCRNA_H5AD`, `FGFR3_SPATIAL_H5AD`, `FGFR3_PROJECT_ROOT`, and `FGFR3_OUT`.

## Spatial Round 2

```powershell
python scripts/run_targeted_v2.py `
  --handoff data/derived `
  --spatial-h5ad D:/data/GSE171351_combined_visium.h5ad `
  --out data/derived/spatial_validation_v2
```

This run is strict leave-one-patient-out and records the score-ablation direction
reversal. It is a descriptive sensitivity layer, not an independent patient-level
validation study.

