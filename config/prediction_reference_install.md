# Prediction reference data

For production-grade predictions, replace placeholder files under `data/reference/`:

| File | Purpose | Suggested source |
|------|---------|------------------|
| `gnomad_v41_placeholder.tsv` → `gnomad_v41.tsv` | Population AF (BA1, BS1, BS2, PM2) | gnomAD v4 constraint / AF export |
| `clinvar_placeholder.tsv` → `clinvar.tsv` | ClinVar classifications (PP5, BP6) | ClinVar VCF/TSV release |
| `revel_placeholder.tsv` → `revel.tsv` | REVEL scores (PP3, BP4) | dbNSFP or REVEL transcript file |

Update `REFERENCE_PATHS` in `global.R` after installing full files.

The home page **Reference readiness** banner shows whether files look like placeholders (<100 rows).
