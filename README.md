# ClinicalVariantR — Variant Prediction Platform

**ClinicalVariantR** is an R Shiny application for ACMG/AMP evidence-based germline variant **prediction** with per-variant curator workflow.

Repository: [https://github.com/safarafique/ClinicalVariantR](https://github.com/safarafique/ClinicalVariantR)

## Bioconductor submission docs

Package metadata and submission guides live under `inst/Bioconductor/`:

| Document | Purpose |
|----------|---------|
| [`inst/Bioconductor/SUBMISSION_CHECKLIST.md`](inst/Bioconductor/SUBMISSION_CHECKLIST.md) | Required files + readiness status |
| [`inst/Bioconductor/PACKAGE_CONVERSION.md`](inst/Bioconductor/PACKAGE_CONVERSION.md) | Convert Shiny-app layout → Bioconductor package |
| [`inst/Bioconductor/ISSUE_TEMPLATE_DRAFT.md`](inst/Bioconductor/ISSUE_TEMPLATE_DRAFT.md) | Draft for BiocContributions issue |

Core package files: `DESCRIPTION` (version `0.99.0`), `NAMESPACE`, `LICENSE`, `NEWS.md`, `man/`, `vignettes/ClinicalVariantR.Rmd`, `tests/`, `inst/CITATION`.

Install and launch (Bioconductor Shiny style — return app, user runs it).
On Windows, **restart R first** if the package is already loaded:

```r
# install.packages("remotes")
# Preferred while developing from this clone:
remotes::install_local("E:/ACGM/ClinicalVariantR", force = TRUE, upgrade = "never")
# Or from GitHub after the latest launcher is pushed:
# remotes::install_github("safarafique/ClinicalVariantR")
library(ClinicalVariantR)
app <- ClinicalVariantR()
if (interactive()) shiny::runApp(app, launch.browser = TRUE)

# Fast path without reinstall:
# shiny::runApp("E:/ACGM/ClinicalVariantR", launch.browser = TRUE)
```

**Important:** documentation scaffold is ready; full package conversion (paths to `inst/`, remove `source()` chains, clean `R CMD check` / `BiocCheck`) must finish before opening a BiocContributions issue.

## Pipelines

| Pipeline | Inputs | Scope |
|----------|--------|--------|
| **Group A (Clinical prediction)** | VCF + Clinical CSV + Pedigree CSV | 28 criteria + per-variant curation |
| **Group B (Automated prediction)** | VCF only | 18 automated criteria |
| **Group C (Gene panel)** | VCF + gene list | Group B engine, panel-filtered |

**Prediction mode (v2.5+):** HPO/OMIM phenotype matching (PP4), trio genotype parsing (PS2/PM3/BP2), ClinGen AutoPVS1, login + encrypted uploads, formal validation report.

**Prediction mode (v2.4+):** stricter PVS1 (LoF panel required), evidence strength labels, prediction limitations per variant, reference readiness checks.

See `config/prediction_reference_install.md` for full gnomAD/ClinVar/REVEL setup.

### Authentication and PHI storage (v2.5)

- Authentication is disabled for local development by default (`AUTH_ENABLED <- FALSE` in `R/auth_storage.R`). For secured demos, provide users through `config/auth_users.csv` or the `CLINICALVARIANTR_USER` / `CLINICALVARIANTR_PASSWORD` environment variables.
- Override with environment variables `CLINICALVARIANTR_USER` and `CLINICALVARIANTR_PASSWORD` when no user file is present.
- Uploads are encrypted at rest (AES-256-GCM via OpenSSL) under `logs/secure_uploads/`; set `CLINICALVARIANTR_ENCRYPTION_KEY` (32 chars) for production.
- Access events are logged to `logs/access_audit.csv`.

### Validation report

```bash
cd ClinicalVariantR
Rscript scripts/generate_validation_report.R
```

Expert-classified variants live in `data/validation/clinicalvariantr_validation_set.tsv` (expand toward 500+ for production sign-off).

Features include VCF requirement validation (green/red readiness), variant preview, **complete large-VCF streaming analysis** (no row cap), optional **bcftools** integration on Ubuntu/WSL, color-coded results, audit logging, and CSV export.

## Quick start

### Install and open the app

```r
install.packages(c("shiny", "bslib", "DT", "data.table", "readr", "remotes"))
# Session -> Restart R if ClinicalVariantR is already loaded
remotes::install_local("E:/ACGM/ClinicalVariantR", force = TRUE, upgrade = "never")
library(ClinicalVariantR)
app <- ClinicalVariantR()
if (interactive()) shiny::runApp(app, launch.browser = TRUE)
```

### Run from a local clone (no reinstall)

```r
shiny::runApp("E:/ACGM/ClinicalVariantR", launch.browser = TRUE)
```

### Ubuntu / WSL (recommended for large VEP VCFs)

```bash
cd /mnt/e/ACGM/ClinicalVariantR
bash scripts/ubuntu_setup.sh          # bcftools + R CLI packages
# Or install R CLI deps only:
Rscript scripts/install_r_cli_deps.R  # data.table, readr, jsonlite → ~/R/.../library
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

Verify scripts (after `install_r_cli_deps.R` completes):

```bash
Rscript scripts/verify_group_a_28.R
Rscript scripts/verify_group_b_c.R
Rscript scripts/generate_validation_report.R
```

### Expert review (after prediction)

Each pipeline includes an **Expert review checklist** accordion with worklist counts and filtered CSV exports:

- **Export LP+ only** — Pathogenic and Likely Pathogenic variants (priority review).
- **Export expert worklist** — LP+ plus VUS with ≥2 pathogenic evidence; excludes PM2-only VUS.

Use the full prediction report for audit; use the worklist export to start manual sign-out review.

## Project structure

For a fuller file-by-file orientation, see `CODEBASE_MAP.md`.

```
ClinicalVariantR/
├── app.R
├── global.R
├── ui.R
├── server.R
├── R/                  # ACMG logic, VCF parsing, validation, streaming
├── data/reference/     # gnomAD, ClinVar, REVEL placeholders
├── data/samples/       # Example VCF and CSV files
├── scripts/            # ubuntu_setup.sh, validate_logic.R
├── logs/               # Runtime audit logs (gitignored)
└── www/                # Custom CSS
```

## Complete VCF analysis options

| Option | Default | Purpose |
|--------|---------|---------|
| Analyze entire VCF | ON | Process every row (no limit) |
| FILTER=PASS only | OFF | Optional QC pre-filter |
| Minimum QUAL | 0 | Drop low-quality rows |
| Use bcftools | auto | Ubuntu/WSL speed boost |
| Chunk size | 10,000 | Variants per memory batch |

## Outputs

**Best input for expert review:** download the ClinicalVariantR **CSV export** after analysis completes.

| Export | When to use |
|--------|-------------|
| `ClinicalVariantR_Expert_Worklist_*.csv` | **Start here** — LP+ and high-priority VUS for sign-out review |
| `ClinicalVariantR_LP_Plus_*.csv` | Pathogenic / Likely Pathogenic only |
| `ClinicalVariantR_Prediction_Report_*.csv` | Full audit trail — all variants and evidence |

Key columns for expert review: `variant_id`, `gene`, `classification`, `criteria_met`, `criteria_rationale`, `prediction_limitations`, `evidence_json`, `confidence_score`, `evidence_strength`.

Also written automatically under `logs/` during streaming runs. Audit trail: `logs/analysis_log.csv`.

## Logic check

```r
setwd("path/to/ClinicalVariantR")
source("global.R")
test_acmg_logic_engine()
```

## License

MIT (add `LICENSE` file if required for Bioconductor submission)
