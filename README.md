# ACMGamp

**ACMGamp** is an R Shiny application for ACMG/AMP evidence-based germline variant interpretation in clinical genomics research.

Repository: [https://github.com/safarafique/ACMGamp](https://github.com/safarafique/ACMGamp)

## Description

ACMGamp provides a modular Shiny platform for variant pathogenicity assessment with two pipelines:

| Pipeline | Inputs | Criteria |
|----------|--------|----------|
| **Group A (Full)** | VCF + Clinical Logs CSV + Pedigree CSV | 28 ACMG/AMP criteria (automated + manual curation) |
| **Group B (Rapid)** | VCF only | Automated ACMG criteria |

Features include VCF requirement validation (green/red readiness), variant preview, **complete large-VCF streaming analysis** (no row cap), optional **bcftools** integration on Ubuntu/WSL, color-coded results, audit logging, and CSV export.

## Quick start

### Install R packages

```r
install.packages(c("shiny", "bslib", "DT", "data.table", "readr"))
```

### Run the app

```r
shiny::runApp("path/to/ACMGamp")
```

### Ubuntu / WSL (recommended for large VEP VCFs)

```bash
cd ACMGamp
bash scripts/ubuntu_setup.sh   # installs bcftools
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

## Project structure

```
ACMGamp/
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

- `Final_Clinical_Report.csv` — full classification report (Download button)
- `logs/analysis_log.csv` — timestamped audit trail

## Logic check

```r
setwd("path/to/ACMGamp")
source("global.R")
test_acmg_logic_engine()
```

## License

MIT (add `LICENSE` file if required for Bioconductor submission)
