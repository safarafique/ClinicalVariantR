# ClinicalVariantR codebase map

Use this map as the first stop for navigation. For Shiny module ownership, see
`R/shiny/SHINY_MODULES.md`.

## Runtime entry points

| Path | Purpose |
|------|---------|
| `app.R` | Classic source-tree launcher for `shiny::runApp(".")`. |
| `global.R` | Loads runtime packages, constants, paths, and app modules in dependency order. |
| `ui.R` | Thin navbar composition; UI panels live under `R/shiny/ui/`. |
| `server.R` | Thin server composition; observer modules live under `R/shiny/`. |
| `R/ClinicalVariantR.R` | Exported Bioconductor-style constructors: `ClinicalVariantR()` and `ClinicalVariantRApp()`. |
| `global_cli.R` | Lightweight bootstrap for command-line scripts without Shiny UI packages. |

## Core R modules

| Area | Main files |
|------|------------|
| ACMG criteria and scoring | `R/acmg_engine.R`, `R/acmg_pipeline.R`, `R/acmg_vcf_criteria.R`, `R/classify_variant.R`, `R/rule_config.R` |
| Variant parsing and matching | `R/vcf_unified_parser.R`, `R/vcf_stream.R`, `R/vcf_validate.R`, `R/variant_key.R`, `R/variant_rescore.R` |
| Evidence sources | `R/reference_data.R`, `R/prediction_config.R`, `R/clinvar_pathogenic_criteria.R`, `R/autopvs1.R`, `R/ps4_case_control.R`, `R/hpo_omim_phenotype.R`, `R/trio_genotypes.R`, `R/clinical_context_criteria.R` |
| Reporting and export | `R/evidence_report.R`, `R/expert_review_export.R`, `R/reproducibility.R`, `R/audit.R`, `R/benchmark.R`, `R/intervar_compare.R` |
| UI helpers | `R/ui_helpers.R`, `R/gene_filter.R`, `R/parse_inputs.R`, `R/auth_storage.R` |

## Shiny modules

| Path | Purpose |
|------|---------|
| `R/shiny/context.R` | Shared session context and reactive state. |
| `R/shiny/shared_server.R` | Shared labels, navbar helpers, and result placeholders. |
| `R/shiny/upload_server.R` | VCF upload, validation, and deferred secure storage. |
| `R/shiny/analysis_server.R` | Shared analysis runner used by workflow tabs. |
| `R/shiny/group_a_server.R` | Clinical workflow with VCF, clinical CSV, pedigree, and manual evidence. |
| `R/shiny/group_b_server.R` | Rapid automated VCF-only workflow. |
| `R/shiny/group_c_server.R` | Gene-panel workflow. |
| `R/shiny/results_server.R` | Result tables, filters, downloads, and exports. |
| `R/shiny/explorer_server.R` | Evidence Explorer tab. |
| `R/shiny/audit_server.R` | Audit log tab. |
| `R/shiny/auth_server.R` | Optional authentication modal. |
| `R/shiny/ui/*.R` | Navbar panel builders and theme helpers. |

## Data and configuration

| Path | Purpose |
|------|---------|
| `config/` | ACMG thresholds, criteria, disease profiles, reference-install guide, optional auth users. |
| `data/reference/` | Placeholder and curated reference TSVs for evidence lookups. |
| `data/gene_panels/` | Gene panel and mechanism lists used by criteria. |
| `data/samples/` | Small example VCF/CSV inputs. |
| `data/validation/` | Expert-labelled validation set. |
| `www/custom.css` | App styling. |

## Scripts and checks

| Path | Purpose |
|------|---------|
| `scripts/verify_*.R` | Focused criteria and workflow verification scripts. |
| `scripts/generate_*benchmark*.R` | Benchmark fixture generation. |
| `scripts/run_*benchmark*.R` | Benchmark execution and comparison. |
| `scripts/compare_*.R`, `standalone/compare_clinicalvariantr_intervar.R` | InterVar/reference comparison tools. |
| `scripts/install_r_cli_deps.R`, `scripts/install_reference_data.R` | Local setup helpers. |
| `tests/testthat/` | Package tests for constructors, variant IDs, parsing, and reference matching. |

## Documentation

| Path | Purpose |
|------|---------|
| `README.md` | User-facing overview, launch instructions, workflows, and quick checks. |
| `vignettes/ClinicalVariantR.Rmd` | Bioconductor introduction vignette. |
| `vignettes/ClinicalVariantR-run-and-test.Rmd` | Run and test guide for app workflows and CLI checks. |
| `R/shiny/SHINY_MODULES.md` | Detailed Shiny module ownership map. |
| `inst/Bioconductor/` | Submission checklist, conversion plan, and issue draft. |
| `NEWS.md` | Release notes. |
| `man/` | Generated help pages. |

## Current conversion boundaries

The repository is still partly source-tree Shiny app and partly Bioconductor
package scaffold. Before submission, finish the layout work in
`inst/Bioconductor/PACKAGE_CONVERSION.md`: move package assets under `inst/`,
replace source chains with package functions, and verify clean `R CMD check`
plus `BiocCheck`.
