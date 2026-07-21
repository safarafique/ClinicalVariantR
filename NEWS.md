# ClinicalVariantR 0.99.2 (2026-07-20)

* Fix analysis crash on variants with pathogenic VEP `CLIN_SIG` in CSQ:
  `rbind` failed when merging ClinVar protein DB (extra `source` column) with
  the CSQ catalog (`numbers of columns of arguments do not match`).
* R CMD check: `.Rbuildignore` now uses `^data/` (and similar) so development
  non-R assets under `data/` are excluded; package samples remain in
  `inst/extdata/`. Roxygen `@importFrom` tags keep NAMESPACE Imports wired.

# ClinicalVariantR 0.99.1 (2026-07-17)

* R CMD check / R-universe WARNING cleanup:
  - Ship non-R assets under `inst/extdata/` (and ignore package-root `data/`).
  - ASCII-only package `R/` sources (em/en dashes and similar replaced).
  - Declare `digest` / `openssl` Imports used by auth storage.
  - Wire `NAMESPACE` imports for `bslib`, `VariantAnnotation`, and `methods`.
  - Declare optional `shinyjs` in Suggests (Group B run-button enable/disable).
  - Add `.Rbuildignore` for development-only root Shiny/config/scripts trees.

# ClinicalVariantR 0.99.0 (2026-07-08)

* Initial Bioconductor submission scaffold (`x.99.y` versioning).
* Interactive Shiny platform for ACMG/AMP germline variant classification from
  VEP / SnpEff / ANNOVAR-annotated VCFs.
* Three workflows: full clinical (Group A), automated rapid (Group B), and
  gene-panel (Group C).
* Streaming whole-VCF analysis with optional bcftools acceleration.
* Structured criterion-level evidence export and reproducibility metadata.
* Package entry points: `ClinicalVariantR()` / `ClinicalVariantRApp()` return a Shiny app object
  (user launches with `shiny::runApp()`).
* Vignettes: introduction and run/test guide covering app launch, sample data,
  UI checks, and CLI/unit tests.
