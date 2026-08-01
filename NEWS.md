# ClinicalVariantR 0.99.4 (2026-08-01)

* Parallel VCF scoring: up to 5 chunks concurrently (`parallel::mclapply` on
  Linux/WSL; PSOCK on Windows), with per-chunk error isolation and failed-chunk
  logging. I/O stays sequential; CSV writes stay on the main process.
* Auto chunk size from system RAM/CPU (and optional dataset size); configurable
  via UI or `CLINICALVARIANTR_CHUNK_SIZE` / `CLINICALVARIANTR_PARALLEL_CHUNKS`.
* Progress bar shows percent complete / variants processed; startup messaging
  simplified to "Running...".
* Fast analysis start: skip full VCF pre-count (use index or quick estimate).
* Add `R/parallel_pipeline.R` and `scripts/benchmark_parallel_chunks.R`.
* BiocCheck cleanup carried forward: avoid `cat` variable name, drop debug
  `message`/`sprintf` noise, CITATION/DESCRIPTION notes as needed.

# ClinicalVariantR 0.99.3 (2026-07-22)

* NAMESPACE: use selective `@importFrom` only (no wholesale `import(shiny)`,
  `import(DT)`, or `import(jsonlite)`) so loading no longer warns about
  replacing `renderDataTable`, `dataTableOutput`, or `validate`.
* `.Rbuildignore`: exclude duplicate `R/shiny_*` / `R/ui_*` helpers; the app
  runs from `inst/shinyapp/` and attaches only shiny + bslib at launch.
* Align app/engine version strings (`APP_VERSION`, `ACMG_PRO_ENGINE`) with
  package `DESCRIPTION` version **0.99.3**.
* Vignettes: Bioconductor vs GitHub install/reinstall paths and app step summary.
* Check NOTES: ignore `LICENSE.md`; ship only `ClinicalVariantR*.Rmd` vignettes.
* Remove unused `crosstalk` dependency; keep dependency bootstrap only under
  `scripts/` (not shipped in the Bioconductor tarball).
* Fix package Rd wording that broke PDF manual generation; keep `DESCRIPTION`
  free of blank lines between fields.
* BiocCheck: drop `inst/CITATION` until a preprint/publication DOI exists;
  expand `Description` to three sentences; avoid `paste` in a condition message.

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
