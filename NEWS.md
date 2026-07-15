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
* Vignettes: introduction (`ACMGamp.Rmd`) and run/test guide
  (`ACMGamp-run-and-test.Rmd`) covering app launch, sample data, UI checks, and
  CLI/unit tests.
