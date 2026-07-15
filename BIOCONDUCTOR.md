# Bioconductor docs for ClinicalVariantR

Start here:

1. **[inst/Bioconductor/SUBMISSION_CHECKLIST.md](inst/Bioconductor/SUBMISSION_CHECKLIST.md)** — what Bioconductor requires and what is already prepared  
2. **[inst/Bioconductor/PACKAGE_CONVERSION.md](inst/Bioconductor/PACKAGE_CONVERSION.md)** — convert this Shiny app tree into a package  
3. **[inst/Bioconductor/ISSUE_TEMPLATE_DRAFT.md](inst/Bioconductor/ISSUE_TEMPLATE_DRAFT.md)** — how to open the submission issue  

## Vignettes

| File | Topic |
|------|--------|
| `vignettes/ACMGamp.Rmd` | Introduction / overview |
| `vignettes/ACMGamp-run-and-test.Rmd` | App launch, sample VCFs, Group A/B/C UI tests, CLI & unit tests |

Build locally (from package root, with BiocStyle installed):

```r
# BiocManager::install("BiocStyle")
devtools::build_vignettes()
# or:
rmarkdown::render("vignettes/ACMGamp-run-and-test.Rmd")
```

Official portals:

- New submissions: https://github.com/Bioconductor/BiocContributions/issues  
- Guidelines: https://contributions.bioconductor.org/  
- Shiny rules: https://contributions.bioconductor.org/shiny.html  
