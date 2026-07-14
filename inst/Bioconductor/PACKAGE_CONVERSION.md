# ACMGamp → Bioconductor package conversion plan

Use this checklist after reading `SUBMISSION_CHECKLIST.md`. Goal: convert the
current Shiny **application** tree into a standards-compliant Bioconductor
**software package** without changing scientific logic.

## Target tree

```
ACMGamp/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE / LICENSE.md
├── NEWS.md
├── README.md
├── app.R                    # optional thin launcher only
├── R/
│   ├── ACMGamp-package.R
│   ├── ACMGamp.R            # exported ACMGamp() / ACMGampApp()
│   ├── utils_*.R            # parsers, engines, helpers
│   ├── interface_*.R        # UI builders (from R/shiny/ui/)
│   ├── observers_*.R        # observeEvent modules
│   └── outputs_*.R          # render* / download* modules
├── man/
├── vignettes/
├── tests/testthat/
└── inst/
    ├── CITATION
    ├── www/                 # from www/
    ├── config/              # from config/
    ├── extdata/             # small sample VCF + placeholders from data/
    └── scripts/             # optional small helper scripts only
```

## Step-by-step

### A. Path resolution

1. Replace hard-coded `file.path("config", ...)` / `file.path("data", ...)` with
   `system.file("config", ..., package = "ACMGamp")` and
   `system.file("extdata", ..., package = "ACMGamp")`.
2. Write logs to `tempdir()` or user-configurable paths, never into the package
   library tree.
3. Load `www` CSS via `system.file("www", "custom.css", package = "ACMGamp")`.

### B. Stop sourcing files

1. Remove `lapply(module_files, source, ...)` from `global.R`.
2. All logic under `R/` must be package functions (roxygen + NAMESPACE).
3. Keep only `ACMGamp()` constructing `shiny::shinyApp(ui, server)`.

### C. Relocate assets

| Move from | Move to |
|-----------|---------|
| `config/*` | `inst/config/` |
| `data/samples/*` (small) | `inst/extdata/samples/` |
| `data/reference/*` placeholders | `inst/extdata/reference/` |
| `data/gene_panels/*` | `inst/extdata/gene_panels/` |
| `data/validation/*` (if small) | `inst/extdata/validation/` |
| `www/*` | `inst/www/` |
| `scripts/*` | GitHub-only branch **or** `inst/scripts/` if ≤5 MB and needed |

Do **not** commit large clinical VCFs, InterVar databases, or `results/`.

### D. Shiny file naming (recommended)

Map current modules gradually:

| Current | Suggested |
|---------|-----------|
| `R/shiny/ui/*_ui.R` | `R/interface_*.R` |
| `R/shiny/*_server.R` | `R/observers_*.R` / `R/outputs_*.R` |
| classification / VCF engines | `R/utils_*.R` |

Reviewers allow modules; naming scheme improves BiocCheck friendliness.

### E. `app.R` (optional)

After conversion, root `app.R` may contain only:

```r
shiny::runApp(ACMGamp::ACMGamp())
```

Do not put business logic in `app.R`.

### F. Auth defaults

1. Remove or sanitize `config/auth_users.csv` default passwords.  
2. Prefer `AUTH_ENABLED = FALSE` for Bioconductor demos, or environment-based credentials.  
3. Document encryption keys via env vars only.

### G. Tests to add first (high value, non-UI)

- Combining rules (`combine_acmg_evidence` / rules audit cases).  
- PM2/BS2 mutual exclusion.  
- VEP CSQ parsing on a 5–20 variant `inst/extdata` fixture.  
- Classification tier labels on the curated 20-variant benchmark (if size OK).

### H. Version policy

- First Bioconductor submission: `0.99.0`.  
- Each rebuild during review: bump to `0.99.1`, `0.99.2`, …  
- After acceptance on devel: Bioconductor will set release versioning.

### I. Repository hygiene

Default branch must be package-only. Move CI configs, manuscript drafts,
benchmark raw dumps, and lab notes to another branch (e.g. `devtools`).

### J. Acceptance criterion before opening the issue

```text
[ ] DESCRIPTION Version 0.99.z + biocViews
[ ] system.file paths work after install
[ ] ACMGamp() returns shiny.appobj; vignette builds
[ ] R CMD build / check clean
[ ] BiocCheck + BiocCheckGitClone clean (or justified notes)
[ ] No file > 5 MB; tarball < 10 MB
[ ] Maintainer email + ORCID finalized
```

When all boxes are checked, open the BiocContributions issue using
`ISSUE_TEMPLATE_DRAFT.md`.
