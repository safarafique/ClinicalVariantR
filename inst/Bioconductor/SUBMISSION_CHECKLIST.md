# Bioconductor submission (ClinicalVariantR)

This document lists **required package files**, the **current status** of this
repository, and **remaining conversion work** before opening a submission issue
at [Bioconductor/BiocContributions](https://github.com/Bioconductor/BiocContributions/issues).

Official guides:

- [Bioconductor package submissions](https://contributions.bioconductor.org/bioconductor-package-submissions.html)
- [General package development](https://contributions.bioconductor.org/general.html)
- [Shiny apps](https://contributions.bioconductor.org/shiny.html)
- [Submitter guide (R-universe process)](https://github.com/Bioconductor/BiocContributions/blob/devel/docs/submitters.md)
- [Package review checklist](https://contributions.bioconductor.org/review-checklist.html)

---

## 1. Required documentation / metadata (created in this scaffold)

| Item | Path | Status |
|------|------|--------|
| Package metadata | `DESCRIPTION` | ✅ `Version: 0.99.0`, `biocViews`, MIT |
| Namespace | `NAMESPACE` | ✅ minimal exports (`ClinicalVariantR`, `ClinicalVariantRApp`) |
| License | `LICENSE`, `LICENSE.md` | ✅ MIT |
| News | `NEWS.md` | ✅ |
| Citation | `inst/CITATION` | ✅ maintainer ORCID/email recorded |
| Package man page | `man/ClinicalVariantR-package.Rd` | ✅ |
| Exported function man pages | `man/ClinicalVariantR.Rd` | ✅ |
| Vignette | `vignettes/ClinicalVariantR.Rmd` | ✅ draft (needs BiocStyle) |
| Unit test scaffold | `tests/testthat.R`, `tests/testthat/` | ✅ minimal |
| Bioconductor checklist | `inst/Bioconductor/SUBMISSION_CHECKLIST.md` | ✅ this file’s companion |
| Submission issue draft | `inst/Bioconductor/ISSUE_TEMPLATE_DRAFT.md` | ✅ |
| Conversion plan | `inst/Bioconductor/PACKAGE_CONVERSION.md` | ✅ |

**Maintainer (current):** Safa Rafique `<safa.sandhu@gmail.com>`.
**Maintainer ORCID:** <https://orcid.org/0000-0003-2646-8106>.
**Authors:** Safa Rafique, Naeem Mahmood, Muhammad Farooq Sabar.
Add remaining author ORCID IDs and affiliations in `DESCRIPTION` / vignette if available before submit.

---

## 2. Precheck validation (automatic when you open the GitHub issue)

Your default branch must satisfy:

1. Public GitHub URL  
2. `DESCRIPTION` + `vignettes/` present  
3. Fields: `Package`, `Version`, `biocViews`  
4. **Package name == GitHub repository name** (case-sensitive) → repo is `ClinicalVariantR` ✅  
5. Version **`x.99.y`** → `0.99.0` ✅  
6. **No `Remotes:`** / **No `Additional_repositories`** ✅  
7. No file **> 5 MB**  
8. No Git LFS  

Then comment exactly: `/accept-policies`

---

## 3. Critical conversion still required (NOT optional)

Bioconductor reviewers treat ClinicalVariantR as a **software package**, not a loose
Shiny project. The following must be finished **before** clean `R CMD check`
and `BiocCheck`:

### 3.1 Package layout

| Current (app layout) | Required (package layout) |
|----------------------|---------------------------|
| `config/`, `data/`, `www/` at package root | Move into `inst/` (e.g. `inst/config`, `inst/extdata`, `inst/www`) |
| `global.R` + `source("R/...")` | Package functions only; remove `source()` chains |
| `ui.R` / `server.R` at root | Functions under `R/` returning UI/server; optional thin `app.R` |
| `scripts/`, `standalone/`, large bench VCFs | Keep **out of default branch** or under `inst/scripts` if small |
| `logs/`, `results/` | Never ship; already gitignored |

See `PACKAGE_CONVERSION.md` for a step-by-step move plan.

### 3.2 Shiny policy ([Chapter 18](https://contributions.bioconductor.org/shiny.html))

- ✅ Entry points return `shinyApp` (`ClinicalVariantR()` / `ClinicalVariantRApp()`).  
- ⚠️ Do **not** call `shiny::runApp()` inside package functions (except optional interactive `launch = TRUE` — prefer documenting user-side `runApp`).  
- ⚠️ Prefer Bioconductor naming: `interface_*.R`, `observers_*.R`, `outputs_*.R`, `utils_*.R` (or migrate existing `R/shiny/` modules toward that convention).  
- Man-page examples must wrap launch in `if (interactive())`.

### 3.3 Checks that must be clean

```r
# From a parent directory of the package source:
R CMD build ClinicalVariantR
R CMD check ClinicalVariantR_0.99.0.tar.gz
# Bioconductor-specific:
BiocCheck::BiocCheckGitClone("ClinicalVariantR")
BiocCheck::BiocCheck("ClinicalVariantR_0.99.0.tar.gz", `new-package` = TRUE)
```

Expectations: **no ERROR / WARNING**; justify any NOTE.

Also:

- Source tarball **< 10 MB**  
- `R CMD check --no-build-vignettes` **< 10 minutes**  
- Individual files **≤ 5 MB**  
- Memory for vignette/examples/tests **< ~8 GB**

### 3.4 Tests and documentation

- Expand `tests/testthat/` for non-reactive engine functions (combining rules,
  parsers, classification helpers).  
- Use `# nocov` around pure Shiny observers if needed.  
- Optionally add `shinytest2` for UI smoke tests.  
- Document **exported** functions with runnable examples; mark internals
  `@keywords internal`.

### 3.5 Auth / default credentials

`config/auth_users.csv` contains default passwords. For Bioconductor:

- Do **not** ship real credentials.  
- Prefer example users only in docs, or generate on first run with env vars.  
- Document `AUTH_ENABLED` clearly; default safe for a local demo.

---

## 4. Submission steps (after conversion + clean checks)

1. Ensure default branch contains **only package code**.  
2. Open issue titled **`ClinicalVariantR`** on
   [BiocContributions](https://github.com/Bioconductor/BiocContributions/issues)
   using the official template; paste
   `https://github.com/safarafique/ClinicalVariantR`.  
3. Fix any precheck failures.  
4. Comment `/accept-policies`.  
5. Link / push to Bioconductor staging as instructed.  
6. Bump `0.99.z` for each rebuild after fixes.  
7. Respond point-by-point to the assigned reviewer.

Draft text: `ISSUE_TEMPLATE_DRAFT.md`.

---

## 5. Honest readiness assessment

| Area | Ready? |
|------|--------|
| Required doc files scaffold | **Yes** (this delivery) |
| Valid Bioconductor software package | **Not yet** — layout still Shiny-app-centric |
| Clean `R CMD check` / `BiocCheck` | **Not yet** — run after conversion |
| Ready to open BiocContributions issue | **Only after** §3 conversion + clean builds |

**Bottom line:** documentation and metadata for submission are prepared.
Complete `PACKAGE_CONVERSION.md`, then verify checks, then submit.
