# ClinicalVariantR Shiny module map

Edit **one file per feature** instead of the full `server.R` / `ui.R`.

## Quick reference — where to change what

| You want to change… | Edit this file |
|---------------------|----------------|
| **Group B upload / validation / Run Analysis** | `R/shiny/group_b_server.R` |
| **Group B layout (inputs, buttons, panels)** | `R/shiny/ui/group_b_ui.R` |
| Group A clinical workflow | `R/shiny/group_a_server.R` + `R/shiny/ui/group_a_ui.R` |
| Group C gene panel | `R/shiny/group_c_server.R` + `R/shiny/ui/group_c_ui.R` |
| Login / sign-in modal | `R/shiny/auth_server.R` |
| VCF upload helpers (all groups) | `R/shiny/upload_server.R` |
| Analysis pipeline runner | `R/shiny/analysis_server.R` |
| Results tables, status, downloads (A/B/C) | `R/shiny/results_server.R` |
| Evidence Explorer tab | `R/shiny/explorer_server.R` + `R/shiny/ui/explorer_ui.R` |
| Audit log tab | `R/shiny/audit_server.R` + `R/shiny/ui/audit_ui.R` |
| Navbar, session label, shared helpers | `R/shiny/shared_server.R` |
| Reactive state (new `reactiveVal`) | `R/shiny/context.R` |
| Theme, landing cards | `R/shiny/ui/theme.R` |
| Home page | `R/shiny/ui/home_ui.R` |
| VCF validation rules | `R/vcf_validate.R` |
| ACMG scoring engine | `R/acmg_engine.R`, `R/acmg_pipeline.R` |
| CSS (readiness badge, run button) | `www/custom.css` |

## Entry points (do not put logic here)

- `app.R` — launches the app
- `server.R` — registers modules in order (~15 lines)
- `ui.R` — composes nav panels (~18 lines)
- `global.R` — loads R modules + Shiny modules

## Server registration order

```
init_shiny_context
  → shared_server (helpers, navbar)
  → upload_server (load_vcf_on_upload, defer_secure_upload)
  → analysis_server (run_complete_analysis)
  → audit_server (load_audit)
  → group_a_server (sets refresh_group_a_validation)
  → group_b_server
  → group_c_server
  → results_server
  → explorer_server
  → auth_server (login; calls reprocess_pending_uploads)
```

## Module pattern

Each `register_*_server(ctx)` receives a shared context `ctx` with:

- `ctx$input`, `ctx$output`, `ctx$session`
- Reactive values: `ctx$vcf_validation_b`, `ctx$report_b_full`, etc.
- Shared functions: `ctx$load_vcf_on_upload`, `ctx$run_complete_analysis`, …

Add a new tab: create `R/shiny/my_tab_server.R`, `R/shiny/ui/my_tab_ui.R`, source in `global.R`, call `register_my_tab_server(ctx)` from `server.R`.
