#' Expert review checklist UI and filtered report exports.
#' @noRd

parse_criteria_met <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(as.character(x)))) return(character())
  parts <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  parts[nzchar(parts)]
}

is_pm2_only_vus <- function(classification, criteria_met) {
  if (!identical(as.character(classification), "VUS")) return(FALSE)
  met <- parse_criteria_met(criteria_met)
  length(met) == 1L && toupper(met[[1L]]) == "PM2"
}

filter_report_pathogenic_plus <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  tiers <- c("Pathogenic", "Likely Pathogenic")
  out <- df[df$classification %in% tiers, , drop = FALSE]
  rownames(out) <- NULL
  out
}

filter_report_expert_worklist <- function(
    df,
    exclude_pm2_only_vus = TRUE,
    include_vus_with_pathogenic_evidence = TRUE) {
  if (is.null(df) || nrow(df) == 0L) return(df)

  path_tiers <- c("Pathogenic", "Likely Pathogenic")
  keep <- df$classification %in% path_tiers

  if (isTRUE(include_vus_with_pathogenic_evidence)) {
    path_n <- if ("pathogenic_evidence_count" %in% names(df)) {
      vapply(df$pathogenic_evidence_count, scalar_int, integer(1L))
    } else {
      rep(0L, nrow(df))
    }
    keep <- keep | (df$classification == "VUS" & !is.na(path_n) & path_n >= 2L)
  }

  out <- df[keep, , drop = FALSE]
  if (nrow(out) == 0L) return(out)

  if (isTRUE(exclude_pm2_only_vus)) {
    drop <- vapply(seq_len(nrow(out)), function(i) {
      is_pm2_only_vus(out$classification[i], out$criteria_met[i])
    }, logical(1L))
    out <- out[!drop, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

expert_worklist_stats <- function(df) {
  empty <- list(
    total = 0L, pathogenic = 0L, likely_pathogenic = 0L,
    lp_plus = 0L, vus = 0L, pm2_only_vus = 0L, worklist = 0L
  )
  if (is.null(df) || nrow(df) == 0L) return(empty)

  pm2_only <- sum(vapply(seq_len(nrow(df)), function(i) {
    is_pm2_only_vus(df$classification[i], df$criteria_met[i])
  }, logical(1L)))

  list(
    total = nrow(df),
    pathogenic = sum(df$classification == "Pathogenic", na.rm = TRUE),
    likely_pathogenic = sum(df$classification == "Likely Pathogenic", na.rm = TRUE),
    lp_plus = sum(df$classification %in% c("Pathogenic", "Likely Pathogenic"), na.rm = TRUE),
    vus = sum(df$classification == "VUS", na.rm = TRUE),
    pm2_only_vus = pm2_only,
    worklist = nrow(filter_report_expert_worklist(df))
  )
}

expert_review_checklist_static_ui <- function(show_curation = FALSE) {
  curation_block <- if (isTRUE(show_curation)) {
    tags$li(
      tags$strong("Group A curation:"),
      " apply PS2/PS3/PP1/PP4/PS4/PM6/PP2 checkboxes only after external proof, then ",
      tags$em("Apply curation & reclassify.")
    )
  } else {
    tags$li(
      tags$strong("Manual criteria:"),
      " PS2, PS3, PP1, PP4, PM3, PM6, BS3, BS4, BP2, BP5 require off-platform or Group A curation."
    )
  }

  tagList(
    div(
      class = "alert alert-info py-2 small mb-3",
      icon("user-md"),
      " Predictions are a ",
      tags$strong("worklist"),
      ", not a signed clinical report. Use this checklist before reporting any variant."
    ),
    tags$ol(
      class = "small mb-3 ps-3",
      tags$li(tags$strong("Triage:"), " start with Pathogenic / Likely Pathogenic; review VUS with >=2 pathogenic evidence points."),
      tags$li(tags$strong("Technical QC:"), " confirm in IGV - coverage, zygosity, mapping, genome build, canonical transcript."),
      tags$li(tags$strong("Population:"), " verify gnomAD/ClinVar independently; do not rely on PM2 when AF is missing."),
      tags$li(tags$strong("ACMG evidence:"), " open evidence detail and confirm each triggered criterion (PVS1, PS1, PM5, PP3, PP5/BP6, etc.)."),
      tags$li(tags$strong("Clinical fit:"), " phenotype matches gene (PP4), inheritance mode, gene on panel/indication."),
      tags$li(tags$strong("Limitations:"), " read ", tags$code("prediction_limitations"), " - resolve or document before sign-out."),
      tags$li(tags$strong("Exclude:"), " Benign/Likely Benign and PM2-only VUS unless lab policy says otherwise."),
      curation_block,
      tags$li(tags$strong("Sign-out:"), " second reviewer per laboratory SOP; document final class and report text.")
    )
  )
}

#' Prominent handoff banner: CSV export is the best input for expert review.
#' @noRd
expert_csv_handoff_ui <- function() {
  div(
    class = "alert alert-success border mb-3",
    h6(class = "alert-heading mb-2", icon("file-csv"), " Best input for expert review: ClinicalVariantR CSV export"),
    tags$ol(
      class = "small mb-2 ps-3",
      tags$li("Run ", tags$strong("Run Analysis"), " and wait for results."),
      tags$li(
        "Download ",
        tags$strong("Export expert worklist (CSV)"),
        " to start sign-out review - or ",
        tags$strong("Download full prediction report (CSV)"),
        " for audit and LIMS import."
      ),
      tags$li(
        "Open the CSV in Excel, R, or your LIMS. Key columns: ",
        tags$code("variant_id"), ", ",
        tags$code("classification"), ", ",
        tags$code("criteria_met"), ", ",
        tags$code("prediction_limitations"), ", ",
        tags$code("evidence_json"), "."
      )
    ),
    tags$p(
      class = "small text-muted mb-0",
      "Filename pattern: ",
      tags$code("ClinicalVariantR_Expert_Worklist_*.csv"),
      " or ",
      tags$code("ClinicalVariantR_Prediction_Report_*.csv"),
      "."
    )
  )
}

render_expert_worklist_stats_ui <- function(df) {
  stats <- expert_worklist_stats(df)
  if (stats$total == 0L) {
    return(div(class = "text-muted small", "Run analysis to see expert worklist counts."))
  }

  div(
    class = "border rounded p-2 mb-3 bg-light small",
    tags$strong("Worklist counts"),
    tags$ul(
      class = "mb-0 ps-3",
      tags$li(format(stats$total, big.mark = ","), " variants analyzed"),
      tags$li(
        tags$span(class = "text-danger", stats$lp_plus, " Pathogenic / Likely Pathogenic"),
        " (priority review)"
      ),
      tags$li(stats$vus, " VUS (", stats$pm2_only_vus, " PM2-only - usually deprioritize)"),
      tags$li(
        tags$strong(stats$worklist), " expert worklist candidates ",
        "(LP+ plus VUS with >=2 pathogenic evidence; PM2-only VUS excluded)"
      )
    )
  )
}

write_expert_report_csv <- function(df, file, filter_fn) {
  filtered <- filter_fn(df)
  utils::write.csv(filtered, file, row.names = FALSE)
  invisible(nrow(filtered))
}

register_expert_worklist_downloads <- function(
    output,
    suffix,
    report_full_reactive,
    session) {

  lp_id <- paste0("download_", suffix, "_lp_plus")
  wl_id <- paste0("download_", suffix, "_worklist")

  output[[lp_id]] <- downloadHandler(
    filename = function() {
      paste0(
        "ClinicalVariantR_LP_Plus_", toupper(suffix), "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"
      )
    },
    content = function(file) {
      df <- report_full_reactive()
      req(!is.null(df), nrow(df) > 0L)
      n <- write_expert_report_csv(df, file, filter_report_pathogenic_plus)
      if (n == 0L) {
        showNotification("No Pathogenic / Likely Pathogenic variants to export.", type = "warning", session = session)
      }
    }
  )

  output[[wl_id]] <- downloadHandler(
    filename = function() {
      paste0(
        "ClinicalVariantR_Expert_Worklist_", toupper(suffix), "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"
      )
    },
    content = function(file) {
      df <- report_full_reactive()
      req(!is.null(df), nrow(df) > 0L)
      n <- write_expert_report_csv(df, file, filter_report_expert_worklist)
      if (n == 0L) {
        showNotification("No expert worklist candidates matched filters.", type = "warning", session = session)
      }
    }
  )
}
