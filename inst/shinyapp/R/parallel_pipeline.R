# ---------------------------------------------------------------------------
# Parallel VCF scoring pipeline (production)
#
# Strategy
# --------
# 1) VCF I/O stays sequential (single reader / bcftools pipe) to avoid
#    overlapping disk seeks and keep parse order deterministic.
# 2) ACMG scoring is CPU-bound. Up to MAX_PARALLEL_CHUNKS (5) batches run
#    concurrently via parallel::mclapply (fork on Linux/WSL) or a PSOCK
#    cluster on Windows.
# 3) CSV / audit writes happen only on the main process after workers return
#    (append is not concurrent-safe).
#
# Chunk-size trade-offs
# ---------------------
# - Larger chunks => less scheduling overhead, higher peak RAM
#   (~chunk_size * parallel_workers * bytes_per_scored_row).
# - Smaller chunks => smoother progress / lower RAM, more overhead.
# - Auto sizing uses RAM + CPU cores + optional dataset size estimate.
# - Override: function args, UI, CLINICALVARIANTR_CHUNK_SIZE,
#   CLINICALVARIANTR_PARALLEL_CHUNKS.
# ---------------------------------------------------------------------------

#' Detect total system RAM in GiB (best-effort; NA if unknown).
#' @noRd
detect_system_ram_gb <- function() {
  if (file.exists("/proc/meminfo")) {
    lines <- tryCatch(readLines("/proc/meminfo", n = 5L, warn = FALSE), error = function(e) character())
    hit <- grep("^MemTotal:", lines, value = TRUE)
    if (length(hit) == 1L) {
      kb <- suppressWarnings(as.numeric(sub(".*?([0-9]+).*", "\\1", hit)))
      if (!is.na(kb) && kb > 0) return(kb / (1024^2))
    }
  }
  if (.Platform$OS.type == "windows") {
    ps <- tryCatch(
      system2(
        "powershell",
        c("-NoProfile", "-Command",
          "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,2)"),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character()
    )
    gb <- suppressWarnings(as.numeric(trimws(ps[1])))
    if (!is.na(gb) && gb > 0) return(gb)
  }
  if (nzchar(Sys.which("sysctl"))) {
    out <- tryCatch(
      system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    )
    bytes <- suppressWarnings(as.numeric(trimws(out[1])))
    if (!is.na(bytes) && bytes > 0) return(bytes / (1024^3))
  }
  NA_real_
}

#' Usable CPU workers (leave one core for OS/Shiny when possible).
#' @noRd
detect_usable_cores <- function() {
  n <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  n <- suppressWarnings(as.integer(n))
  if (is.na(n) || n < 1L) return(1L)
  max(1L, n - 1L)
}

#' Resolve concurrency: min(requested, max 5, usable cores).
#' @noRd
resolve_parallel_chunks <- function(parallel_chunks = NULL) {
  max_n <- if (exists("MAX_PARALLEL_CHUNKS", inherits = TRUE)) {
    as.integer(get("MAX_PARALLEL_CHUNKS", inherits = TRUE))
  } else {
    5L
  }
  default_n <- if (exists("DEFAULT_PARALLEL_CHUNKS", inherits = TRUE)) {
    as.integer(get("DEFAULT_PARALLEL_CHUNKS", inherits = TRUE))
  } else {
    5L
  }
  if (is.null(parallel_chunks) || length(parallel_chunks) < 1L || is.na(parallel_chunks[[1L]])) {
    env_val <- Sys.getenv("CLINICALVARIANTR_PARALLEL_CHUNKS", unset = "")
    parallel_chunks <- if (nzchar(env_val)) as.integer(env_val) else default_n
  }
  n <- suppressWarnings(as.integer(parallel_chunks[[1L]]))
  if (is.na(n) || n < 1L) n <- default_n
  as.integer(max(1L, min(n, max_n, detect_usable_cores())))
}

#' Optimal chunk size from RAM, concurrency, and optional dataset size.
#' Approximate peak RAM ~ parallel_chunks * chunk_size * ~4KB scored rows.
#' @noRd
auto_chunk_size_for_system <- function(parallel_chunks = 5L,
                                       ram_gb = NULL,
                                       n_variants = NULL) {
  env_val <- Sys.getenv("CLINICALVARIANTR_CHUNK_SIZE", unset = "")
  if (nzchar(env_val)) {
    n <- suppressWarnings(as.integer(env_val))
    if (!is.na(n) && n >= 500L) return(as.integer(n))
  }

  parallel_chunks <- max(1L, as.integer(parallel_chunks %||% 5L))
  if (is.null(ram_gb) || length(ram_gb) < 1L || is.na(ram_gb[[1L]])) {
    ram_gb <- detect_system_ram_gb()
  }
  ram_gb <- suppressWarnings(as.numeric(ram_gb[[1L]]))

  chunk <- if (is.na(ram_gb)) {
    10000L
  } else if (ram_gb < 8) {
    2500L
  } else if (ram_gb < 12) {
    5000L
  } else if (ram_gb < 20) {
    8000L
  } else if (ram_gb < 40) {
    15000L
  } else if (ram_gb < 64) {
    25000L
  } else {
    40000L
  }

  if (parallel_chunks < 5L) {
    chunk <- as.integer(min(50000L, round(chunk * (5 / parallel_chunks))))
  }

  n_variants <- suppressWarnings(as.integer(n_variants %||% NA_integer_))
  if (!is.na(n_variants) && n_variants > 0L) {
    if (n_variants < chunk) {
      chunk <- max(500L, n_variants)
    } else if (n_variants > 5e6 && !is.na(ram_gb) && ram_gb >= 32) {
      chunk <- as.integer(min(50000L, round(chunk * 1.25)))
    }
  }
  max(500L, as.integer(chunk))
}

#' Resolve chunk size (auto unless explicit override).
#' @noRd
resolve_chunk_size <- function(chunk_size = NULL,
                               parallel_chunks = 5L,
                               auto = TRUE,
                               n_variants = NULL) {
  parallel_chunks <- resolve_parallel_chunks(parallel_chunks)
  if (isTRUE(auto) || is.null(chunk_size) || length(chunk_size) < 1L ||
      is.na(chunk_size[[1L]])) {
    return(auto_chunk_size_for_system(parallel_chunks, n_variants = n_variants))
  }
  n <- suppressWarnings(as.integer(chunk_size[[1L]]))
  if (is.na(n) || n < 500L) {
    return(auto_chunk_size_for_system(parallel_chunks, n_variants = n_variants))
  }
  max(500L, n)
}

#' UI: auto chunk size + capped parallel workers.
#' @noRd
performance_tuning_ui <- function(suffix) {
  ram_gb <- detect_system_ram_gb()
  workers <- resolve_parallel_chunks(5L)
  auto_n <- auto_chunk_size_for_system(workers, ram_gb = ram_gb)
  ram_label <- if (is.na(ram_gb)) "unknown" else paste0(round(ram_gb, 1), " GB")
  cores <- detect_usable_cores()
  tagList(
    checkboxInput(
      paste0("auto_chunk_", suffix),
      "Auto chunk size from system RAM/CPU (recommended)",
      TRUE
    ),
    numericInput(
      paste0("chunk_size_", suffix),
      "Chunk size (variants per batch)",
      value = auto_n,
      min = 500,
      step = 500
    ),
    tags$p(
      class = "text-muted small mb-1",
      sprintf(
        "Parallel workers: %d (max 5; usable cores %d). Detected RAM: %s. Recommended chunk: %s.",
        workers, cores, ram_label, format(auto_n, big.mark = ",")
      )
    ),
    tags$p(
      class = "text-muted small",
      "I/O streams sequentially; up to 5 chunks score in parallel. Failed chunks are logged and skipped."
    )
  )
}

#' Append one failed-chunk CSV row (main process only).
#' @noRd
log_failed_chunk <- function(path, chunk_id, n_rows, err_msg) {
  if (is.null(path) || !nzchar(as.character(path))) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    chunk_id = as.integer(chunk_id),
    n_rows = as.integer(n_rows %||% NA_integer_),
    error = as.character(err_msg %||% ""),
    stringsAsFactors = FALSE
  )
  exists_already <- file.exists(path)
  tryCatch(
    utils::write.table(
      row,
      file = path,
      sep = ",",
      row.names = FALSE,
      col.names = !exists_already,
      append = exists_already,
      qmethod = "double"
    ),
    error = function(e) message("Unable to write failed-chunk log: ", conditionMessage(e))
  )
  invisible(NULL)
}

#' Score one chunk with isolation; never throws.
#' @noRd
score_one_chunk_safe <- function(job,
                                 mode,
                                 manual_inputs,
                                 manual_by_variant,
                                 clinical_context,
                                 pedigree_context,
                                 refs,
                                 session_id,
                                 profile_id,
                                 run_metadata) {
  chunk_id <- job$id %||% NA_integer_
  n_rows <- if (!is.null(job$df) && is.data.frame(job$df)) nrow(job$df) else NA_integer_
  tryCatch(
    {
      report <- run_acmg_pro_chunk(
        variants_df = job$df,
        mode = mode,
        manual_inputs = manual_inputs,
        manual_by_variant = manual_by_variant,
        clinical_context = clinical_context,
        pedigree_context = pedigree_context,
        refs = refs,
        session_id = session_id,
        profile_id = profile_id,
        run_metadata = run_metadata,
        write_audit = FALSE
      )
      list(ok = TRUE, chunk_id = chunk_id, report = report, error = NULL, n_rows = n_rows)
    },
    error = function(e) {
      list(
        ok = FALSE,
        chunk_id = chunk_id,
        report = NULL,
        error = conditionMessage(e),
        n_rows = n_rows
      )
    }
  )
}

#' Score a wave of chunks with controlled parallelism (max 5).
#' Returns list(results, n_ok, n_failed). CSV writes stay on the caller.
#' @noRd
score_chunk_jobs <- function(
    jobs,
    mode,
    manual_inputs,
    manual_by_variant,
    clinical_context,
    pedigree_context,
    refs,
    session_id,
    profile_id,
    run_metadata,
    workers = 1L,
    failed_log = NULL) {

  if (length(jobs) == 0L) {
    return(list(results = list(), n_ok = 0L, n_failed = 0L))
  }

  workers <- resolve_parallel_chunks(workers)
  workers <- max(1L, min(as.integer(workers), length(jobs)))

  run_safe <- function(job) {
    score_one_chunk_safe(
      job = job,
      mode = mode,
      manual_inputs = manual_inputs,
      manual_by_variant = manual_by_variant,
      clinical_context = clinical_context,
      pedigree_context = pedigree_context,
      refs = refs,
      session_id = session_id,
      profile_id = profile_id,
      run_metadata = run_metadata
    )
  }

  results <- if (workers <= 1L || length(jobs) == 1L) {
    lapply(jobs, run_safe)
  } else if (.Platform$OS.type != "windows") {
    parallel::mclapply(
      jobs,
      run_safe,
      mc.cores = workers,
      mc.preschedule = TRUE,
      mc.cleanup = TRUE
    )
  } else {
    cl <- tryCatch(parallel::makeCluster(workers), error = function(e) NULL)
    if (is.null(cl)) {
      lapply(jobs, run_safe)
    } else {
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      export_ok <- tryCatch({
        env <- .GlobalEnv
        nms <- ls(envir = env, all.names = FALSE)
        keep <- vapply(nms, function(nm) {
          if (nm %in% c("input", "output", "session", "server", "ui", "ctx")) return(FALSE)
          obj <- get(nm, envir = env, inherits = FALSE)
          is.function(obj) || is.data.frame(obj) || inherits(obj, "data.table") ||
            is.list(obj) || is.atomic(obj)
        }, logical(1))
        parallel::clusterExport(cl, nms[keep], envir = env)
        parallel::clusterExport(
          cl,
          c("mode", "manual_inputs", "manual_by_variant", "clinical_context",
            "pedigree_context", "refs", "session_id", "profile_id", "run_metadata",
            "run_safe", "score_one_chunk_safe"),
          envir = environment()
        )
        TRUE
      }, error = function(e) FALSE)
      if (isTRUE(export_ok)) {
        tryCatch(
          parallel::parLapply(cl, jobs, run_safe),
          error = function(e) lapply(jobs, run_safe)
        )
      } else {
        lapply(jobs, run_safe)
      }
    }
  }

  n_ok <- 0L
  n_failed <- 0L
  for (res in results) {
    if (isTRUE(res$ok)) {
      n_ok <- n_ok + 1L
    } else {
      n_failed <- n_failed + 1L
      msg <- res$error %||% "unknown chunk failure"
      message(sprintf(
        "Chunk %s failed (%s rows): %s — continuing",
        as.character(res$chunk_id %||% "?"),
        as.character(res$n_rows %||% "?"),
        msg
      ))
      log_failed_chunk(failed_log, res$chunk_id, res$n_rows, msg)
    }
  }
  list(results = results, n_ok = n_ok, n_failed = n_failed)
}
