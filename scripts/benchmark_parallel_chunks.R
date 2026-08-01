#!/usr/bin/env Rscript
# Micro-benchmark: sequential vs parallel chunk scoring (max 5 workers).
#
# Usage (from package root, after sourcing the Shiny stack or installing deps):
#   Rscript scripts/benchmark_parallel_chunks.R [n_variants] [chunk_size]
#
# Notes:
# - Uses synthetic missense-like rows; measures score_chunk_jobs only (not VCF I/O).
# - Parallel benefit is largest on Linux/WSL (fork). Windows PSOCK has higher startup cost.

args <- commandArgs(trailingOnly = TRUE)
n_variants <- if (length(args) >= 1L) as.integer(args[[1]]) else 5000L
chunk_size <- if (length(args) >= 2L) as.integer(args[[2]]) else 1000L

root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
shiny_root <- if (dir.exists(file.path(root, "inst/shinyapp"))) {
  file.path(root, "inst/shinyapp")
} else {
  root
}
setwd(shiny_root)

# Minimal bootstrap of scoring stack used by the app.
source("global.R", local = FALSE)

make_jobs <- function(n, size) {
  n_chunks <- as.integer(ceiling(n / size))
  lapply(seq_len(n_chunks), function(i) {
    start <- (i - 1L) * size + 1L
    end <- min(n, i * size)
    k <- end - start + 1L
    df <- data.frame(
      chrom = rep("1", k),
      pos = as.integer(seq.int(start, end)),
      ref = rep("A", k),
      alt = rep("G", k),
      gene = rep("BRCA1", k),
      consequence = rep("missense_variant", k),
      gnomad_af = rep(1e-5, k),
      revel_score = rep(0.2, k),
      stringsAsFactors = FALSE
    )
    df$variant_id <- paste(df$chrom, df$pos, df$ref, df$alt, sep = "-")
    list(df = df, id = i)
  })
}

jobs <- make_jobs(n_variants, chunk_size)
message(sprintf(
  "Benchmark: %d variants, chunk_size=%d, n_jobs=%d, RAM~%s GB, usable_cores=%d",
  n_variants, chunk_size, length(jobs),
  {
    ram <- detect_system_ram_gb()
    if (is.na(ram)) "unknown" else as.character(round(ram, 1))
  },
  detect_usable_cores()
))

run_once <- function(workers) {
  t0 <- proc.time()[["elapsed"]]
  out <- score_chunk_jobs(
    jobs = jobs,
    mode = "rapid",
    manual_inputs = list(),
    manual_by_variant = list(),
    clinical_context = NULL,
    pedigree_context = NULL,
    refs = NULL,
    session_id = "bench",
    profile_id = if (exists("DEFAULT_PROFILE_ID")) DEFAULT_PROFILE_ID else "general_germline",
    run_metadata = NULL,
    workers = workers,
    failed_log = NULL
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  list(elapsed = elapsed, n_ok = out$n_ok, n_failed = out$n_failed)
}

seq_res <- run_once(1L)
par_workers <- resolve_parallel_chunks(5L)
par_res <- run_once(par_workers)

speedup <- if (isTRUE(par_res$elapsed > 0)) seq_res$elapsed / par_res$elapsed else NA_real_
cat("\n=== ClinicalVariantR chunk scoring benchmark ===\n")
cat(sprintf("Sequential (workers=1): %.2fs  ok=%d fail=%d\n",
            seq_res$elapsed, seq_res$n_ok, seq_res$n_failed))
cat(sprintf("Parallel   (workers=%d): %.2fs  ok=%d fail=%d\n",
            par_workers, par_res$elapsed, par_res$n_ok, par_res$n_failed))
cat(sprintf("Speedup: %.2fx\n", speedup))
cat("Trade-off: higher workers => more CPU use and RAM (~workers * chunk_size rows in flight).\n")
