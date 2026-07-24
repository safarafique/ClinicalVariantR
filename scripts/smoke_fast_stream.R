# Smoke test for buffered VCF streaming optimizations.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
root <- normalizePath(file.path(script_dir, ".."))
setwd(root)
source("global.R", local = FALSE)

stopifnot(exists("VCF_LINE_BUFFER"), identical(VCF_LINE_BUFFER, 50000L))

tiny <- list(
  data.frame(a = 1L, b = "x", stringsAsFactors = FALSE),
  data.frame(a = 2L, b = "y", stringsAsFactors = FALSE)
)
bound <- rbind_parsed_rows(tiny)
stopifnot(nrow(bound) == 2L, identical(bound$a, c(1L, 2L)))

vcf <- "data/samples/example_variants.vcf"
stopifnot(file.exists(vcf))

n <- 0L
t0 <- Sys.time()
st <- stream_vcf_chunks(
  vcf,
  chunk_size = 50L,
  processor = function(df, id) {
    n <<- n + nrow(df)
  }
)
cat(
  "stream rows=", st$rows_analyzed, " n=", n,
  " secs=", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3),
  "\n", sep = ""
)
stopifnot(identical(st$rows_analyzed, n), n > 0L)

out <- tempfile(fileext = ".csv")
t1 <- Sys.time()
res <- analyze_complete_vcf(
  vcf,
  mode = "rapid",
  output_csv = out,
  chunk_size = 50L,
  write_audit = FALSE,
  use_bcftools = FALSE
)
cat(
  "classified=", res$rows_classified,
  " secs=", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 3),
  "\n", sep = ""
)
stopifnot(res$rows_classified > 0L, file.exists(out))
rep <- load_full_analysis_report(out)
stopifnot(!is.null(rep), nrow(rep) == res$rows_classified)

cat("OK\n")
