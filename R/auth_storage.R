#' Authentication and encrypted storage for clinical PHI uploads.
#' @noRd

AUTH_USERS_PATH <- file.path("config", "auth_users.csv")
SECURE_UPLOAD_DIR <- file.path("logs", "secure_uploads")
AUTH_ENABLED <- FALSE
# Skip encrypting very large uploads on receive (analysis uses Shiny temp path).
SECURE_UPLOAD_ENCRYPT_MAX_MB <- as.numeric(Sys.getenv("CLINICALVARIANTR_ENCRYPT_MAX_MB", unset = "25"))
AES_GCM_IV_BYTES <- 12L

load_auth_users <- function(path = AUTH_USERS_PATH) {
  if (!file.exists(path)) return(data.frame(username = character(), password_hash = character()))
  utils::read.csv(path, stringsAsFactors = FALSE)
}

hash_password <- function(password, salt = "") {
  if (!requireNamespace("digest", quietly = TRUE)) {
    return(paste0("PLAIN:", password))
  }
  digest::digest(paste0(salt, password), algo = "sha256")
}

verify_user_password <- function(username, password, path = AUTH_USERS_PATH) {
  if (!isTRUE(AUTH_ENABLED)) return(TRUE)
  users <- load_auth_users(path)
  if (nrow(users) == 0L) {
    env_user <- Sys.getenv("CLINICALVARIANTR_USER", unset = "")
    env_pass <- Sys.getenv("CLINICALVARIANTR_PASSWORD", unset = "")
    if (nzchar(env_user) && nzchar(env_pass)) {
      return(identical(username, env_user) && identical(password, env_pass))
    }
    return(TRUE)
  }
  row <- users[tolower(users$username) == tolower(username), , drop = FALSE]
  if (nrow(row) == 0L) return(FALSE)
  expected <- row$password_hash[1]
  if (grepl("^PLAIN:", expected)) {
    return(identical(password, sub("^PLAIN:", "", expected)))
  }
  identical(hash_password(password), expected) ||
    identical(hash_password(password, salt = username), expected)
}

get_encryption_key <- function() {
  key <- Sys.getenv("CLINICALVARIANTR_ENCRYPTION_KEY", unset = "")
  if (nzchar(key)) return(charToRaw(substr(key, 1, 32)))
  if (requireNamespace("digest", quietly = TRUE)) {
    return(charToRaw(substr(digest::digest(Sys.info()[["nodename"]], algo = "sha256"), 1, 32)))
  }
  charToRaw(paste0("ClinicalVariantR-default-key-", Sys.info()[["user"]]))
}

encrypt_upload_file <- function(src_path, dest_path = NULL) {
  dir.create(SECURE_UPLOAD_DIR, recursive = TRUE, showWarnings = FALSE)
  if (is.null(dest_path)) {
    dest_path <- file.path(
      SECURE_UPLOAD_DIR,
      paste0(basename(src_path), ".", format(Sys.time(), "%Y%m%d_%H%M%S"), ".enc")
    )
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    raw_data <- readBin(src_path, "raw", file.info(src_path)$size)
    key <- get_encryption_key()
    iv <- openssl::rand_bytes(AES_GCM_IV_BYTES)
    ct <- openssl::aes_gcm_encrypt(raw_data, key, iv)
    out <- c(as.raw(iv), ct)
    writeBin(out, dest_path)
    return(list(path = dest_path, encrypted = TRUE, method = "AES-256-GCM"))
  }
  file.copy(src_path, dest_path, overwrite = TRUE)
  list(path = dest_path, encrypted = FALSE, method = "copy-with-warning")
}

decrypt_upload_file <- function(enc_path, dest_path) {
  if (!requireNamespace("openssl", quietly = TRUE)) {
    file.copy(enc_path, dest_path, overwrite = TRUE)
    return(dest_path)
  }
  raw_data <- readBin(enc_path, "raw", file.info(enc_path)$size)
  iv <- raw_data[seq_len(AES_GCM_IV_BYTES)]
  ct <- raw_data[-seq_len(AES_GCM_IV_BYTES)]
  key <- get_encryption_key()
  plain <- openssl::aes_gcm_decrypt(ct, key, iv)
  writeBin(plain, dest_path)
  dest_path
}

secure_store_shiny_upload <- function(shiny_fileinfo, label = "upload") {
  if (is.null(shiny_fileinfo) || is.null(shiny_fileinfo$datapath)) return(NULL)
  size_mb <- file.info(shiny_fileinfo$datapath)$size / (1024^2)
  if (is.finite(SECURE_UPLOAD_ENCRYPT_MAX_MB) && size_mb > SECURE_UPLOAD_ENCRYPT_MAX_MB) {
    return(list(
      label = label,
      original_name = shiny_fileinfo$name,
      secure_path = shiny_fileinfo$datapath,
      encrypted = FALSE,
      method = sprintf("temp-only (%.1f MB; encrypt skipped on upload)", size_mb),
      stored_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
    ))
  }
  result <- encrypt_upload_file(shiny_fileinfo$datapath)
  list(
    label = label,
    original_name = shiny_fileinfo$name,
    secure_path = result$path,
    encrypted = result$encrypted,
    method = result$method,
    stored_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
  )
}

append_access_audit <- function(username, action, details = "", path = file.path("logs", "access_audit.csv")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
    username = username %||% "anonymous",
    action = action,
    details = details,
    stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    utils::write.csv(row, path, row.names = FALSE)
  } else {
    utils::write.table(row, path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, quote = TRUE)
  }
  invisible(path)
}
