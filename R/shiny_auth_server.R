#' Login modal and authentication observers.
#' @noRd
register_auth_server <- function(ctx) {
  input <- ctx$input
  session <- ctx$session
  authorized <- ctx$authorized
  auth_user <- ctx$auth_user

  observe({
    if (!isTRUE(AUTH_ENABLED) || authorized()) return()
    showModal(modalDialog(
      title = tagList(icon("lock"), " ClinicalVariantR Login"),
      p("Sign in to upload clinical data. Uploads are encrypted at rest when OpenSSL is available."),
      p(class = "text-muted small", "Use configured credentials. Select files after sign-in, or re-select if you chose a file before login."),
      textInput("login_user", "Username", placeholder = "admin"),
      passwordInput("login_pass", "Password"),
      footer = actionButton("login_btn", "Sign in", class = "btn-primary"),
      easyClose = FALSE,
      fade = FALSE
    ))
  })

  observeEvent(input$login_btn, {
    user <- trimws(input$login_user %||% "")
    pass <- input$login_pass %||% ""
    if (verify_user_password(user, pass)) {
      authorized(TRUE)
      auth_user(user)
      append_access_audit(user, "login")
      removeModal()
      showNotification(paste("Signed in as", user), type = "message")
      ctx$reprocess_pending_uploads()
    } else {
      showNotification("Invalid username or password", type = "error")
    }
  })

  invisible(ctx)
}
