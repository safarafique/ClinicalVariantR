# ACMGamp — thin server entry (logic lives in R/shiny/*_server.R)

server <- function(input, output, session) {
  ctx <- init_shiny_context(input, output, session)

  register_shared_server(ctx)
  register_upload_server(ctx)
  register_analysis_server(ctx)
  register_audit_server(ctx)

  register_group_a_server(ctx)
  register_group_b_server(ctx)
  register_group_c_server(ctx)

  register_results_server(ctx)
  register_explorer_server(ctx)
  register_auth_server(ctx)
}
