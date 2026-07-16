\name{ClinicalVariantR}
\alias{ClinicalVariantR}
\alias{ClinicalVariantRApp}
\title{Build the ClinicalVariantR Shiny application}
\description{
Returns a Shiny app object for ACMG/AMP germline variant classification.
Following Bioconductor Shiny guidelines, this function \strong{returns} the app
and does not launch it. Call \code{shiny::runApp()} yourself.
}
\usage{
ClinicalVariantR(...)

ClinicalVariantRApp(...)
}
\arguments{
\item{...}{Reserved for future options (currently ignored).}
}
\value{
A \code{shiny.appobj} created by \code{shiny::shinyApp()}.
}
\examples{
\dontrun{
app <- ClinicalVariantR()
if (interactive()) {
    shiny::runApp(app)
}
}
}
\seealso{
\code{\link{ClinicalVariantRApp}}
}
