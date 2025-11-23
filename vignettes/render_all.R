files <- c(
    "CovCombR-vignette.Rmd",
    "iris-example.Rmd",
    "combining-grms.Rmd",
    "power-analysis-3x3.Rmd",
    "statistical-methods.Rmd",
    "PowerAnalysis.Rmd"
)

for (f in files) {
    cat(sprintf("\nRendering %s...\n", f))

    # Render to HTML
    tryCatch(
        {
            rmarkdown::render(f, output_format = "html_document", quiet = TRUE)
            cat(sprintf("  - HTML: Success\n"))
        },
        error = function(e) {
            cat(sprintf("  - HTML: Failed (%s)\n", e$message))
        }
    )

    # Render to PDF
    tryCatch(
        {
            rmarkdown::render(f, output_format = "pdf_document", quiet = TRUE)
            cat(sprintf("  - PDF: Success\n"))
        },
        error = function(e) {
            cat(sprintf("  - PDF: Failed (%s)\n", e$message))
        }
    )
}
