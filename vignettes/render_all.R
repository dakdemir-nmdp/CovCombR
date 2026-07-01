files <- list.files(pattern = "\\.Rmd$")

for (f in files) {
    cat(sprintf("\nRendering %s...\n", f))

    # Render using each vignette's own declared output format(s) in its YAML
    # header, matching what R CMD build / devtools::build_vignettes() does.
    tryCatch(
        {
            rmarkdown::render(f, quiet = TRUE)
            cat(sprintf("  - Success\n"))
        },
        error = function(e) {
            cat(sprintf("  - Failed (%s)\n", e$message))
        }
    )
}
