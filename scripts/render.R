args <- commandArgs(trailingOnly = TRUE)
rmd <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else "NP_Crisis.Rmd"

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Missing package 'rmarkdown'. Install it (or run scripts/restore-renv.R) and retry.", call. = FALSE)
}

rmarkdown::render(rmd, output_format = "html_document")

