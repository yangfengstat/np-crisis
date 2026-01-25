args <- commandArgs(trailingOnly = TRUE)
rmd <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else "NP_Crisis.Rmd"

out <- if (length(args) >= 2 && nzchar(args[[2]])) {
  args[[2]]
} else if (grepl("\\.Rmd$", rmd, ignore.case = TRUE)) {
  sub("\\.Rmd$", ".knit.md", rmd, ignore.case = TRUE)
} else {
  paste0(rmd, ".knit.md")
}

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("Missing package 'knitr'. Install it (or run scripts/restore-renv.R) and retry.", call. = FALSE)
}

knitr::knit(rmd, output = out)

