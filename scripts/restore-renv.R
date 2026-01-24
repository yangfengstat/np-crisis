if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::restore(prompt = FALSE)

