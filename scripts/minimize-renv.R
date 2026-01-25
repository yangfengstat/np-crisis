if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

project <- getwd()

# Packages used directly in the analysis (incl. foreach worker packages).
required_pkgs <- c(
  "data.table", "doParallel", "dplyr", "e1071", "foreach", "foreign", "ggplot2",
  "ggpubr", "imputeMissings", "nproc", "randomGLM", "readr", "readxl",
  "rfUtilities", "ROCR", "ROSE", "rpart", "tidyr", "zoo"
)

# Prefer copying from the user's global library (no downloads).
renv::hydrate(project = project, packages = required_pkgs, prompt = FALSE)

# Remove packages not required by this project (best-effort).
renv::clean(project = project, prompt = FALSE)

# Rebuild a minimal lockfile based on usage, plus the required_pkgs above.
renv::snapshot(project = project, packages = required_pkgs, prompt = FALSE)
