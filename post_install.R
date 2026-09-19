# Three packages that environment.yaml cannot supply. Run inside the activated
# rna-thermo environment:
#
#   conda activate rna-thermo
#   Rscript post_install.R
#
# 1. canvasXpress   - not packaged on conda-forge at all. Pure R plus bundled
#                     JavaScript assets, so it installs from CRAN without a
#                     toolchain.
# 2. OmnipathR      - bioconda has no build for R 4.4: versions jump from
#                     3.10.1 (r43) straight to 3.18.4 (r45), skipping the
#                     Bioconductor 3.20 release (3.14.x). Installed from
#                     Bioconductor instead.
# 3. GSVA           - bioconda has no osx-arm64 build at any version (linux-64,
#                     linux-aarch64 and osx-64 only), so it compiles from source.
# 4. clinfun        - conda-forge has no osx-arm64 build either (linux-64,
#                     osx-64, win-64 only). NeedsCompilation: yes, including
#                     Fortran, hence gfortran in environment.yaml.
#
# All four are used only by R_script_Fig_3.R and Cross_validation_MFE.R. If you
# are not running those two scripts you can skip this file entirely -- the
# Figure 1, 2, 4 and 5 scripts need nothing from here.

options(repos = c(CRAN = "https://cloud.r-project.org"))

message("R version in use: ", getRversion())
if (getRversion() < "4.4.0" || getRversion() >= "4.5.0") {
  warning(
    "Expected R 4.4.x (environment.yaml pins 4.4.2). Bioconductor 3.20 ",
    "packages are built against R 4.4; other versions may misbehave."
  )
}

# Attempt each install independently. A failure is reported and recorded but
# does not abort the run: the packages here are unrelated to each other, and
# stopping on the first problem hides the state of the rest.
failures <- character()

install_if_missing <- function(pkg, installer) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, ": already present (", utils::packageVersion(pkg), ")")
    return(invisible(TRUE))
  }
  message(pkg, ": installing")
  ok <- tryCatch({ installer(); TRUE }, error = function(e) {
    message(pkg, ": install threw - ", conditionMessage(e)); FALSE
  })
  if (!ok || !requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, ": FAILED")
    failures <<- c(failures, pkg)
    return(invisible(FALSE))
  }
  message(pkg, ": installed (", utils::packageVersion(pkg), ")")
}

install_if_missing("canvasXpress", function() {
  install.packages("canvasXpress")
})

install_if_missing("OmnipathR", function() {
  BiocManager::install("OmnipathR", version = "3.20", ask = FALSE, update = FALSE)
})

install_if_missing("GSVA", function() {
  BiocManager::install("GSVA", version = "3.20", ask = FALSE, update = FALSE)
})

install_if_missing("clinfun", function() {
  install.packages("clinfun")
})

message("\nVerifying every package the scripts reference:")
required <- c(
  # Bioconductor
  "biomaRt", "clusterProfiler", "decoupleR", "derfinder", "DESeq2",
  "enrichplot", "GenomicRanges", "GEOquery", "GSVA", "IRanges", "pathview",
  "rtracklayer", "SummarizedExperiment", "OmnipathR",
  # CRAN
  "data.table", "dplyr", "ggplot2", "ggpubr", "ggrepel", "GGally", "gplots",
  "gridExtra", "patchwork", "pheatmap", "RColorBrewer", "scales", "purrr",
  "readr", "rlang", "stringr", "tibble", "tidyr", "tidyverse", "clinfun",
  "enrichR", "factoextra", "httr", "jsonlite", "canvasXpress"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("MISSING (", length(missing), "): ", paste(missing, collapse = ", "))
  quit(status = 1)
}
message("All ", length(required), " packages available.")
