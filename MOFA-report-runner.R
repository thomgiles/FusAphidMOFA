#!/usr/bin/env Rscript

## Render the dynamic UoN DRS MOFA HTML report.
##
## Examples:
##   Rscript MOFA-report-runner.R
##   Rscript MOFA-report-runner.R --run-analysis
##   Rscript MOFA-report-runner.R --output MOFA-report.html

args <- commandArgs(trailingOnly = TRUE)

has_flag <- function(flag) {
  flag %in% args
}

value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

if (has_flag("--help") || has_flag("-h")) {
  cat(
    "Usage: Rscript MOFA-report-runner.R [options]\n\n",
    "Options:\n",
    "  --run-analysis       Run MOFA-CORE.R workflow functions before rendering the report.\n",
    "  --force              Rebuild output steps even when expected files already exist.\n",
    "  --new-normalisation  Use the newer normalisation path instead of the publication default.\n",
    "  --input-file FILE    Workbook to use when --run-analysis is supplied.\n",
    "  --full-tables        Render full CSV/TSV table contents, split into 50-row tabs.\n",
    "  --all-images         Include every image in the report.\n",
    "  --max-images N       Include at most N images in the report. Default: all images.\n",
    "  --output FILE        Output HTML file. Default: MOFA-report.html.\n",
    "  --help, -h           Show this help text.\n",
    sep = ""
  )
  quit(status = 0)
}

run_analysis <- has_flag("--run-analysis")
skip_existing <- !has_flag("--force")
old_normalisation <- !has_flag("--new-normalisation")
full_tables <- has_flag("--full-tables")
input_file <- value_after("--input-file", "")
output_file <- value_after("--output", "MOFA-report.html")

max_images_arg <- value_after("--max-images")
if (has_flag("--all-images") || is.null(max_images_arg)) {
  max_images <- Inf
} else {
  max_images <- suppressWarnings(as.numeric(max_images_arg))
  if (!is.finite(max_images) || max_images <= 0) {
    stop("--max-images must be a positive number, or use --all-images.", call. = FALSE)
  }
}

report_packages <- c("rmarkdown", "knitr")
missing_report_packages <- report_packages[
  !vapply(report_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_report_packages) > 0) {
  stop(
    "Missing report package(s) in the current R library: ",
    paste(missing_report_packages, collapse = ", "),
    ". Install them before rendering the report.",
    call. = FALSE
  )
}

if (!run_analysis) {
  message("Rendering MOFA report...")
  message("  run_analysis:  ", run_analysis)
  message("  skip_existing: ", skip_existing)
  message("  old_normalisation: ", old_normalisation)
  message("  input_file:    ", if (nzchar(input_file)) input_file else "<auto>")
  message("  full_tables:   ", full_tables)
  message("  max_images:    ", if (is.infinite(max_images)) "Inf" else max_images)
  message("  output_file:   ", output_file)
}

rmarkdown::render(
  input = "MOFA-report.Rmd",
  output_file = output_file,
  params = list(
    run_analysis = run_analysis,
    skip_existing = skip_existing,
    input_workbook = input_file,
    old_normalisation = old_normalisation,
    full_tables = full_tables,
    max_images = max_images
  ),
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

if (!run_analysis) {
  message("Report written to: ", output_file)
}
