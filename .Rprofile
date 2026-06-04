local({
  args <- commandArgs(trailingOnly = TRUE)
  startup_args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", startup_args, value = TRUE)
  script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else ""
  is_report_runner <- identical(basename(script_file), "MOFA-report-runner.R")
  activate_renv <- !is_report_runner || "--run-analysis" %in% args

  if (activate_renv && file.exists("renv/activate.R")) {
    ..md5.. <- NULL
    source("renv/activate.R")
    rm("..md5..")
  }
})
