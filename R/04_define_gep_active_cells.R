# Define active cells for each GEP using Gaussian mixture models.

stopifnot(requireNamespace("yaml", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))

source("R/00_utils_gep_annotation.R")

config_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(config_file)) {
  config_file <- "config/example_config.yml"
}

config <- yaml::read_yaml(config_file)$gep_active_cells

split0_h <- read_usage_matrix(config$paths$split0_usages)
communities_other <- readRDS(config$paths$community_index)

activity_source <- config$parameters$activity_source
seeds <- unlist(config$parameters$seeds)
fit_strategy <- config$parameters$fit_strategy
max_tries <- config$parameters$max_tries
posterior_cutoff <- config$parameters$posterior_cutoff
if (is.null(posterior_cutoff) || is.na(posterior_cutoff)) {
  posterior_cutoff <- NULL
}

active_cells <- list()
thresholds <- list()

for (gep_id in seq_along(communities_other)) {
  program_ids <- communities_other[[gep_id]]$Split0
  activity <- get_gep_activity(split0_h, program_ids, activity_source)
  names(activity) <- rownames(split0_h)

  classified <- fit_gmm_active_cells(
    activity = activity,
    seeds = seeds,
    fit_strategy = fit_strategy,
    max_tries = max_tries,
    posterior_cutoff = posterior_cutoff
  )

  if (is.null(classified)) {
    warning(sprintf("GMM failed for GEP %s", gep_id))
    next
  }

  thresholds[[gep_id]] <- unique(classified[, c(
    "threshold",
    "active_component",
    "best_seed",
    "bic"
  )])
  thresholds[[gep_id]]$GEP <- gep_id
  thresholds[[gep_id]]$n_active_cells <- sum(classified$active)

  active_subset <- classified[classified$active, c("cell", "activity", "component", "posterior")]
  active_subset$GEP <- gep_id
  active_cells[[gep_id]] <- active_subset[, c("GEP", "cell", "activity", "component", "posterior")]

  message(sprintf("Finished GEP %s", gep_id))
}

thresholds <- do.call(rbind, thresholds)
thresholds <- thresholds[, c("GEP", "threshold", "active_component", "best_seed", "bic", "n_active_cells")]
active_cells <- do.call(rbind, active_cells)

dir.create(config$paths$output_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  thresholds,
  file.path(config$paths$output_dir, "gep_gmm_thresholds.csv")
)
readr::write_csv(
  active_cells,
  file.path(config$paths$output_dir, "gep_gmm_active_cells.csv")
)
