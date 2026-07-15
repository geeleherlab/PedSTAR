# Score GOSH benchmark label-matching outputs.

stopifnot(requireNamespace("yaml", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))

source("R/00_utils_gep_annotation.R")

config_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(config_file)) {
  config_file <- "config/example_config.yml"
}

config <- yaml::read_yaml(config_file)$gosh_benchmark

truth <- readr::read_csv(config$paths$human_labels, show_col_types = FALSE)
predictions <- readr::read_csv(config$paths$label_matching_predictions, show_col_types = FALSE)

accuracy_by_attempt <- score_label_matching_attempts(predictions, truth)
summary <- data.frame(
  mean_accuracy = mean(accuracy_by_attempt$accuracy),
  sd = sd(accuracy_by_attempt$accuracy),
  n = nrow(accuracy_by_attempt)
)
summary$se <- summary$sd / sqrt(summary$n)
summary$ci_low <- summary$mean_accuracy - 1.96 * summary$se
summary$ci_high <- summary$mean_accuracy + 1.96 * summary$se

dir.create(config$paths$results_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  accuracy_by_attempt,
  file.path(config$paths$results_dir, "gosh_accuracy_by_attempt.csv")
)
readr::write_csv(
  summary,
  file.path(config$paths$results_dir, "gosh_accuracy_summary.csv")
)

print(summary)
