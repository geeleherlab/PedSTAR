# Calculate pairwise Spearman co-expression between GEP activity profiles.

stopifnot(requireNamespace("yaml", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))

source("R/00_utils_gep_annotation.R")

config_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(config_file)) {
  config_file <- "config/example_config.yml"
}

config <- yaml::read_yaml(config_file)$gep_coexpression

split0_h <- read_usage_matrix(config$paths$split0_usages)
communities_other <- readRDS(config$paths$community_index)
activity_source <- config$parameters$activity_source
if (is.null(activity_source)) {
  activity_source <- "first_split0_program"
}

n_geps <- length(communities_other)
gep_activity <- matrix(NA_real_, nrow = nrow(split0_h), ncol = n_geps)
rownames(gep_activity) <- rownames(split0_h)
colnames(gep_activity) <- paste0("GEP", seq_len(n_geps))

for (gep_id in seq_len(n_geps)) {
  program_ids <- communities_other[[gep_id]]$Split0

  if (activity_source == "first_split0_program") {
    gep_activity[, gep_id] <- split0_h[, program_ids[1]]
  } else if (activity_source == "mean_split0_programs") {
    gep_activity[, gep_id] <- rowMeans(split0_h[, program_ids, drop = FALSE], na.rm = TRUE)
  } else {
    stop("Unsupported activity_source. Use 'first_split0_program' or 'mean_split0_programs'.")
  }
}

results <- vector("list", n_geps * n_geps)
idx <- 1

for (gep1 in seq_len(n_geps)) {
  pvalues <- numeric(n_geps)
  rhos <- numeric(n_geps)

  for (gep2 in seq_len(n_geps)) {
    test <- cor.test(
      gep_activity[, gep1],
      gep_activity[, gep2],
      method = "spearman",
      exact = FALSE
    )

    rhos[gep2] <- unname(test$estimate)
    pvalues[gep2] <- test$p.value
  }

  for (gep2 in seq_len(n_geps)) {
    results[[idx]] <- data.frame(
      GEP1 = gep1,
      GEP2 = gep2,
      Rho_spearman = rhos[gep2],
      Pvalues_spearman = pvalues[gep2],
      Padj_spearman = NA_real_
    )
    idx <- idx + 1
  }

  message(sprintf("Finished GEP %s of %s", gep1, n_geps))
}

coexpression <- do.call(rbind, results)
coexpression$Padj_spearman <- p.adjust(coexpression$Pvalues_spearman, method = "BH")

coexpression_clean <- coexpression[coexpression$GEP1 != coexpression$GEP2, ]
coexpression_clean$pair_id <- paste0(
  pmin(coexpression_clean$GEP1, coexpression_clean$GEP2),
  "_",
  pmax(coexpression_clean$GEP1, coexpression_clean$GEP2)
)
coexpression_clean <- coexpression_clean[!duplicated(coexpression_clean$pair_id), ]
coexpression_clean <- coexpression_clean[order(coexpression_clean$Rho_spearman, decreasing = TRUE), ]
coexpression_clean$pair_id <- NULL
coexpression_clean$Padj_spearman <- p.adjust(coexpression_clean$Pvalues_spearman, method = "BH")

dir.create(config$paths$output_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  coexpression,
  file.path(config$paths$output_dir, "gep_spearman_coexpression_all_pairs.csv")
)
readr::write_csv(
  coexpression_clean,
  file.path(config$paths$output_dir, "gep_spearman_coexpression_unique_pairs.csv")
)
