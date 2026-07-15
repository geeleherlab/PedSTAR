# Summarize DepMap expression and protein data for selected candidate genes.

stopifnot(requireNamespace("yaml", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))

config_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(config_file)) {
  config_file <- "config/example_config.yml"
}

config <- yaml::read_yaml(config_file)$depmap

summarize_depmap_table <- function(path, output_prefix, config) {
  depmap <- readr::read_csv(path, show_col_types = FALSE)

  grouped <- depmap |>
    dplyr::mutate(
      disease_group = dplyr::case_when(
        .data$primary_disease == "Ewing Sarcoma" ~ "EWS",
        .data$primary_disease == "Osteosarcoma" ~ "OS",
        .data$primary_disease == "Rhabdomyosarcoma" ~ "RMS",
        .data$primary_disease %in% config$brain_disease_labels ~ "Brain",
        grepl("Sarcoma|sarcoma", .data$primary_disease) ~ "Other_sarcomas",
        TRUE ~ "Others"
      )
    )

  genes <- config$genes
  missing_genes <- setdiff(genes, colnames(grouped))
  if (length(missing_genes) > 0) {
    stop(sprintf("Missing genes in DepMap table: %s", paste(missing_genes, collapse = ", ")))
  }

  summary <- grouped |>
    dplyr::group_by(.data$disease_group) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(genes),
        list(
          mean = function(x) mean(x, na.rm = TRUE),
          median = function(x) median(x, na.rm = TRUE)
        )
      ),
      n = dplyr::n(),
      .groups = "drop"
    )

  readr::write_csv(
    grouped,
    file.path(config$paths$output_dir, sprintf("%s_grouped.csv", output_prefix))
  )
  readr::write_csv(
    summary,
    file.path(config$paths$output_dir, sprintf("%s_summary.csv", output_prefix))
  )
}

dir.create(config$paths$output_dir, recursive = TRUE, showWarnings = FALSE)

summarize_depmap_table(config$paths$expression_table, "depmap_candidate_gene_expression", config)
summarize_depmap_table(config$paths$protein_table, "depmap_candidate_gene_protein", config)
