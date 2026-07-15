# Summarize GTEx expression for selected candidate genes.

stopifnot(requireNamespace("yaml", quietly = TRUE))
stopifnot(requireNamespace("data.table", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("tidyr", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))

config_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(config_file)) {
  config_file <- "config/example_config.yml"
}

config <- yaml::read_yaml(config_file)$gtex_expression

expr <- data.table::fread(config$paths$gene_tpm_gct, skip = 2)
anno <- read.delim(config$paths$sample_attributes)

genes_keep <- config$genes

expr_sub <- expr |>
  dplyr::filter(.data$Description %in% genes_keep)

df_long <- expr_sub |>
  tidyr::pivot_longer(
    cols = -dplyr::all_of(c("Name", "Description")),
    names_to = "SAMPID",
    values_to = "TPM"
  ) |>
  dplyr::left_join(
    anno |> dplyr::select(dplyr::all_of(c("SAMPID", "SMTSD"))),
    by = "SAMPID"
  )

tissue_summary <- df_long |>
  dplyr::group_by(.data$Description, .data$SMTSD) |>
  dplyr::summarise(
    median_TPM = median(.data$TPM, na.rm = TRUE),
    mean_TPM = mean(.data$TPM, na.rm = TRUE),
    .groups = "drop"
  )

tissue_matrix <- tissue_summary |>
  dplyr::select(dplyr::all_of(c("Description", "SMTSD", "median_TPM"))) |>
  tidyr::pivot_wider(names_from = .data$SMTSD, values_from = .data$median_TPM)

dir.create(config$paths$output_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  tissue_summary,
  file.path(config$paths$output_dir, "gtex_candidate_gene_tissue_summary.csv")
)
readr::write_csv(
  tissue_matrix,
  file.path(config$paths$output_dir, "gtex_candidate_gene_median_tpm_matrix.csv")
)
