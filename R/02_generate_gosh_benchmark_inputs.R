# Generate GOSH benchmark GEP input files.

stopifnot(requireNamespace("yaml", quietly = TRUE))

source("R/00_utils_gep_annotation.R")

config_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(config_file)) {
  config_file <- "config/example_config.yml"
}

config <- yaml::read_yaml(config_file)$gosh_benchmark

comms <- readRDS(config$paths$gep_loadings)
if (!is.null(config$dataset_key)) {
  comms <- comms[[config$dataset_key]]
}

communities_other <- readRDS(config$paths$community_index)
split0_h <- read_usage_matrix(config$paths$split0_usages)
singler_scores <- readRDS(config$paths$single_r_label_scores)
singler_delta <- readRDS(config$paths$single_r_delta_values)
gene_set_databases <- load_gene_sets(config$gene_sets)

cell_ids <- sub(".*?(.*).h5ad", "\\1", rownames(singler_delta))
rownames(singler_delta) <- cell_ids
rownames(singler_scores) <- cell_ids

output_dir <- config$paths$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (gep_id in seq_len(ncol(comms))) {
  gep_vector <- comms[, gep_id]
  program_ids <- communities_other[[gep_id]]$Split0

  top_cells <- get_top_cells(
    split0_h = split0_h,
    program_ids = program_ids,
    top_cells = config$parameters$top_cells
  )

  enrichment_list <- summarize_enrichments(
    gep_vector = gep_vector,
    gene_set_databases = gene_set_databases,
    top_n = config$parameters$top_genes,
    n_terms = config$parameters$top_enrichment_terms
  )

  mean_scores <- summarize_single_r_scores(singler_scores, top_cells)

  write_gep_prompt_input(
    gep_vector = gep_vector,
    gep_id = gep_id,
    file = file.path(output_dir, sprintf("gosh_gep_%s_input.txt", gep_id)),
    top_genes = config$parameters$top_genes,
    enrichment_list = enrichment_list,
    top_cell_count = length(top_cells),
    mean_scores = mean_scores
  )

  message(sprintf("Wrote GOSH GEP %s", gep_id))
}
