# Generate one GEP annotation input file per program.

stopifnot(requireNamespace("yaml", quietly = TRUE))

source("R/00_utils_gep_annotation.R")

config_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(config_file)) {
  config_file <- "config/example_config.yml"
}

config <- yaml::read_yaml(config_file)

comms <- readRDS(config$paths$gep_loadings)
metadata <- readRDS(config$paths$metadata)
communities_other <- readRDS(config$paths$community_index)
split0_h <- readRDS(config$paths$split0_usages)
singler_scores <- readRDS(config$paths$single_r_label_scores)
singler_delta <- readRDS(config$paths$single_r_delta_values)

numbat_metadata <- NULL
if (!is.null(config$paths$numbat_metadata) && file.exists(config$paths$numbat_metadata)) {
  numbat_metadata <- readRDS(config$paths$numbat_metadata)
}

gene_set_databases <- load_gene_sets(config$gene_sets)

cell_ids <- sub(".*?(.*).h5ad", "\\1", rownames(singler_delta))
rownames(singler_delta) <- cell_ids
rownames(singler_scores) <- cell_ids

output_dir <- config$paths$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_geps <- ncol(comms)
for (gep_id in seq_len(n_geps)) {
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
  disease_freq <- summarize_disease_origin(top_cells, metadata)
  numbat_scores <- summarize_numbat_scores(top_cells, numbat_metadata)

  write_gep_prompt_input(
    gep_vector = gep_vector,
    gep_id = gep_id,
    file = file.path(output_dir, sprintf("gep_%s_input.txt", gep_id)),
    top_genes = config$parameters$top_genes,
    enrichment_list = enrichment_list,
    top_cell_count = length(top_cells),
    mean_scores = mean_scores,
    disease_freq = disease_freq,
    numbat_scores = numbat_scores
  )

  message(sprintf("Wrote GEP %s", gep_id))
}
