# Helper functions for generating GEP annotation prompt inputs.

load_gene_sets <- function(paths) {
  stopifnot(requireNamespace("GSA", quietly = TRUE))

  lapply(paths, GSA::GSA.read.gmt)
}

read_usage_matrix <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }

  stopifnot(requireNamespace("vroom", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))

  usage <- vroom::vroom(path, show_col_types = FALSE)
  tibble::column_to_rownames(as.data.frame(usage), 1)
}

calculate_gene_set_enrichment <- function(gep_vector, gene_set_db, top_n = 50) {
  gep_vector <- sort(gep_vector, decreasing = TRUE)
  top_genes <- names(gep_vector)[seq_len(min(top_n, length(gep_vector)))]
  total_genes <- length(gep_vector)

  results <- lapply(seq_along(gene_set_db$genesets), function(i) {
    genes_in_set <- intersect(gene_set_db$genesets[[i]], names(gep_vector))
    genes_found <- intersect(top_genes, genes_in_set)
    in_set_top <- length(genes_found)
    in_set_not_top <- length(genes_in_set) - in_set_top
    not_set_top <- length(top_genes) - in_set_top
    not_set_not_top <- total_genes - in_set_top - in_set_not_top - not_set_top

    fisher_out <- fisher.test(matrix(
      c(in_set_top, in_set_not_top, not_set_top, not_set_not_top),
      nrow = 2
    ))

    data.frame(
      Term = gene_set_db$geneset.names[[i]],
      FDR = NA_real_,
      `Genes Found` = in_set_top,
      `Gene Set Size` = length(genes_in_set),
      p_value = fisher_out$p.value,
      odds_ratio = unname(fisher_out$estimate),
      lower_ci = fisher_out$conf.int[[1]],
      genes = paste(genes_found, collapse = ", "),
      check.names = FALSE
    )
  })

  enrichment <- do.call(rbind, results)
  enrichment$FDR <- p.adjust(enrichment$p_value, method = "BH")
  enrichment[order(enrichment$lower_ci, decreasing = TRUE), ]
}

summarize_enrichments <- function(gep_vector, gene_set_databases, top_n = 50, n_terms = 5) {
  lapply(gene_set_databases, function(db) {
    enrichment <- calculate_gene_set_enrichment(gep_vector, db, top_n = top_n)
    enrichment <- enrichment[seq_len(min(n_terms, nrow(enrichment))), ]
    rownames(enrichment) <- NULL
    enrichment[, c("Term", "FDR", "Genes Found", "Gene Set Size")]
  })
}

get_top_cells <- function(split0_h, program_ids, top_cells = 100) {
  cells <- unlist(lapply(program_ids, function(program_id) {
    activity <- split0_h[, program_id]
    names(activity) <- rownames(split0_h)
    names(sort(activity, decreasing = TRUE))[seq_len(min(top_cells, length(activity)))]
  }))

  unique(cells)
}

summarize_single_r_scores <- function(singler_scores, cells) {
  colMeans(singler_scores[cells, , drop = FALSE], na.rm = TRUE)
}

summarize_disease_origin <- function(cells, metadata, sample_col = "Sample_ID", disease_col = "Disease") {
  samples <- sub(".*:", "", cells)
  disease <- metadata[[disease_col]][match(samples, metadata[[sample_col]])]
  prop.table(table(disease))
}

summarize_numbat_scores <- function(cells, numbat_metadata, min_cells = 10) {
  if (is.null(numbat_metadata)) {
    return(NULL)
  }

  numbat_cells <- sub("^(.*):(.*)$", "\\2_\\1", cells)
  matched <- numbat_metadata[numbat_metadata$cell.id %in% numbat_cells, ]
  matched <- matched[!is.na(matched$numbat.cell), ]

  if (nrow(matched) < min_cells) {
    return(NULL)
  }

  data.frame(
    SNV = c(mean(matched$numbat.cell, na.rm = TRUE), sd(matched$numbat.cell, na.rm = TRUE)),
    CNV = c(mean(matched$numbat.expression, na.rm = TRUE), sd(matched$numbat.expression, na.rm = TRUE)),
    Joint = c(mean(matched$numbat.joint, na.rm = TRUE), sd(matched$numbat.joint, na.rm = TRUE))
  )
}

score_label_matching_attempts <- function(predictions, truth, gep_col = "GEP", truth_col = "Paper") {
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))

  merged <- dplyr::left_join(truth, predictions, by = gep_col)
  attempt_cols <- grep("^Attempt", names(merged), value = TRUE)
  if (length(attempt_cols) == 0) {
    stop("No prediction columns found. Expected columns named Attempt*.")
  }

  tidyr::pivot_longer(
    merged,
    cols = dplyr::all_of(attempt_cols),
    names_to = "attempt",
    values_to = "prediction"
  ) |>
    dplyr::summarise(
      accuracy = mean(.data$prediction == .data[[truth_col]], na.rm = TRUE),
      .by = "attempt"
    )
}

write_gep_prompt_input <- function(gep_vector,
                                   gep_id,
                                   file,
                                   top_genes = 50,
                                   enrichment_list = NULL,
                                   top_cell_count = NULL,
                                   mean_scores = NULL,
                                   disease_freq = NULL,
                                   numbat_scores = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  top_gene_values <- sort(gep_vector, decreasing = TRUE)
  top_gene_values <- top_gene_values[seq_len(min(top_genes, length(top_gene_values)))]

  con <- file(file, open = "w")
  on.exit(close(con))

  writeLines(sprintf("GEP %d:\nTop %d genes with loadings:", gep_id, length(top_gene_values)), con)
  for (i in seq_along(top_gene_values)) {
    writeLines(sprintf("- %s: %s", names(top_gene_values)[i], round(top_gene_values[i], 4)), con)
  }

  if (!is.null(enrichment_list)) {
    writeLines("\nTop enrichment terms:", con)
    for (db_name in names(enrichment_list)) {
      writeLines(sprintf("\nFrom %s:", db_name), con)
      db <- enrichment_list[[db_name]]
      for (i in seq_len(nrow(db))) {
        writeLines(sprintf(
          "- %s | FDR: %.3g | Genes Found: %s/%s",
          gsub("_", " ", db$Term[i]),
          as.numeric(db$FDR[i]),
          db$`Genes Found`[i],
          db$`Gene Set Size`[i]
        ), con)
      }
    }
  }

  if (!is.null(top_cell_count)) {
    writeLines(sprintf("\nNumber of top cells activating the GEP: %s", top_cell_count), con)
  }

  if (!is.null(mean_scores)) {
    writeLines("\nMean singleR raw scores across top activated cells:", con)
    for (cell_type in names(sort(mean_scores, decreasing = TRUE))) {
      writeLines(sprintf("- %s: %s", cell_type, round(mean_scores[[cell_type]], 4)), con)
    }
  }

  if (!is.null(disease_freq)) {
    writeLines("\nDisease origin for top activated cells:", con)
    for (disease_type in names(sort(disease_freq, decreasing = TRUE))) {
      writeLines(sprintf("- %s: %s", disease_type, round(disease_freq[[disease_type]], 4)), con)
    }
  }

  if (!is.null(numbat_scores)) {
    writeLines("\nMean NumBat probability of malignancy across top activated cells:", con)
    for (type in colnames(numbat_scores)) {
      writeLines(sprintf(
        "- %s: %s (sd: %s)",
        type,
        round(numbat_scores[1, type], 3),
        round(numbat_scores[2, type], 3)
      ), con)
    }
  }
}
