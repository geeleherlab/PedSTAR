# R Scripts

Scripts for GEP prompt generation and summary analyses.

- `00_utils_gep_annotation.R`: shared helper functions for GEP summaries and prompt-input files.
- `01_generate_gep_annotation_inputs.R`: writes one PedSTAR annotation input text file per GEP.
- `02_generate_gosh_benchmark_inputs.R`: writes benchmark input text files for the GOSH GEPs.
- `03_score_gosh_label_matching.R`: scores AI label-matching attempts against human labels.
- `04_calculate_gep_coexpression.R`: calculates pairwise Spearman co-expression between GEP activity profiles.
- `05_summarize_gtex_expression.R`: summarizes GTEx tissue expression for selected candidate genes.
- `06_summarize_depmap_expression.R`: summarizes DepMap RNA and protein data for selected candidate genes.
