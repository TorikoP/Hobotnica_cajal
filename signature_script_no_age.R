#!/usr/bin/env Rscript
library(dplyr)
source("/tank/projects/vpalagina_hobotnica/hobotnica/custom_function.R")
path_sig_matrix <- "/tank/projects/vpalagina_hobotnica/hobotnica/PhenoAgeV2_res_imputed_rand_signature/sig_matrix_no_age"
substract_age <- "/tank/projects/vpalagina_hobotnica/hobotnica/substracted_age_PhenoAgeV2_imputed_NA"
output_dir <- "/tank/projects/vpalagina_hobotnica/hobotnica/PhenoAgeV2_res_imputed_rand_signature/distrib_h_score_no_age"
existing_results_path <- "/tank/projects/vpalagina_hobotnica/hobotnica/PhenoAgeV2_res_imputed_rand_signature/PhenoAgeV2_H_scores_no_age.csv"


folder_list <- list.dirs(path_sig_matrix, full.names = TRUE, recursive = FALSE)
file_list_PhenoAge <- list.files(path = substract_age, full.names = TRUE)

result_table <- data.frame(
  Dataset_ID = character(), 
  H_score_no_age = numeric(), 
  p_value_no_age = numeric(), 
  stringsAsFactors = FALSE
)
if (file.exists(existing_results_path)) {
  existing_results <- read.csv(existing_results_path, stringsAsFactors = FALSE)
  processed_datasets <- existing_results$Dataset_ID
} else {
  existing_results <- data.frame(
    Dataset_ID = character(), 
    H_score_no_age = numeric(), 
    p_value_no_age = numeric(), 
    stringsAsFactors = FALSE
  )
  processed_datasets <- character()
}
for (i in 1:length(file_list_PhenoAge)) {
    file1 <- file_list_PhenoAge[i]
    sample_id1 <- basename(file1)

    if (sample_id1 %in% processed_datasets) {
        message(paste("Skipping", sample_id1, "- already processed."))
        next
    }
    data <- read.csv(file1, sep = ",", header = TRUE, row.names = 1)
    annotation <- data$Condition
    distMatrix <- kendall_dist(data)
    h_score <- Hobotnica(distMatrix, annotation)
    matching_folder <- folder_list[basename(folder_list) == tools::file_path_sans_ext(sample_id1)]

    if (length(matching_folder) > 0) {
        folder <- matching_folder[1]
        H_score_distrib <- ParallelHobotnica(folder, annotation)
        pval <- (1 + sum(H_score_distrib >= h_score)) / length(H_score_distrib)
        txt_output_path <- file.path(output_dir, paste0("H_score_distrib_", tools::file_path_sans_ext(sample_id1), ".txt"))
        write.table(H_score_distrib, txt_output_path, row.names = FALSE, col.names = FALSE)
        
        existing_results <- existing_results %>%
          add_row(Dataset_ID = sample_id1, H_score_no_age = h_score, p_value_no_age = pval)
    } else {
        warning(paste("No matching folder found for", sample_id1))
    }
    write.csv(existing_results, existing_results_path, row.names = FALSE)
}