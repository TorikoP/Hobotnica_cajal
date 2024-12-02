RandomSignature <- function(dataset, n_sites, n_permutation) {
  
  annotation <- dataset %>% select(Condition)
  data <- dataset %>% select(-Condition, -Age)

  h_scores <- vector("list", n_permutation)
  column_names_list <- vector("list", n_permutation)

  for (i in 1:n_permutation) {
    random_columns <- sample(1:ncol(data), n_sites, replace = FALSE)
    subset_matrix <- data[, random_columns, drop = FALSE]
    distMatrix <- kendall_dist(subset_matrix)
    h_score <- Hobotnica(distMatrix, annotation)
    
    h_scores[[i]] <- h_score
    column_names_list[[i]] <- colnames(subset_matrix)
  }
  return(list(h_scores = h_scores, column_names = column_names_list))
}

##graph for randome signature visualization 
plot_signature_distrib <- function(distribution_of_h_scores, h_score, pval) {
random_h_scores <- as.numeric(distribution_of_h_scores)
real_h_score <- as.numeric(round(h_score, 3))
pval <- pval

breaks <- seq(0, 1, by = 0.001)

plot <- ggplot(data = data.frame(h_scores = random_h_scores), aes(x = h_scores)) +
  geom_histogram(binwidth = 0.002, fill = "skyblue", color = "black", alpha = 0.9) +
  
  # Highlight the real h_score
  geom_vline(aes(xintercept = real_h_score), color = "red", linetype = "dashed", size = 0.7) +
  
  # Add labels and title
  labs(title = "Distribution of Random h_scores with Real h_score Highlighted",
       x = "Distribution of randome H-scores",
       y = "Frequency") +
  
  # Annotate the real h_score
  annotate("text", x = real_h_score - 0.1, 
           y = max(hist(random_h_scores, breaks = breaks, plot = FALSE)$counts),
           label = paste("H-score=", real_h_score, "\np-value =", pval), color = "red", vjust = + 1.0)
           
return(plot)
}

kendall_dist <- function(ds_samples){
    if ("Condition" %in% colnames(ds_samples) || "Age" %in% colnames(ds_samples)) {
    # Remove the columns 'Condition' and 'Age' if they exist
    matrix <- ds_samples %>% select(-Condition, -Age)
  } else {
        matrix <- ds_samples
    }
    distMatrix <- Dist(matrix, method = "kendall", nbproc = 12) #32
        # Print the resulting distance matrix
return(distMatrix)
}

PCA_my_plot <- function(data, filename, phrase){

  # Select columns for PCA, excluding 'Condition' and 'Age'
  pca_data <- data %>% select(-Condition, -Age)

  # Convert all columns to numeric
  pca_data <- pca_data %>% mutate_all(as.numeric)

  # Perform PCA
  pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

  # Plot the PCA result with percentage of removed columns as subtitle
   p <- autoplot(pca_result, data = data, colour = 'Condition') +
    labs(title = filename, subtitle = paste0(phrase)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  return(p)
}

Hobotnica <- function(distMatrix, annotation){
    if (typeof(annotation) == "list") {
        annotation <- as.vector(unlist(annotation))
    } else {
        annotation <- as.vector(annotation)
    }
    rank.m <- as.matrix(distMatrix) # transform distance matrix to matrix object
    rank.m[lower.tri(rank.m)] <- rank(rank.m[lower.tri(rank.m)]) # transform distances to ranks
    rank.m[upper.tri(rank.m)] <- rank(rank.m[upper.tri(rank.m)]) #

    inclass_sum <- 0
    classes <- unique(annotation) # unique classes
    Ns <- vector()

    for (i  in 1:length(classes)){

        clas <- classes[i]
        class_samples <- which(annotation == clas)
        l_tmp <- length(class_samples)
        Ns[i] <- l_tmp
        tmp_sum_inclass <- sum(rank.m[class_samples,class_samples]) # sum of ranks, describing in-class distances
        inclass_sum <- inclass_sum + tmp_sum_inclass


    }
    Ns_sum <- sum(Ns)
    biggest_bossible_rank <-  Ns_sum * (Ns_sum - 1)/2
    number_of_unique_inclass_elements <-  sum(Ns * (Ns-1))/2
    maximal_value <- number_of_unique_inclass_elements * (2*biggest_bossible_rank - number_of_unique_inclass_elements + 1)
    minimal_value <- number_of_unique_inclass_elements* (1 + number_of_unique_inclass_elements)

    normalization_factor <- maximal_value - minimal_value
    return (max(0, 1 - (inclass_sum - minimal_value)/normalization_factor ))

}

GenerateRandomSignature <- function(dataset_path, n_sites, n_permutations, sample_id) {
  
  library(reticulate) # for python 
  use_python("/home/vpalagina/miniconda3/envs/my_env/bin/python", required = TRUE)
  source_python("read_pickle.py")

  dataset <- read_pickle_file(dataset_path)
  annotation <- dataset %>% select(Condition)
  data <- dataset %>% select(-Condition, -Age)
  ann_and_age <- dataset %>% select(Condition, Age)
  
  for (i in 1:n_permutations) {
    random_columns <- sample(colnames(data), n_sites, replace = FALSE)
    subset_matrix <- data[, random_columns, drop = FALSE]
    sig_matrix <- cbind(subset_matrix, ann_and_age)
    column_names_file <- sprintf("/tank/projects/vpalagina_hobotnica/hobotnica/PhenoAgeV2_res_imputed_rand_signature/signature/%s/signature_%03d.txt", sample_id, i)
    subset_matrix_file <- sprintf("/tank/projects/vpalagina_hobotnica/hobotnica/PhenoAgeV2_res_imputed_rand_signature/sig_matrix/%s/subset_matrix_%03d.csv", sample_id, i)
    
    # Ensure the directories exist
    dir.create(dirname(column_names_file), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(subset_matrix_file), recursive = TRUE, showWarnings = FALSE)
    
    # Write random column names and subset matrix
    writeLines(random_columns, con = column_names_file)
    write.csv(sig_matrix, file = subset_matrix_file, row.names = TRUE)
    }
  }

  library(foreach)
library(doParallel)
library(amap)

ParallelHobotnica <- function(folder_csv_matrices, annotation){

n_cores <- 12                 #sum = 48 cores detectCores() - 40  
cl <- makeCluster(n_cores)
registerDoParallel(cl)

csv_files <- list.files(folder_csv_matrices, full.names = TRUE)
H_scores <- foreach(i = 1:length(csv_files), .packages = c("amap", "dplyr"),
                    .export = c("kendall_dist", "Hobotnica"),
                    .combine = c) %dopar% {
  
  subset_matrix <- read.csv(csv_files[i], row.names = 1)
  DistMatrix <- kendall_dist(subset_matrix)
  H_result <- Hobotnica(DistMatrix, annotation)
  return(H_result)
}
stopCluster(cl)
return(H_scores)

}