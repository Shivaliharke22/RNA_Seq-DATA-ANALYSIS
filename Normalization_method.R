############## NORMALIZATION METHODS ####################
library(edgeR)
library(tidyverse)
library(rlang)
library(biomaRt)
library(janitor)
library(GEOquery)
library(limma)
library(preprocessCore)
### 1. Quantile normalization Normalize a matrix of Agilent expression values (already summarized)
# Install and load the preprocessCore package

# Ensure numeric
patient_data_matrix <- apply(patient_data_matrix, 2, as.numeric)

# remove the first column (GeneSymbol)
patient_data_matrix <- as.matrix(patient_data_clean[, -1])

# Perform Quantile Normalization
patient_data_normalized_matrix <- normalize.quantiles(patient_data_matrix)
# Convert back to dataframe
patient_data_normalized_df <- as.data.frame(patient_data_normalized_matrix)
# Restore column and row names
colnames(patient_data_normalized_df) <- colnames(patient_data_clean)[-1]


#### ZSCORE NORMALIZATION ###########


# Apply Z-score normalization per lncRNA (row)
z_score_row <- function(row) {
  s <- sd(row)
  if (s == 0) return(rep(0, length(row)))  # or NA
  (row - mean(row)) / s
}

# Apply Z-score normalization per row (lncRNA)
z_normalized_matrix <- t(apply(filtered_counts[,2:9], 1, z_score_row))
min(z_normalized_matrix)
max(z_normalized_matrix)
#convert matrix into Dataframe 

z_normalized_matrix <- cbind(GeneSymbol, z_normalized_matrix)
#z_normalized_matrix <- as.matrix(z_normalized_df)


# Convert all elements of the matrix to numeric
#z_normalized_matrix[,2:9] <- apply(z_normalized_matrix[,2:9], 2, as.numeric)


#######log2 Normalization ##########

# If patient_data_matrix is a matrix:
log2_normalized_matrix <- log2(filtered_counts[,2:9] + 1)

# add column GeneSymbol in normalized matrix
log2_normalized_matrix <- cbind(GeneSymbol, log2_normalized_matrix)
dim(log2_normalized_matrix)
