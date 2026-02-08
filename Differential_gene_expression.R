working_Data_df = as.data.frame(working_Data)
working_Data_df <- working_Data_df %>%
  select(
    HGNC_Symbol, p1Pve, p1Nve ,p2Pve , p2Nve , p3Pve ,p3Nve , control1 ,control2,            # keep this as first column
    -c(`Old_GeneSymbol`, `New_GeneSymbol`, `Ensembl_gene_ID`, `ProbeName`)  # columns to remove
  )

# Convert to matrix for normalization
working_Data_matrix <- as.matrix(working_Data[, 2:9])
rownames(working_Data_matrix) <- working_Data_df$HGNC_Symbol

#Create Sample Info

sample_info <- data.frame(
  row.names = c("p1Pve", "p1Nve" ,"p2Pve" , "p2Nve" , "p3Pve" ,"p3Nve" , "control1" ,"control2"),
  Condition = c("PVE", "NVE", "PVE", "NVE", "PVE", "NVE", "Control", "Control")
)

# Assume first column = gene IDs, rest = counts
HGNC_Symbol <- working_Data_df[[1]]

# Convert all other columns to numeric except the 1st
counts_only <- data.frame(lapply(working_Data_df[ , -1], function(x) as.numeric(as.character(x))))
# Put gene IDs back
clean_numeric <- cbind(HGNC_Symbol = HGNC_Symbol, counts_only)


# Convert to matrix
working_data_matrix_int <- as.matrix(counts_only)

# 3. Force storage mode to integer (important for DESeq2/edgeR)
mode(working_data_matrix_int) <- "integer"

# 4. Add gene IDs as rownames
rownames(working_data_matrix_int) <- HGNC_Symbol

library(DESeq2)
                                 
DGE_Data <- DESeqDataSetFromMatrix(countData = working_data_matrix_int,
                                   colData = sample_info,
                                   design = ~ Condition)

#Run DESeq2 Analysis
DGE_Data <- DESeq(DGE_Data) 
#Create the results object
result <- results(DGE_Data)

result_clean <- result[!is.na(result$padj), ]
sum(is.na(result$`padj`))

# Convert to dataframe for easy handling
result_clean <- as.data.frame(result)

# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1  # means 2-fold change

# Add regulation column
result_clean$regulation <- "NotSig"
result_clean$regulation[result_clean$padj < padj_cutoff & result_clean$log2FoldChange > log2fc_cutoff] <- "Up"
result_clean$regulation[result_clean$padj < padj_cutoff & result_clean$log2FoldChange < -log2fc_cutoff] <- "Down"


# Upregulated genes
up_genes <- subset(result_clean, padj < padj_cutoff & log2FoldChange > log2fc_cutoff)
up_genes <- subset(result_clean, regulation == "Up")
# Downregulated genes
down_genes <- subset(result_clean, padj < padj_cutoff & log2FoldChange < -log2fc_cutoff)
#Extract all rows where the regulation column is "down."
down_genes <- subset(result_clean, regulation == "Down")

combined_df_up <- rbind(up_genes, down_genes)


library(tidyverse)

# VOLCANO Plots
ggplot(result_clean, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value") +
  theme_minimal()
