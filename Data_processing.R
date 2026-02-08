library(DESeq2)
library(dplyr)
library(tidyverse)
library(Biostrings)

# data loading in the R environment
read.csv("~/Documents/IIT-Work/lnc_RNAseq/working_data.csv", header = T) -> exp_data

# Download manually from GEO, or programmatically:
url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL19nnn/GPL19023/annot/GPL19023.annot.gz"

# Download file
download.file(url, destfile = "GPL19023.annot.txt", mode = "wb")

# Download manually from GEO, or programmatically:
url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL19nnn/GPL19023/annot/GPL19023.annot.gz"

# Download file
download.file(url, destfile = "GPL19023.annot.txt", mode = "wb")

# Read into R
anno <- read.delim("GPL19023.annot.gz", stringsAsFactors = FALSE)

# Inspect columns
head(colnames(anno))

# Read annotation file
anno <- read.delim("GPL_19023.txt", sep = "\t", stringsAsFactors = FALSE)

# Make DNAStringSet
probes <- DNAStringSet(anno$SEQUENCE)
names(probes) <- anno$ID

##################### Download hg19 genome (FASTA) if not already###############
#BASH Command
#wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#gunzip hg19.fa.gz

#Run BLAT
#Install Blat
#blat hg19.fa GPL19023_probes.fa probes_hg19.psl
################################################################################


library(GenomicRanges)
library(rtracklayer)

# Import probe alignment results (after converting .psl → .bed)
#(converting .psl → .bed) #BASH
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslToBed
#chmod +x pslToBed
#sudo mv pslToBed /usr/local/bin/

#pslToBed input.psl output.bed
probes_gr <- import("probe_hg19.bed")
#linux terminal command
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
#gunzip gencode.v19.annotation.gtf.gz

# Import GTF annotation (Ensembl/Gencode for hg19)
gtf <- import("gencode.v19.annotation.gtf")

# Keep only "gene" features
genes_gr <- gtf[gtf$type == "gene"]

# Find overlaps
hits <- findOverlaps(probes_gr, genes_gr)

# Map probes to genes
# Build mapping table
probe_to_gene <- data.frame(
  ProbeID    = mcols(probes_gr)$name[queryHits(hits)],   # probe IDs from your alignment
  GeneID     = mcols(genes_gr)$gene_id[subjectHits(hits)],   # Ensembl gene ID
  GeneSymbol = mcols(genes_gr)$gene_name[subjectHits(hits)]  # HGNC gene symbol
)

# remove duplicate and Keep only the first occurrence of each ProbeID
probe_to_gene_unique <- probe_to_gene[!duplicated(probe_to_gene$ProbeID), ]

#1. Find common ProbeIDs (intersection)
common_probes <- intersect(exp_data$ProbeName, probe_to_gene$ProbeID)

#######check duplicates#######################
exp_data$ProbeName[duplicated(exp_data$ProbeName)] %>% length()

# remove duplicates and keep only the first occurrence of each ProbeName
exp_data_unique <- exp_data[!duplicated(exp_data$ProbeName), ]

library(dplyr)

# Make probe_to_gene unique
probe_to_gene_unique <- probe_to_gene %>%
  distinct(ProbeID, .keep_all = TRUE)

# Make exp_data unique
exp_data_unique <- exp_data %>%
  distinct(ProbeName, .keep_all = TRUE)

# Now join
output_data <- exp_data_unique %>%
  left_join(probe_to_gene_unique, by = c("ProbeName" = "ProbeID"))

#  Remove version numbers from GeneID
valid_output_data <- output_data %>%
  mutate(GeneID = sub("\\..*$", "", GeneID))

############################################################################################
############################################################################################

library(biomaRt)
library(httr)
library(jsonlite)

#load the file properly
mart_hgnc_v2 <- read.delim(
  "~/Documents/IIT-Work/lnc_RNAseq/mart.103.v2",
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE
)

# Function to get HGNC symbols via Ensembl REST API with batching and retries
get_hgnc_batch_retry <- function(ensembl_ids, batch_size = 100, max_retries = 3, delay = 1) {
  
  ensembl_ids <- clean_ensembl_ids(ensembl_ids)  # remove versions
  chunks <- split(ensembl_ids, ceiling(seq_along(ensembl_ids)/batch_size))
  all_results <- list()
  
  for (i in seq_along(chunks)) {
    ids <- chunks[[i]]
    success <- FALSE
    attempt <- 1
    
    while(!success && attempt <= max_retries) {
      Sys.sleep(delay)  # short wait
      r <- POST(
        url = "https://rest.ensembl.org/xrefs/id",
        content_type_json(),
        accept_json(),
        body = toJSON(list(ids = ids), auto_unbox = TRUE)
      )
      
      if (r$status_code == 200) {
        tmp <- fromJSON(content(r, "text", encoding = "UTF-8"))
        res <- do.call(rbind, lapply(tmp, function(x) {
          if(length(x) > 0) {
            data.frame(ensembl_gene_id = x$id, hgnc_symbol = x$display_id)
          } else {
            data.frame(ensembl_gene_id = x$id, hgnc_symbol = NA)  # keep ID if no match
          }
        }))
        all_results[[i]] <- res
        success <- TRUE
      } else{
        message("Batch ", i, " failed on attempt ", attempt, ". Retrying...")
        attempt <- attempt + 1
        Sys.sleep(delay * attempt)  # exponential backoff
      }
    }
    
    if (!success) {
      message("Batch ", i, " failed after ", max_retries, " attempts. IDs: ", paste(ids, collapse = ", "))
      all_results[[i]] <- data.frame(ensembl_gene_id = ids, hgnc_symbol = NA)
    }
  }
  
  do.call(rbind, all_results)
}

# Clean Ensembl IDs: remove version numbers
clean_ensembl_ids <- function(ids) {
  sub("\\..*$", "", ids)
}

# Example usage:
valid_ids <- EnsembleIds$EnsembleIds  # your vector of Ensembl IDs
hgnc.list <- get_hgnc_batch_retry(valid_ids)

#merge output_data and mart_hgnc_v2 
final_data <- valid_output_data %>%
  right_join(mart_hgnc_v2, by = "GeneID")
#remove rows containing Na values
final_data_clean <- final_data %>%
  drop_na()


#remove unnecessary columns
final_data_clean <- final_data %>%
  select(-c(
    `HGNC ID`,
    `Status`,
    `Previous symbols`,
    `Alias symbols`,
    `Chromosome`,
    `Accession numbers`,
    `RefSeq IDs`,
    `Previous name`,
    `Alias names`
  ))

final_data_clean <- final_data_clean %>%
  rename(
    Old_GeneSymbol = `GeneSymbol.x`,
    New_GeneSymbol = `GeneSymbol.y`,
    HGNC_Symbol    = `Approved symbol`,
    p1Pve          = `X.149PE...raw.`,
    p1Nve          = `X.149PE...raw..1`,
    p2Pve          = `X.225PE...raw.`,
    p2Nve          = `X.225PE...raw..1`,
    p3Pve          = `X.320PE...raw.`,
    p3Nve          = `X.320PE...raw..1`,
    control1       = `X.Control_1..raw.`,
    control2       = `X.Control_2..raw.`
  )

#FINAL CHECk
sum(duplicated(final_data_clean$ProbeName))

final_data_clean_nodup <- final_data_clean %>%
  distinct(HGNC_Symbol, .keep_all = TRUE)

sum(is.na(final_data_clean_nodup$`HGNC_Symbol`))

sum(final_data_clean_nodup$`HGNC_Symbol` == "", na.rm = TRUE)

Test_data <- final_data_clean_nodup %>%
  filter(if_all(everything(), ~ . != ""))



























