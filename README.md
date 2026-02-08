# RNA_Seq-DATA-ANALYSIS
A comprehensive RNA-seq data analysis pipeline for transcriptome profiling. Normalization, differential gene expression analysis, and Data visualization are all part of the workflow. It uses patterns of gene expression to produce biologically significant discoveries.

**Mantle Cell Lymphoma (MCL)**
It is a rare and aggressive subtype of non-Hodgkin B-cell lymphoma that originates in the mantle zone of lymphoid follicles. It is characterized by the overexpression of CCND1 (Cyclin D1)

**OBJECTIVE**
     1. Raw RNA-seq data preprocessing and quality assessment
     2. Normalization of count data using DeSeq2
     3. Differential gene expression analysis
     4. Identification of upregulated and downregulated genes
     5. Visualization and biological interpretation

**File Description**

**Data_processing.R**
 Contains scripts for initial data preprocessing, including data import, cleaning, formatting, gene annotation, and preparation of count matrices for downstream RNA-seq analysis.


**Differential_gene_expression.R**
Implements differential gene expression analysis using the edgeR package. This script performs normalization, dispersion estimation, statistical testing, and identifies significantly upregulated and downregulated genes.


**Normalization_method.R**
Includes normalization procedures applied to RNA-seq data, such as TMM normalization, log transformation, and exploratory data analysis to assess sample distributions before and after normalization.


**README.md**
Provides an overview of the project, dataset description, analysis workflow, and interpretation of results.


**mart.103**
BioMart annotation file downloaded from Ensembl (release 103). It is used to map Ensembl gene IDs to HGNC gene symbols and to classify genes based on biotype (e.g., lncRNA, protein-coding).


**working_Data_df.csv**
Processed expression dataset after filtering, annotation, and normalization. This file is used as the main input for differential expression analysis.


**working_data.csv**
An intermediate processed dataset generated during the preprocessing steps before final filtering.


**upregulated_gene.csv**
List of significantly upregulated genes identified in Mantle Cell Lymphoma samples compared to controls.


**downregulated_gene.csv**
List of significantly downregulated genes identified in Mantle Cell Lymphoma samples compared to controls.


**volcano.png**
Volcano plot visualizing log2 fold change versus statistical significance of differentially expressed genes.


**norm.png**
Visualization of expression value distributions before and after normalization.
