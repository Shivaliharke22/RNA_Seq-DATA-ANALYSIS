# RNA-seq Analysis of Mantle Cell Lymphoma

## Repository Overview
This repository contains an RNA-seq data analysis pipeline
to study gene expression changes in Mantle Cell Lymphoma.

---

## File Description

| File | Description |
|-----|-------------|
| Data_processing.R | Data preprocessing and annotation |
| Differential_gene_expression.R | edgeR-based DE analysis |
| Normalization_method.R | Normalization techniques |
| mart.103 | Ensembl BioMart annotation file |
| upregulated_gene.csv | Upregulated genes |
| downregulated_gene.csv | Downregulated genes |

---

## Analysis Workflow
1. Data preprocessing
2. Normalization
3. Differential expression
4. Visualization



**Mantle Cell Lymphoma (MCL)**
It is a rare and aggressive subtype of non-Hodgkin B-cell lymphoma that originates in the mantle zone of lymphoid follicles. It is characterized by the overexpression of CCND1 (Cyclin D1)

**OBJECTIVE**
     
	 1. Raw RNA-seq data preprocessing and quality assessment
	
     2. Normalization of count data using DeSeq2
	
     3. Differential gene expression analysis
	
     4. Identification of upregulated and downregulated genes
	
     5. Visualization and biological interpretation



**upregulated_gene.csv**: 

    List of significantly upregulated genes identified in Mantle Cell Lymphoma samples compared to controls.

**downregulated_gene.csv**:

    List of significantly downregulated genes identified in Mantle Cell Lymphoma samples compared to controls.

**volcano.png**:

    Volcano plot visualizing log2 fold change versus statistical significance of differentially expressed genes.

**norm.png**:

    Visualization of expression value distributions before and after normalization.
