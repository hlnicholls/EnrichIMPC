# EnrichIMPC
EnrichIMPC is an R package for performing phenotype enrichment analysis using the International Mouse Phenotyping Consortium (IMPC) dataset. It allows users to identify overrepresented mouse phenotypes associated with their input gene list, filtered by user-defined phenotype keywords.

## Overview
This tool is intended to support gene set interpretation and sensitivity analyses by:
 - Searching IMPC phenotypes using user-provided search terms.
 - Identifying matches between the user gene list and phenotype hits.
 - Performing enrichment tests (based on GeneOverlap) to compare phenotype associations in the user genes vs. the full IMPC dataset.
 - Adjusting p-values and providing summary tables and publication-quality bar plots.

## Installation
Clone the package and use devtools::load_all() or install from your local directory:

```
# Clone and install
devtools::load_all("/path/to/EnrichIMPC")
```

### Dependencies
```
install.packages(c("data.table", "dplyr", "ggplot2", "ggpubr", "magrittr", "splitstackshape", "stringr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GeneOverlap")
```

## Example Usage
```
# Load your package
devtools::load_all("/path/to/EnrichIMPC")

# Load or prepare your inputs:
# - user_gene_list: a dataframe with a column named 'Gene' (uppercase preferred)
# - search_terms: a character vector of phenotypes of interest

# Example:
user_gene_list <- read.csv("user_genes.csv")  # must include 'Gene' column
search_terms <- c("cardiovascular", "heart", "blood pressure")

# Run the enrichment analysis
enrichIMPC(user_gene_list, search_terms, enrichment_path = "output/")
```

## Output
- ```IMPC_results_table.csv```: tabular enrichment results including p-values and FDR-adjusted values.
- ```Genes_ranked_IMPC_enrichment.png```: bar plot of top 20 phenotypes by gene count and enrichment.

