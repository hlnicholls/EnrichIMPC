
# run_example.R

devtools::load_all("/Users/hannahnicholls/GitHub/EnrichIMPC")  # change to your actual path

data("impc_df")

# User-provided gene list
user_gene_list <- data.frame(Gene = c("ACE", "MYH7", "TTN", "LMNA", "BAG3", "TNNT2", 'FLNC', 'FHOD3', 'MYH6'))

# Phenotype search terms
search_terms <- c('cardio', 'heart', 'cardiac', 'muscle', 'adrenal', 'obesity')

# Run full enrichment workflow
enrichIMPC(user_gene_list, search_terms, enrichment_path = "IMPC_results_table.csv")
