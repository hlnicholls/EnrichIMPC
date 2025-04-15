# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(stats)
library(splitstackshape)
library(GeneOverlap)

#' Filter IMPC dataframe by phenotype search terms
#'
#' @param impc_df Dataframe of IMPC results
#' @param search_terms Character vector of terms to search within `Phenotype Hits`
#' @return Filtered dataframe with only relevant phenotypes
filter_phenotypes <- function(impc_df, search_terms) {
  search_pattern <- paste(search_terms, collapse = "|")
  specific_impc_df <- impc_df %>%
    filter(stringr::str_detect(`Phenotype Hits`, search_pattern))
  specific_impc_df$Gene <- toupper(specific_impc_df$Gene)
  return(specific_impc_df)
}

#' Filter the IMPC dataset using a user-provided gene list
#'
#' @param user_gene_list Dataframe with a column 'Gene'
#' @param impc_df_subset Filtered IMPC dataframe
#' @return Subset of rows where Gene matches user input
filter_by_user_genes <- function(user_gene_list, impc_df_subset) {
  filtered <- impc_df_subset[impc_df_subset$Gene %in% user_gene_list$Gene, ]
  select(filtered, Gene, `# Phenotype Hits`, `Phenotype Hits`)
}

#' Split phenotype annotations and retain shared phenotypes
#'
#' @param impc_df Full IMPC dataset
#' @param mouse_user_genes Subset of IMPC matching user genes
#' @return A list of split dataframes used in enrichment
split_phenotypes <- function(impc_df, mouse_user_genes) {
  impc_df_phenotypes <- cSplit(impc_df, "Phenotype Hits", sep = "::", direction = "long")
  impc_df_phenotypes_user_input <- cSplit(mouse_user_genes, "Phenotype Hits", sep = "::", direction = "long")
  colnames(impc_df_phenotypes_user_input)[2] <- 'gene_group_size'
  colnames(impc_df_phenotypes)[3] <- 'dataset_size'
  impc_df_phenotypes <- subset(impc_df_phenotypes, `Phenotype Hits` %in% impc_df_phenotypes_user_input$`Phenotype Hits`)
  list(impc_df_phenotypes = impc_df_phenotypes, impc_df_phenotypes_user_input = impc_df_phenotypes_user_input)
}

#' Perform enrichment tests per phenotype
#'
#' @param impc_df_phenotypes_user_input User gene-phenotype pairs
#' @param impc_df_phenotypes Background gene-phenotype pairs
#' @return Enrichment results with p-values and adjusted p-values
run_enrichment <- function(impc_df_phenotypes_user_input, impc_df_phenotypes) {
  d1_split <- split(impc_df_phenotypes_user_input, impc_df_phenotypes_user_input$`Phenotype Hits`)
  d2_split <- split(impc_df_phenotypes, impc_df_phenotypes$`Phenotype Hits`)
  stopifnot(all(names(d1_split) == names(d2_split)))

  enrichment_tests <- Map(function(d1, d2) {
    go.obj <- newGeneOverlap(d1$Gene, d2$Gene, genome.size = 20000)
    return(testGeneOverlap(go.obj))
  }, d1_split, d2_split)

  results <- tibble(Phenotype = names(enrichment_tests), enrichment_tests = enrichment_tests) %>%
    rowwise() %>%
    mutate(
      across(enrichment_tests,
             .fns = list(tested = getTested, pval = getPval),
             .names = '{.fn}')
    ) %>%
    select(-enrichment_tests)

  results$pval_adjusted <- p.adjust(results$pval, "fdr")
  results
}

#' Collapse enrichment results to remove duplicate genes per phenotype
#'
#' @param results Enrichment results
#' @param impc_df_phenotypes Background gene-phenotype pairs
#' @param user_gene_list User gene list
#' @return Deduplicated, merged results
collapse_results <- function(results, impc_df_phenotypes, user_gene_list) {
  df2 <- select(impc_df_phenotypes,  `MGI Gene Id`,  `Phenotype Hits`, Gene)
  df2 <- filter(df2, Gene %in% user_gene_list$Gene)
  setnames(df2, 'Phenotype Hits', 'Phenotype')
  results2 <- merge(results, df2, by='Phenotype', all.x=TRUE)
  results2 <- select(results2, -tested)

  collapsed_results <- results2 %>%
    group_by(Phenotype) %>%
    summarise(
      pval = first(pval),
      pval_adjusted = first(pval_adjusted),
      `MGI Gene Id` = toString(unique(`MGI Gene Id`)),
      Gene = toString(unique(Gene))
    ) %>%
    ungroup()
  collapsed_results
}

#' Write final enrichment results to CSV
#'
#' @param results Enrichment results
#' @param enrichment_path Output path (not used but retained for compatibility)
write_results <- function(results, enrichment_path) {
  fwrite(results, 'IMPC_results_table.csv')
}

#' Create bar plot for top 20 phenotypes with gene count and p-value labels
#'
#' @param mouse_user_genes User genes with phenotype hits
#' @param results Enrichment results
generate_plot <- function(mouse_user_genes, results) {
  user_mouse_count <- filter(mouse_user_genes, `# Phenotype Hits` > 0)
  user_mouse_count_split <- cSplit(user_mouse_count, "Phenotype Hits", sep = "::", direction = "long")
  user_mouse_count_split <- select(user_mouse_count_split, `Phenotype Hits`)
  colnames(user_mouse_count_split)[1] <- 'Phenotype'

  mouse_count_pheno <- aggregate(list(numdup=rep(1,nrow(user_mouse_count_split))), user_mouse_count_split, length)
  mouse_count_pheno <- mouse_count_pheno %>% dplyr::arrange(desc(numdup))
  mouse_count_pheno <- mouse_count_pheno[1:20,]
  colnames(mouse_count_pheno)[2] <- 'Count'
  mouse_count_pheno$Phenotype <- factor(mouse_count_pheno$Phenotype, levels = mouse_count_pheno$Phenotype)

  dt <- merge(mouse_count_pheno, results, by='Phenotype', all.x=TRUE)
  dt$pval_adjusted <- signif(dt$pval_adjusted, digits = 5)

  png("Genes_ranked_IMPC_enrichment.png", width = 14, height = 8, units = "in", res = 300)

  top_20_phenotypes_plot <- ggplot(data=dt, aes(x=Phenotype, y=Count)) +
    geom_bar(stat="identity", fill="darkslateblue") +
    theme(text = element_text(family = "Arial"),
          axis.text.x = element_text(color = "grey20", size = 16),
          axis.text.y = element_text(color = "grey20", size = 14),
          axis.title.x = element_text(color = "grey20", size = 16),
          axis.title.y = element_text(color = "grey20", size = 16)) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          legend.title=element_text(size=14), legend.text=element_text(size=14)) +
    coord_flip() +
    ggtitle("Top 20 Most Frequent IMPC Mouse Phenotypes")  +
    ylab("Gene Count") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text(aes(label=pval_adjusted), position=position_dodge(width=0.9), vjust=0.25, hjust = 1.1, colour = "white")

  print(top_20_phenotypes_plot)
  dev.off()
}

#' Run full enrichment pipeline on IMPC data
#'
#' @param user_gene_list Dataframe with user gene list (column 'Gene')
#' @param search_terms Terms to filter phenotypes
#' @param enrichment_path Output path for writing results (not currently used)
enrichIMPC <- function(user_gene_list, search_terms, enrichment_path) {
  mouse_user_genes <- filter_by_user_genes(user_gene_list, impc_df)
  impc_df_specific_terms <- filter_phenotypes(impc_df, search_terms)
  mouse_user_genes_specific_terms <- filter_by_user_genes(user_gene_list, impc_df_specific_terms)

  splits <- split_phenotypes(impc_df, mouse_user_genes)
  results <- run_enrichment(splits$impc_df_phenotypes_user_input, splits$impc_df_phenotypes)
  collapsed <- collapse_results(results, splits$impc_df_phenotypes, user_gene_list)
  write_results(collapsed, enrichment_path)
  generate_plot(mouse_user_genes, results)
}

