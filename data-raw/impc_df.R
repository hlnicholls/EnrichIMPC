library(readr)

impc_df <- read_csv("data-raw/phenotypeHitsPerGene.csv.gz")
impc_df$Gene <- toupper(impc_df$`Gene Symbol`)
usethis::use_data(impc_df, compress = "xz", overwrite = TRUE)
