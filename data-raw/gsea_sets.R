## code to prepare `gsea_sets` dataset goes here
library("clusterProfiler")
library("dplyr")

gsea_sets <- clusterProfiler::read.gmt("inst/extdata/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_symbol.gmt")
gsea_sets$term <- gsub("%.*$", "", gsea_sets$term) 

usethis::use_data(gsea_sets, overwrite = TRUE)