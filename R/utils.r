#' @title Olink Dist Plot
#' @description Create a Distance plot
#' @param npx NPX dataset (output from OlinkAnalyze::read_NPX)
#' @param panel Panel Name
#' @param filename filename, the suffix will determine the filetype
#' @param width Figure width,  Default: 7
#' @param height Figure height,  Default: 7
#' @return distance plot
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  olink_dist(plasma_npx, "Olink CARDIOVASCULAR III", filename="test.png")
#'  }
#' }
#' @seealso
#'  \code{\link[OlinkAnalyze]{olink_dist_plot}}
#' @rdname olink_dist
#' @export
#' @import OlinkAnalyze
#' @import ggplot2
#' @import dplyr
olink_dist <- function(npx, panel, filename, width = 7, height = 7) {
  olink_dist_plot(npx %>% dplyr::filter(Panel == panel)) +
  theme(text = element_text(size=6),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
  ggsave(filename, width = width, height = height)
}

#' @title gsea_olink_run
#' @description Run GSEA analysis from the output of an Olink statistical test using the estimate as the ranking metric
#' @param x output from an Olink statistical test
#' @param gs Gene set collection, Default: gsea_sets
#' @param minGSSize minimum gene set size to be considered, Default: 2
#' @param contrastName Name of contrast for filename purposes 
#' @param save save results, logical, Default: TRUE
#' @param saveFolder path to folder where results are to be save if save=TRUE
#' @return returns a cluterProfiler result object, suitable for further visualization
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[tibble]{enframe}}
#'  \code{\link[openxlsx]{write.xlsx}}
#' @rdname gsea_olink_run
#' @export 
#' @importFrom tibble deframe
#' @import clusterProfiler
#' @importFrom openxlsx write.xlsx
gsea_olink_run <- function(x, gs=gsea_sets, minGSSize = 2, contrastName = compname, save=TRUE, saveFolder) {
  # Create genelist
  geneList <- x %>% 
    arrange(-estimate) %>% 
    select(Assay, estimate) %>% 
    tibble::deframe()
  # PErform enrichment
  gsea_out <- clusterProfiler::GSEA(geneList,
             TERM2GENE=gs,
             exponent = 1,
             minGSSize = minGSSize,
             maxGSSize = 500,
             eps = 1e-10,
             pvalueCutoff = 1,
             pAdjustMethod = "BH",
             verbose = FALSE,
             seed = 120,
             by = "fgsea")
  
  if(save==TRUE){
    openxlsx::write.xlsx(gsea_out, paste0(saveFolder, contrastName, "_GSEA_Results_Full.xlsx"))
  }
  
  gsea_out
}


#' @title gsea_olink_viz
#' @description Create several visualizations starting from a cluterProfiler result object
#' @param gsea_res clusterProfiler result object
#' @param saveFolder path to folder where results are to be saved
#' @param nterms Number of terms to be visuzlized for several visualizations, Default: 10
#' @param contrastName Name of contrast for filename purposes , Default: compname
#' @return Several visualizations are saved to disk in the folder provided
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[enrichplot]{c("dotplot,compareClusterResult-method", "dotplot", "dotplot")}},\code{\link[enrichplot]{c("pairwise_termsim", "pairwise_termsim")}},\code{\link[enrichplot]{c("emapplot", "emapplot")}},\code{\link[enrichplot]{c("gseaplot2", "gseaplot2")}}
#'  \code{\link[grid]{arrow}}
#' @rdname gsea_olink_viz
#' @export 
#' @importFrom enrichplot dotplot pairwise_termsim emapplot gseaplot2
#' @importFrom grid arrow
gsea_olink_viz <- function(gsea_res, saveFolder, nterms = 10, contrastName = compname){

  viz_obj <- gsea_res %>% 
    as.data.frame() %>% 
    arrange(-NES) %>% 
    {rbind(head(.,nterms), tail(.,nterms))} %>%
    mutate(col = case_when(
      p.adjust < 0.05 ~ "black",
      TRUE ~ "grey"),
      hjust = case_when(
        NES > 1 ~ 1,
        TRUE ~ 0))
  
  # barplot
    bar_viz <- ggplot(viz_obj, aes(reorder(ID, NES), NES, label=ID, fill=col, hjust=hjust)) +
    geom_bar(stat = "identity") +
    scale_fill_identity() +
    geom_text(size = 3, aes(y = 0)) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both"))) +
    ylab(paste(paste0("Downregulated in ",strsplit(compname, split = "_vs_")[[1]][1]), "NES", paste0("Upregulated in ", strsplit(compname, split = "_vs_")[[1]][1]), sep = "              "))
  ggsave(paste0(saveFolder, contrastName, "_barplot.pdf"))
  
  # Enrichment map - needs to be subsetted with asis=T to keep structure
  red_obj <- gsea_res[gsea_res$ID %in% viz_obj$ID, asis=T]
  emap <- enrichplot::pairwise_termsim(red_obj, method = "JC", semData = NULL, showCategory = nrow(red_obj))
  enrichplot::emapplot(emap, showCategory=nrow(red_obj), layout = "kk", color = "NES")
  ggsave(paste0(saveFolder, contrastName, "_emap.pdf"))
  
  # GSEA plots
  dir_gsea <- paste0(saveFolder, "GSEA_PLOTS/")
  dir.create(dir_gsea, showWarnings = FALSE, recursive = TRUE)
  for(i in 1:nrow(viz_obj)){
    enrichplot::gseaplot2(red_obj, geneSetID = i, title = red_obj$Description[i])
    ggsave(paste0(dir_gsea, contrastName, "_gseaplot_", gsub(" |/", "_", substr(red_obj$Description[i], 0, 15)),  ".pdf"))
}
}


