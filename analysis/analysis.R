# Clear global env from earlier session and remove result folder for a clean run when sourcing script - comment out if not necessary
rm(list=ls())
unlink("Results", recursive=TRUE)

library("Analysis5204")
data(list=c("plasma_metadata", "plasma_npx", "serum_metadata", "serum_npx", "gsea_sets", "lum_dat"), package = "Analysis5204")
library("ggplot2")
library("dplyr")
library("SummarizedExperiment")
library("pcaExplorer")
library("emmeans")
library("ComplexHeatmap")

#################################################################################
# Raw data files are loaded automatically when loading the Analysis5204 package.#
#################################################################################

# First analyze plasma only

# TODO #2 First issue attempt


# Three proteins (MCP-1, OPG and uPA) are measured in both the Cardiovascular and the Inflammation panel. As seen in the values correlate well between the panels so only the Inflammation data is kept for these proteins. Has to been done before making the summarizedExperiment because pivot_wider otherwise fails.
doubles <- c("MCP-1", "OPG", "uPA")
plasma_npx <- plasma_npx[!(plasma_npx$Assay %in% doubles & plasma_npx$Panel == "Olink CARDIOVASCULAR III"),]

# Add 4 additional luminex proteins to plasma_npx; they are on a log2 scale (like the NPX data) - similar result of running anova seprately or including wiht other samples
lum_dat$`CCL18/PARC` <- as.numeric(lum_dat$`CCL18/PARC`)
lum_dat$`CA15-3/MUC-1` <- as.numeric(lum_dat$`CA15-3/MUC-1`)
lum_dat$`TIMP-1` <- as.numeric(lum_dat$`TIMP-1`)
lum_dat$C5a <- as.numeric(lum_dat$C5a)
lum_clean <- lum_dat %>% 
  mutate(`CCL18/PARC` = log2(`CCL18/PARC`),
         `CA15-3/MUC-1` = log2(`CA15-3/MUC-1`),
         `TIMP-1` = log2(`TIMP-1`),
         C5a = log2(C5a)) %>%
  mutate(`CCL18_PARC` = `CCL18/PARC`,
         `CA15-3_MUC-1` = `CA15-3/MUC-1`) %>%
  tidyr::pivot_longer(cols = c(`CCL18_PARC`, `CA15-3_MUC-1`, `TIMP-1`, C5a), names_to = "Assay", values_to = "NPX") %>%
  mutate(Panel = "Olink CARDIOVASCULAR III",
         MissingFreq = 0,
         Panel_Version = "Lum",
         PlateID = 1,
         QC_Warning = "Pass",
         LOD = 0,
         Normalization = "Log Normalized",
         UniProt = case_when(
           Assay == "TIMP-1" ~ "P01033",
           Assay == "CCL18_PARC" ~ "P55774",
           Assay == "CA15-3_MUC-1" ~ "P15941",
           Assay == "C5a" ~ "P01031"),
         OlinkID = case_when(
           Assay == "TIMP-1" ~ "OID50001",
           Assay == "CCL18_PARC" ~ "OID50002",
           Assay == "CA15-3_MUC-1" ~ "OID50003",
           Assay == "C5a" ~ "OID50004")) %>%
  rename(SampleID = Sample.ID) %>%
  left_join(plasma_npx %>% dplyr::select(SampleID, Index) %>% distinct()) %>%
  dplyr::select(SampleID, Index, OlinkID,UniProt, Assay,MissingFreq , Panel,Panel_Version, PlateID,QC_Warning, LOD,  NPX, Normalization)

plasma_npx <- rbind(plasma_npx, lum_clean)

# Make a SummarizedExperiment object
npx <- plasma_npx %>% 
  dplyr::select(SampleID, Assay, NPX) %>% 
  tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
  tibble::column_to_rownames("SampleID") %>%
  t()
metadata <- plasma_metadata %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("id") %>%
  mutate(diagnosis = recode_factor(diagnosis, `0` = "MPA", `1` = "GPA")) %>%
  mutate(abs = dplyr::case_when(
    `pr3-anca` == 1 ~ "PR3_pos", 
    `mpo-anca` == 1 ~ "MPO_pos")) 
# Make sure order of assays columns and colData rows is the same!!
plasma <- SummarizedExperiment(assays = list(npx = npx), colData = metadata[colnames(npx),])

# correct classes
colData(plasma)$age <- as.numeric(colData(plasma)$age)
colData(plasma)$crp <- as.numeric(colData(plasma)$crp)
colData(plasma)$sr <- as.numeric(colData(plasma)$sr)
colData(plasma)$esr <- as.numeric(colData(plasma)$esr)
colData(plasma)$das28 <- as.numeric(colData(plasma)$das28)
colData(plasma)$time_in_freezer <- as.numeric(colData(plasma)$time_in_freezer)
colData(plasma)$sex <- as.factor(colData(plasma)$sex)
colData(plasma)$cortisone <- as.factor(colData(plasma)$cortisone)
colData(plasma)$diagnosis <- as.factor(colData(plasma)$diagnosis)
colData(plasma)$disease <- as.factor(colData(plasma)$disease)
colData(plasma)$group <- as.factor(colData(plasma)$group)
colData(plasma)$`pr3-anca` <- as.factor(colData(plasma)$`pr3-anca`)
colData(plasma)$`mpo-anca` <- as.factor(colData(plasma)$`mpo-anca`)
colData(plasma)$abs <- as.factor(colData(plasma)$abs)

# Create a Results directory in the top directory
plasma_qc_dir <- "Results/Plasma/QC/"
dir.create(plasma_qc_dir, recursive = TRUE)

olink_dist(plasma_npx, "Olink CARDIOVASCULAR III", paste0(plasma_qc_dir, "plasma_CVIII_dist.png"), width = 12)
olink_dist(plasma_npx, "Olink INFLAMMATION", paste0(plasma_qc_dir, "plasma_INF_dist.png"), width = 12)

plasma_qc <- OlinkAnalyze::olink_qc_plot(plasma_npx)
ggsave(paste0(plasma_qc_dir, "plasma_qc.png"), plot = plasma_qc)

#One sample sample (VASKA637) is filtered out due to technical concerns (failed sample according to Olink QC rapport) and another one due to relatively high level of NA values (ra1622).
plasma <- plasma[,!colnames(plasma) %in% c("vaska637", "ra1622")]

# Impute 1 NA assays - better than throwing all sample away
assay(plasma) <- t(apply(assay(plasma), 1, function(x) ifelse(is.na(x), median(x, na.rm=T), x)))

# PCA - centering, no scaling
pcaplot(plasma, 
        intgroup = "group", 
        ellipse = FALSE, 
        text_labels = FALSE, 
        ntop = 100, 
        pcX = 1, 
        pcY = 2)
pcaobj_plasma <- prcomp(t(assay(plasma)))

#Screeplot
pcascree(pcaobj_plasma,type="pev",
         title="Proportion of explained proportion of variance - plasma",
         pc_nr = 20)

# Correlate PCs
res_pca_plasma <- correlatePCs(pcaobj_plasma,colData(plasma)[,c("age", "sex", "group", "creatinine", "ckd_epi", "time_in_freezer")])
plotPCcorrs(res_pca_plasma)

# hi loadings
hi_loadings(pcaobj_plasma,topN = 10)


###########################
### Heatmap ###############
###########################

keep <- names(head(sort(apply(assay(plasma), 1, sd), decreasing = TRUE), 100))
hm_data <- t(scale(t(assay(plasma)[keep,])))

ha <- HeatmapAnnotation(group = colData(plasma)$group, col = list(group = c("healthy_controls" = "green",
                                                                         "active_disease" = "red",
                                                                         "aav_remission" = "yellow",
                                                                         "RA" = "blue",
                                                                         "SLE_nephritis" = "pink")))

pdf(paste0(plasma_qc_dir, "Heatmap_AllGroups_100_Variable_Prots.pdf"))
hm <- Heatmap(hm_data, 
        top_annotation = ha,
        name = "Z-score",
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 5),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete"
        )
print(hm)
dev.off()

#########################
### Univariate Analysis #
#########################

#test <- assay(plasma) %>% t() %>% cbind(colData(plasma) %>% as.data.frame() %>% select(group, age, ckd_epi))
#mod1 <- aov(IL6 ~ age + ckd_epi + group, data = test)
#mod1 %>% emmeans(pairwise ~ "group") %>% purrr::pluck("contrasts")


#Create olink-function compatible dataset; reduced to the summarizedExperiment samples for the statistical tests
obj <- plasma_npx %>% 
  filter(SampleID %in% colnames(plasma)) %>% 
  left_join(tibble::rownames_to_column(as.data.frame(colData(plasma))), by = c("SampleID" = "rowname"))

# Collect samplenames for later comparison selecting of names in heatmap and multivar
sampSel <- obj %>% 
  dplyr::select(SampleID,group, abs, diagnosis) %>% 
  distinct() %>% 
  tidyr::pivot_longer(cols = -SampleID)
  

# Group comparisons
grps <- OlinkAnalyze::olink_anova_posthoc(df = obj, 
                                              variable = "group", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "group",
                                              verbose = FALSE)
# GPA - MPA
gpa_mpa <- obj %>% 
  filter(diagnosis %in% c("MPA", "GPA")) %>%
  OlinkAnalyze::olink_anova_posthoc(variable = "diagnosis", 
                                    covariates = c("age", "ckd_epi"),
                                    effect = "diagnosis",
                                    verbose = FALSE)

# pr3 - mpo
pr3_mpo <- obj %>%
  filter(abs %in% c("PR3_pos", "MPO_pos")) %>%
  OlinkAnalyze::olink_anova_posthoc(variable = "abs", 
                                    covariates = c("age", "ckd_epi"),
                                    effect = "abs",
                                    verbose = FALSE)

# Collect all results
univariate_results <- rbind(grps, gpa_mpa, pr3_mpo)

for(comp in unique(univariate_results$contrast)){
  
# Create a Results directory in the top directory 
compname <- gsub(" - ", "_vs_", comp)
plasma_uni_dir <- paste0("Results/Plasma/Comparisons/", compname, "/Univariate/")
dir.create(plasma_uni_dir, recursive = TRUE)

df_sub <- univariate_results %>% filter(contrast == comp)
p1 <- ggplot(df_sub, aes(x=estimate, y = -log10(Adjusted_pval))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) + 
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave(paste0(plasma_uni_dir, compname, "_volcano_anova_posthoc.pdf"), plot = p1)

tab <- df_sub %>%
  as.data.frame() %>%
  dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
  dplyr::slice_head(n=15) %>%
  gt::gt() %>%
  gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
gt::gtsave(tab, paste0(plasma_uni_dir, compname, "_Top15_table.pdf"))

openxlsx::write.xlsx(df_sub, paste0(plasma_uni_dir, compname, "_anova_posthoc_results.xlsx"))

# GSEA
plasma_gsea_dir <- paste0("Results/Plasma/Comparisons/", compname, "/Univariate/GSEA/")
dir.create(plasma_gsea_dir, recursive = TRUE)
gsea_out <- gsea_olink_run(x=df_sub, gs=gsea_sets, save=TRUE, saveFolder = plasma_gsea_dir, contrastName = compname)
gsea_olink_viz(gsea_res = gsea_out, saveFolder = plasma_gsea_dir, nterms = 10, contrastName = compname)


# Heatmap
# collect samples
hm_samps <- sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% pull(SampleID)
# Collect prots
keep_sig <- df_sub %>% dplyr::filter(Adjusted_pval < 0.05) %>% dplyr::pull(Assay)
# Extract data
hm_data <- t(scale(t(assay(plasma)[keep_sig,hm_samps ])))

anots <- sampSel %>% 
  dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% 
  left_join(data.frame(name=hm_samps), by = c("SampleID" = "name")) %>% 
  pull(value) %>%
  droplevels()
ha <- HeatmapAnnotation(group = anots, col = list(group = c("healthy_controls" = "green",
                                                            "active_disease" = "red",
                                                            "aav_remission" = "yellow",
                                                            "RA" = "blue",
                                                            "SLE_nephritis" = "pink",
                                                            "MPO_pos" = "green",
                                                            "PR3_pos" = "red",
                                                            "GPA" = "green",
                                                            "MPA" = "red")))

pdf(paste0(plasma_uni_dir, "Heatmap_Sig_Prots_Prots.pdf"))
hm <- Heatmap(hm_data, 
        top_annotation = ha,
        name = "Z-score",
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 5),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete"
)
print(hm)
dev.off()

}


##############################
#### Network enrichment ######
##############################

stringdb <- readr::read_delim("https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz", delim = " ") %>% 
  dplyr::filter(combined_score > 800) %>%
  dplyr::select(protein1, protein2)
# Convert to gene symbol
conv <- readr::read_delim("https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz", delim="\t") %>%
  dplyr::select(protein_external_id, preferred_name)
net <- stringdb %>% 
  left_join(conv, by = c("protein1" = "protein_external_id")) %>% 
  rename(protein_A = preferred_name)%>%
  left_join(conv, by = c("protein2" = "protein_external_id")) %>% 
  rename(protein_B = preferred_name) %>%
  dplyr::select(protein_A, protein_B) %>%
  as.matrix()

gsea_sets_for_nea <- clusterProfiler::read.gmt("inst/extdata/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_symbol.gmt")
fgs <- gsea_sets_for_nea[grep("MSIGDB", gsea_sets_for_nea$term),]
fgs$term <- droplevels(fgs$term)
fgs <- lapply(split(fgs, f = fgs$term), function(x) x$gene)
names(fgs) <- gsub("%.*$", "", names(fgs))

sigList_up <- lapply(split(univariate_results, univariate_results$contrast), function(x) x %>% dplyr::filter(Adjusted_pval < 0.05 & estimate > 0) %>% mutate(Assay = gsub("-", "", Assay)) %>% pull(Assay)) 
sigList_up <- sigList_up[lengths(sigList_up) > 0]
sigList_down <- lapply(split(univariate_results, univariate_results$contrast), function(x) x %>% dplyr::filter(Adjusted_pval < 0.05 & estimate < 0) %>% mutate(Assay = gsub("-", "", Assay)) %>% pull(Assay)) 
sigList_down <- sigList_down[lengths(sigList_down) > 0]

labels <- unique(c(net[,1], net[,2]))
labels_ordered <- sort(labels)
neat.res.up <- neat::neat(alist = sigList_up, blist = fgs, network = net, nettype ="undirected", nodes = labels_ordered, alpha = 0.05)
neat.res.down <- neat::neat(alist = sigList_down, blist = fgs, network = net, nettype ="undirected", nodes = labels_ordered, alpha = 0.05)

for(contrast in unique(univariate_results$contrast)){
  
  # Create a Results directory in the top directory 
  compname <- gsub(" - ", "_vs_", contrast)
  plasma_nea_dir <- paste0("Results/Plasma/Comparisons/", compname, "/Univariate/NEA/")
  dir.create(plasma_nea_dir, recursive = TRUE)
  
  nea_sub_up <- neat.res.up %>% dplyr::filter(A == contrast) %>% dplyr::filter(conclusion == "Overenrichment") %>% mutate(ratio = nab/expected_nab) %>% arrange(-ratio)
  openxlsx::write.xlsx(nea_sub_up, paste0(plasma_nea_dir, compname, "_neat_up_results.xlsx"))

  nea_sub_down <- neat.res.down %>% dplyr::filter(A == contrast) %>% dplyr::filter(conclusion == "Overenrichment") %>% mutate(ratio = nab/expected_nab) %>% arrange(-ratio)
  openxlsx::write.xlsx(nea_sub_down, paste0(plasma_nea_dir, compname, "_neat_down_results.xlsx"))

  }

#############################
##### Correlation Analysis ##
#############################

# Create a Correlation directory in the top directory 
for(marker in c("crp", "sr")){
plasma_marker_corr_dir <- paste0("Results/Plasma/Correlation/", marker, "/")
dir.create(plasma_marker_corr_dir, recursive = TRUE)

# Calculate the pearson correlation
pearson <- obj %>% 
  filter(group == "active_disease") %>% 
  group_by(Assay) %>%
  do(broom::tidy(cor.test(.$NPX, .[[marker]]))) %>%
  dplyr::select(Assay, pearson.cor = estimate, conf.low, conf.high, p.value, method, alternative) %>%
  ungroup()  %>%
  mutate(pval.adj = p.adjust (p.value, method='BH')) %>%
  arrange(pval.adj)
openxlsx::write.xlsx(pearson, paste0(plasma_marker_corr_dir, marker, "_pearson_correlation.xlsx"))
# PLot the linear regression and R squared value per protein
plots <- obj %>% 
  filter(group == "active_disease") %>% 
  group_by(Assay) %>%
  group_modify(~tibble(plots=list(
    ggplot(.) +
      aes_string(x = "NPX", y = marker) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_bw() +
      ggtitle(.y[[1]]) +
      ggpubr::stat_cor(aes(label = ..r.label..), color = "red", geom = "label")
    )))
plasma_marker_corr_plot_dir <- paste0(plasma_marker_corr_dir, "regression_plots/")
dir.create(plasma_marker_corr_plot_dir, recursive = TRUE)
 for(fig in 1:nrow(plots)){
   ggsave(plot=plots$plots[[fig]], paste0(plasma_marker_corr_plot_dir, plots$Assay[[fig]], "_regression.pdf"))
   }
}

#########################
#### Cortisone ##########
#########################

# Create a Results directory in the top directory 
cortisone_res <- obj %>% 
  filter(group == "active_disease") %>% 
  filter(cortisone %in% c("0", "1")) %>% 
  droplevels() %>%
  mutate(cortisone = recode_factor(cortisone, `1` = "Cortisone", `0` = "No_Cortisone")) %>%
  OlinkAnalyze::olink_anova_posthoc(variable = "cortisone", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "cortisone",
                                              verbose = FALSE)

  # Create a Results directory in the top directory 
  compname <- gsub(" - ", "_vs_", unique(cortisone_res$contrast))
  plasma_cort_uni_dir <- paste0("Results/Plasma/Cortisone/",compname, "/Univariate/")
  dir.create(plasma_cort_uni_dir, recursive = TRUE)
  
  p1 <- ggplot(cortisone_res, aes(x=estimate, y = -log10(Adjusted_pval))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) + 
    geom_vline(xintercept = 0) +
    theme_bw()
  ggsave(paste0(plasma_cort_uni_dir, compname, "_volcano_anova_posthoc.pdf"), plot = p1)
  
  tab <- cortisone_res %>%
    as.data.frame() %>%
    dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
    dplyr::slice_head(n=15) %>%
    gt::gt() %>%
    gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
  gt::gtsave(tab, paste0(plasma_cort_uni_dir, compname, "_Top15_table.pdf"))
  
  openxlsx::write.xlsx(cortisone_res, paste0(plasma_cort_uni_dir, compname, "_anova_posthoc_results.xlsx"))
  
  # GSEA
  plasma_cort_gsea_dir <- paste0("Results/Plasma/Cortisone/", compname, "/Univariate/GSEA/")
  dir.create(plasma_cort_gsea_dir, recursive = TRUE)
  gsea_out <- gsea_olink_run(x=cortisone_res, gs=gsea_sets, save=TRUE, saveFolder = plasma_cort_gsea_dir, contrastName = compname)
  gsea_olink_viz(gsea_res = gsea_out, saveFolder = plasma_cort_gsea_dir, nterms = 10, contrastName = compname)

  
###########################
###### BVAS Correlation ###
###########################

# Create a Correlation directory in the top directory 

    plasma_bvas_corr_dir <- paste0("Results/Plasma/Correlation/bvas/")
    dir.create(plasma_bvas_corr_dir, recursive = TRUE)
    
    # Calculate the pearson correlation
    bvas_list <- obj %>% 
      filter(group == "active_disease") %>% 
      group_by(Assay) %>%
      dplyr::filter(bvas > 0) %>%
      do(broom::tidy(cor.test(.$NPX, .$bvas))) %>%
      dplyr::select(Assay, pearson.cor = estimate, conf.low, conf.high, p.value, method, alternative) %>%
      ungroup()  %>%
      mutate(pval.adj = p.adjust (p.value, method='BH')) %>%
      arrange(pval.adj)
    openxlsx::write.xlsx(bvas_list, paste0(plasma_bvas_corr_dir, "bvas_pearson_correlation.xlsx"))
    # PLot the linear regression and R squared value per protein
    plots <- obj %>% 
      filter(group == "active_disease") %>% 
      group_by(Assay) %>%
      dplyr::filter(bvas > 0) %>%
      group_modify(~tibble(plots=list(
        ggplot(.) +
          aes_string(x = "NPX", y = "bvas") +
          geom_point() +
          geom_smooth(method = "lm") +
          theme_bw() +
          ggtitle(.y[[1]]) +
          ggpubr::stat_cor(aes(label = ..r.label..), color = "red", geom = "label")
      )))
    plasma_marker_corr_plot_dir <- paste0(plasma_bvas_corr_dir, "regression_plots/")
    dir.create(plasma_marker_corr_plot_dir, recursive = TRUE)
    for(fig in 1:nrow(plots)){
      ggsave(plot=plots$plots[[fig]], paste0(plasma_marker_corr_plot_dir, plots$Assay[[fig]], "_bvas_correlation.pdf"))
    }
  

########################
#### Multivariate ######
########################

for(comp in unique(univariate_results$contrast)){
# Create a Results directory in the top directory 
compname <- gsub(" - ", "_vs_", comp)
plasma_multi_dir <- paste0("Results/Plasma/Comparisons/", compname, "/Multivariate/")
dir.create(plasma_multi_dir, recursive = TRUE)

pls_obj <- plasma[,sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% pull(SampleID)]

plsda.res <-splsda(t(assay(pls_obj)),factor(colData(pls_obj)[,which.max(sum(colData(pls_obj) %in% strsplit(comp, " - ")[[1]]))]), ncomp=5, scale=TRUE)
plsda.perf <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, nrepeat = 10)

pdf(paste0(plasma_multi_dir, compname, "_ClassificationError_5comp.pdf"))
plot(plsda.perf)
dev.off()

plotIndiv(plsda.res, comp = c(1:2))
ggsave(paste0(plasma_multi_dir, compname, "_SamplePlot_2comp.pdf"))

openxlsx::write.xlsx(list(loadings=plsda.res$loadings$X, scores=plsda.res$variates$X),
                     paste0(plasma_multi_dir, compname, "_Loadings_Scores.xlsx"),
                     row.names=TRUE)

pdf(paste0(plasma_multi_dir, compname, "_Loadings_Comp1_Top15.pdf"))
plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 1)
dev.off()
pdf(paste0(plasma_multi_dir, compname, "_Loadings_Comp2_Top15.pdf"))
plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 2)
dev.off()

# Colect 10 highest absolute loadings
high_c1 <- plsda.res$loadings$X %>%
  as.data.frame() %>%
  arrange(-abs(comp1)) %>%
  head(5) %>% purrr::pluck(rownames)
high_c2 <- plsda.res$loadings$X %>%
  as.data.frame() %>%
  arrange(-abs(comp2)) %>%
  head(5) %>% purrr::pluck(rownames)
high_loads <- c(high_c1, high_c2)


ind.coord = plsda.res$variates$X[, 1:2]
var.coord = t(cor(ind.coord, plsda.res$X))[high_loads,1:2]
bip <- ggplot()+
  geom_point(data = ind.coord, aes(x=comp1, y=comp2, col = plsda.res$Y))+
  geom_segment(data = var.coord, aes(x=0,y=0,xend=comp1*10,yend=comp2*10),
               arrow=arrow(length=unit(0.1,"cm")), color = "#DCDCDC") +
  theme_bw() +
  geom_text(data = var.coord, aes(x=comp1*10, y=comp2*10, label=rownames(var.coord)),color="#006400")
ggsave(paste0(plasma_multi_dir, compname, "_Biplot.pdf"), plot=bip)
}
    
    
    
####################################################################################################################    
############################## Serum ###############################################################################    
####################################################################################################################    
    
# Three proteins (MCP-1, OPG and uPA) are measured in both the Cardiovascular and the Inflammation panel. As seen in the values correlate well between the panels so only the Inflammation data is kept for these proteins. Has to been done before making the summarizedExperiment because pivot_wider otherwise fails.
doubles <- c("MCP-1", "OPG", "uPA")
serum_npx <- serum_npx[!(serum_npx$Assay %in% doubles & serum_npx$Panel == "Olink CARDIOVASCULAR III"),]
    
    # Make a SummarizedExperiment object
    npx <- serum_npx %>% 
      dplyr::select(SampleID, Assay, NPX) %>% 
      tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
      tibble::column_to_rownames("SampleID") %>%
      t()
    metadata <- serum_metadata %>% 
      as.data.frame() %>% 
      tibble::column_to_rownames("id") %>%
      dplyr::rename(ckd_epi = `ckd-epi`) %>%
      mutate(diagnosis = recode_factor(diagnosis, `0` = "MPA", `1` = "GPA")) %>%
      mutate(abs = dplyr::case_when(
        `pr3-anca` == 1 ~ "PR3_pos", 
        `mpo-anca` == 1 ~ "MPO_pos")) 
    # Make sure order of assays columns and colData rows is the same!!
    serum <- SummarizedExperiment(assays = list(npx = npx), colData = metadata[colnames(npx),])
    
    # correct classes
    colData(serum)$age <- as.numeric(colData(serum)$age)
    colData(serum)$crp <- as.numeric(colData(serum)$crp)
    colData(serum)$sr <- as.numeric(colData(serum)$sr)
    colData(serum)$esr <- as.numeric(colData(serum)$esr)
    colData(serum)$das28 <- as.numeric(colData(serum)$das28)
    colData(serum)$time_in_freezer <- as.numeric(colData(serum)$time_in_freezer)
    colData(serum)$sex <- as.factor(colData(serum)$sex)
    colData(serum)$cortisone <- as.factor(colData(serum)$cortisone)
    colData(serum)$diagnosis <- as.factor(colData(serum)$diagnosis)
    colData(serum)$disease <- as.factor(colData(serum)$disease)
    colData(serum)$group <- as.factor(colData(serum)$group)
    colData(serum)$`pr3-anca` <- as.factor(colData(serum)$`pr3-anca`)
    colData(serum)$`mpo-anca` <- as.factor(colData(serum)$`mpo-anca`)
    colData(serum)$abs <- as.factor(colData(serum)$abs)
    
    # Create a Results directory in the top directory 
    serum_qc_dir <- "Results/Serum/QC/"
    dir.create(serum_qc_dir, recursive = TRUE)
    
    olink_dist(serum_npx, "Olink CARDIOVASCULAR III", paste0(serum_qc_dir, "serum_CVIII_dist.png"), width = 12)
    olink_dist(serum_npx, "Olink INFLAMMATION", paste0(serum_qc_dir, "serum_INF_dist.png"), width = 12)
    
    serum_qc <- OlinkAnalyze::olink_qc_plot(serum_npx)
    ggsave(paste0(serum_qc_dir, "serum_qc.png"), plot = serum_qc)
    
    #One sample sample (VASKA637) is filtered out due to technical concerns (failed sample according to Olink QC rapport) and another one due to relatively high level of NA values (ra1622).
    serum <- serum[,!colnames(serum) %in% c("vaska637", "ra1622")]
    
    # PCA - centering, no scaling
    pcaplot(serum, 
            intgroup = "group", 
            ellipse = FALSE, 
            text_labels = FALSE, 
            ntop = 100, 
            pcX = 1, 
            pcY = 2)
    pcaobj_serum <- prcomp(t(assay(serum)))
    
    #Screeplot
    pcascree(pcaobj_serum,type="pev",
             title="Proportion of explained proportion of variance - serum",
             pc_nr = 20)
    
    # Correlate PCs
    res_pca_serum <- correlatePCs(pcaobj_serum,colData(serum)[,c("age", "sex", "group", "creatinine", "ckd_epi", "time_in_freezer")])
    plotPCcorrs(res_pca_serum)
    
    # hi loadings
    hi_loadings(pcaobj_serum,topN = 10)
    
    
    ###########################
    ### Heatmap ###############
    ###########################
    
    keep <- names(head(sort(apply(assay(serum), 1, sd), decreasing = TRUE), 100))
    hm_data <- t(scale(t(assay(serum)[keep,])))
    
    ha <- HeatmapAnnotation(group = colData(serum)$group, col = list(group = c("healthy_controls" = "green",
                                                                                "active_disease" = "red",
                                                                                "aav_remission" = "yellow",
                                                                                "RA" = "blue",
                                                                                "SLE_nephritis" = "pink")))
    
    pdf(paste0(serum_qc_dir, "Heatmap_AllGroups_100_Variable_Prots.pdf"))
    hm <- Heatmap(hm_data, 
                  top_annotation = ha,
                  name = "Z-score",
                  show_column_names = FALSE,
                  row_names_gp = gpar(fontsize = 5),
                  clustering_distance_rows = "euclidean",
                  clustering_method_rows = "complete",
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "complete"
    )
    print(hm)
    dev.off()
    
    #########################
    ### Univariate Analysis #
    #########################
    
    #test <- assay(serum) %>% t() %>% cbind(colData(serum) %>% as.data.frame() %>% select(group, age, ckd_epi))
    #mod1 <- aov(IL6 ~ age + ckd_epi + group, data = test)
    #mod1 %>% emmeans(pairwise ~ "group") %>% purrr::pluck("contrasts")
    
    
    #Create olink-function compatible dataset; reduced to the summarizedExperiment samples for the statistical tests
    obj <- serum_npx %>% 
      filter(SampleID %in% colnames(serum)) %>% 
      left_join(tibble::rownames_to_column(as.data.frame(colData(serum))), by = c("SampleID" = "rowname"))
    
    # Collect samplenames for later comparison selecting of names in heatmap and multivar
    sampSel <- obj %>% 
      dplyr::select(SampleID,group, abs, diagnosis) %>% 
      distinct() %>% 
      tidyr::pivot_longer(cols = -SampleID)
    
    
    # Group comparisons
    grps <- OlinkAnalyze::olink_anova_posthoc(df = obj, 
                                              variable = "group", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "group",
                                              verbose = FALSE)
    # GPA - MPA
    gpa_mpa <- obj %>% 
      filter(diagnosis %in% c("MPA", "GPA")) %>%
      OlinkAnalyze::olink_anova_posthoc(variable = "diagnosis", 
                                        covariates = c("age", "ckd_epi"),
                                        effect = "diagnosis",
                                        verbose = FALSE)
    
    # pr3 - mpo
    pr3_mpo <- obj %>%
      filter(abs %in% c("PR3_pos", "MPO_pos")) %>%
      OlinkAnalyze::olink_anova_posthoc(variable = "abs", 
                                        covariates = c("age", "ckd_epi"),
                                        effect = "abs",
                                        verbose = FALSE)
    
    # Collect all results
    univariate_results <- rbind(grps, gpa_mpa, pr3_mpo)
    
    for(comp in unique(univariate_results$contrast)){
      
      # Create a Results directory in the top directory 
      compname <- gsub(" - ", "_vs_", comp)
      serum_uni_dir <- paste0("Results/Serum/Comparisons/", compname, "/Univariate/")
      dir.create(serum_uni_dir, recursive = TRUE)
      
      df_sub <- univariate_results %>% filter(contrast == comp)
      p1 <- ggplot(df_sub, aes(x=estimate, y = -log10(Adjusted_pval))) +
        geom_point() +
        geom_hline(yintercept = -log10(0.05)) + 
        geom_vline(xintercept = 0) +
        theme_bw()
      ggsave(paste0(serum_uni_dir, compname, "_volcano_anova_posthoc.pdf"), plot = p1)
      
      tab <- df_sub %>%
        as.data.frame() %>%
        dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
        dplyr::slice_head(n=15) %>%
        gt::gt() %>%
        gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
      gt::gtsave(tab, paste0(serum_uni_dir, compname, "_Top15_table.pdf"))
      
      openxlsx::write.xlsx(df_sub, paste0(serum_uni_dir, compname, "_anova_posthoc_results.xlsx"))
      
      # GSEA
      serum_gsea_dir <- paste0("Results/Serum/Comparisons/", compname, "/Univariate/GSEA/")
      dir.create(serum_gsea_dir, recursive = TRUE)
      gsea_out <- gsea_olink_run(x=df_sub, gs=gsea_sets, save=TRUE, saveFolder = serum_gsea_dir, contrastName = compname)
      gsea_olink_viz(gsea_res = gsea_out, saveFolder = serum_gsea_dir, nterms = 10, contrastName = compname)
      
      
      # Heatmap
      # collect samples
      hm_samps <- sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% pull(SampleID)
      # Collect prots
      keep_sig <- df_sub %>% dplyr::filter(Adjusted_pval < 0.05) %>% dplyr::pull(Assay)
      # Extract data
      hm_data <- t(scale(t(assay(serum)[keep_sig,hm_samps ])))
      
      anots <- sampSel %>% 
        dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% 
        left_join(data.frame(name=hm_samps), by = c("SampleID" = "name")) %>% 
        pull(value) %>%
        droplevels()
      ha <- HeatmapAnnotation(group = anots, col = list(group = c("healthy_controls" = "green",
                                                                  "active_disease" = "red",
                                                                  "aav_remission" = "yellow",
                                                                  "SLE_nephritis" = "pink",
                                                                  "MPO_pos" = "green",
                                                                  "PR3_pos" = "red",
                                                                  "GPA" = "green",
                                                                  "MPA" = "red")))
      
      pdf(paste0(serum_uni_dir, "Heatmap_Sig_Prots_Prots.pdf"))
      hm <- Heatmap(hm_data, 
                    top_annotation = ha,
                    name = "Z-score",
                    show_column_names = FALSE,
                    row_names_gp = gpar(fontsize = 5),
                    clustering_distance_rows = "euclidean",
                    clustering_method_rows = "complete",
                    clustering_distance_columns = "euclidean",
                    clustering_method_columns = "complete"
      )
      print(hm)
      dev.off()
      
    }
    
    
    ##############################
    #### Network enrichment ######
    ##############################
   
    # 
    
    sigList_up <- lapply(split(univariate_results, univariate_results$contrast), function(x) x %>% dplyr::filter(Adjusted_pval < 0.05 & estimate > 0) %>% mutate(Assay = gsub("-", "", Assay)) %>% pull(Assay)) 
    sigList_up <- sigList_up[lengths(sigList_up) > 0]
    sigList_down <- lapply(split(univariate_results, univariate_results$contrast), function(x) x %>% dplyr::filter(Adjusted_pval < 0.05 & estimate < 0) %>% mutate(Assay = gsub("-", "", Assay)) %>% pull(Assay)) 
    sigList_down <- sigList_down[lengths(sigList_down) > 0]
    
    labels <- unique(c(net[,1], net[,2]))
    labels_ordered <- sort(labels)
    neat.res.up <- neat::neat(alist = sigList_up, blist = fgs, network = net, nettype ="undirected", nodes = labels_ordered, alpha = 0.05)
    neat.res.down <- neat::neat(alist = sigList_down, blist = fgs, network = net, nettype ="undirected", nodes = labels_ordered, alpha = 0.05)
    
    for(contrast in unique(univariate_results$contrast)){
      
      # Create a Results directory in the top directory 
      compname <- gsub(" - ", "_vs_", contrast)
      serum_nea_dir <- paste0("Results/Serum/Comparisons/", compname, "/Univariate/NEA/")
      dir.create(serum_nea_dir, recursive = TRUE)
      
      nea_sub_up <- neat.res.up %>% dplyr::filter(A == contrast) %>% dplyr::filter(conclusion == "Overenrichment") %>% mutate(ratio = nab/expected_nab) %>% arrange(-ratio)
      openxlsx::write.xlsx(nea_sub_up, paste0(serum_nea_dir, compname, "_neat_up_results.xlsx"))
      
      nea_sub_down <- neat.res.down %>% dplyr::filter(A == contrast) %>% dplyr::filter(conclusion == "Overenrichment") %>% mutate(ratio = nab/expected_nab) %>% arrange(-ratio)
      openxlsx::write.xlsx(nea_sub_down, paste0(serum_nea_dir, compname, "_neat_down_results.xlsx"))
      
    }
    
    #############################
    ##### Correlation Analysis ##
    #############################
    
    # Create a Correlation directory in the top directory 
    for(marker in c("crp", "sr")){
      serum_marker_corr_dir <- paste0("Results/Serum/Correlation/", marker, "/")
      dir.create(serum_marker_corr_dir, recursive = TRUE)
      
      # Calculate the pearson correlation
      pearson <- obj %>% 
        filter(group == "active_disease") %>% 
        group_by(Assay) %>%
        do(broom::tidy(cor.test(.$NPX, .[[marker]]))) %>%
        dplyr::select(Assay, pearson.cor = estimate, conf.low, conf.high, p.value, method, alternative) %>%
        ungroup()  %>%
        mutate(pval.adj = p.adjust (p.value, method='BH')) %>%
        arrange(pval.adj)
      openxlsx::write.xlsx(pearson, paste0(serum_marker_corr_dir, marker, "_pearson_correlation.xlsx"))
      # PLot the linear regression and R squared value per protein
      plots <- obj %>% 
        filter(group == "active_disease") %>% 
        group_by(Assay) %>%
        group_modify(~tibble(plots=list(
          ggplot(.) +
            aes_string(x = "NPX", y = marker) +
            geom_point() +
            geom_smooth(method = "lm") +
            theme_bw() +
            ggtitle(.y[[1]]) +
            ggpubr::stat_cor(aes(label = ..r.label..), color = "red", geom = "label")
        )))
      serum_marker_corr_plot_dir <- paste0(serum_marker_corr_dir, "regression_plots/")
      dir.create(serum_marker_corr_plot_dir, recursive = TRUE)
      for(fig in 1:nrow(plots)){
        ggsave(plot=plots$plots[[fig]], paste0(serum_marker_corr_plot_dir, plots$Assay[[fig]], "_regression.pdf"))
      }
    }
    
    #########################
    #### Cortisone ##########
    #########################
    
    # Create a Results directory in the top directory 
    cortisone_res <- obj %>% 
      filter(group == "active_disease") %>% 
      filter(cortisone %in% c("0", "1")) %>% 
      droplevels() %>%
      mutate(cortisone = recode_factor(cortisone, `1` = "Cortisone", `0` = "No_Cortisone")) %>%
      OlinkAnalyze::olink_anova_posthoc(variable = "cortisone", 
                                        covariates = c("age", "ckd_epi"),
                                        effect = "cortisone",
                                        verbose = FALSE)
    
    # Create a Results directory in the top directory 
    compname <- gsub(" - ", "_vs_", unique(cortisone_res$contrast))
    serum_cort_uni_dir <- paste0("Results/Serum/Cortisone/",compname, "/Univariate/")
    dir.create(serum_cort_uni_dir, recursive = TRUE)
    
    p1 <- ggplot(cortisone_res, aes(x=estimate, y = -log10(Adjusted_pval))) +
      geom_point() +
      geom_hline(yintercept = -log10(0.05)) + 
      geom_vline(xintercept = 0) +
      theme_bw()
    ggsave(paste0(serum_cort_uni_dir, compname, "_volcano_anova_posthoc.pdf"), plot = p1)
    
    tab <- cortisone_res %>%
      as.data.frame() %>%
      dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
      dplyr::slice_head(n=15) %>%
      gt::gt() %>%
      gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
    gt::gtsave(tab, paste0(serum_cort_uni_dir, compname, "_Top15_table.pdf"))
    
    openxlsx::write.xlsx(cortisone_res, paste0(serum_cort_uni_dir, compname, "_anova_posthoc_results.xlsx"))
    
    # GSEA
    serum_cort_gsea_dir <- paste0("Results/Serum/Cortisone/", compname, "/Univariate/GSEA/")
    dir.create(serum_cort_gsea_dir, recursive = TRUE)
    gsea_out <- gsea_olink_run(x=cortisone_res, gs=gsea_sets, save=TRUE, saveFolder = serum_cort_gsea_dir, contrastName = compname)
    gsea_olink_viz(gsea_res = gsea_out, saveFolder = serum_cort_gsea_dir, nterms = 10, contrastName = compname)
    
    
    ###########################
    ###### BVAS Correlation ###
    ###########################
    
    # Create a Correlation directory in the top directory 
    
    serum_bvas_corr_dir <- paste0("Results/Serum/Correlation/bvas/")
    dir.create(serum_bvas_corr_dir, recursive = TRUE)
    
    # Calculate the pearson correlation
    bvas_list <- obj %>% 
      filter(group == "active_disease") %>% 
      group_by(Assay) %>%
      dplyr::filter(bvas > 0) %>%
      do(broom::tidy(cor.test(.$NPX, .$bvas))) %>%
      dplyr::select(Assay, pearson.cor = estimate, conf.low, conf.high, p.value, method, alternative) %>%
      ungroup()  %>%
      mutate(pval.adj = p.adjust (p.value, method='BH')) %>%
      arrange(pval.adj)
    openxlsx::write.xlsx(bvas_list, paste0(serum_bvas_corr_dir, "bvas_pearson_correlation.xlsx"))
    # PLot the linear regression and R squared value per protein
    plots <- obj %>% 
      filter(group == "active_disease") %>% 
      group_by(Assay) %>%
      dplyr::filter(bvas > 0) %>%
      group_modify(~tibble(plots=list(
        ggplot(.) +
          aes_string(x = "NPX", y = "bvas") +
          geom_point() +
          geom_smooth(method = "lm") +
          theme_bw() +
          ggtitle(.y[[1]]) +
          ggpubr::stat_cor(aes(label = ..r.label..), color = "red", geom = "label")
      )))
    serum_marker_corr_plot_dir <- paste0(serum_bvas_corr_dir, "regression_plots/")
    dir.create(serum_marker_corr_plot_dir, recursive = TRUE)
    for(fig in 1:nrow(plots)){
      ggsave(plot=plots$plots[[fig]], paste0(serum_marker_corr_plot_dir, plots$Assay[[fig]], "_bvas_correlation.pdf"))
    }
    
    
    ########################
    #### Multivariate ######
    ########################
    
    for(comp in unique(univariate_results$contrast)){
      # Create a Results directory in the top directory 
      compname <- gsub(" - ", "_vs_", comp)
      serum_multi_dir <- paste0("Results/Serum/Comparisons/", compname, "/Multivariate/")
      dir.create(serum_multi_dir, recursive = TRUE)
      
      pls_obj <- serum[,sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% pull(SampleID)]
      
      plsda.res <-splsda(t(assay(pls_obj)),factor(colData(pls_obj)[,which.max(sum(colData(pls_obj) %in% strsplit(comp, " - ")[[1]]))]), ncomp=5, scale=TRUE)
      plsda.perf <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, nrepeat = 10)
      
      pdf(paste0(serum_multi_dir, compname, "_ClassificationError_5comp.pdf"))
      plot(plsda.perf)
      dev.off()
      
      plotIndiv(plsda.res, comp = c(1:2))
      ggsave(paste0(serum_multi_dir, compname, "_SamplePlot_2comp.pdf"))
      
      openxlsx::write.xlsx(list(loadings=plsda.res$loadings$X, scores=plsda.res$variates$X),
                           paste0(serum_multi_dir, compname, "_Loadings_Scores.xlsx"),
                           row.names=TRUE)
      
      pdf(paste0(serum_multi_dir, compname, "_Loadings_Comp1_Top15.pdf"))
      plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 1)
      dev.off()
      pdf(paste0(serum_multi_dir, compname, "_Loadings_Comp2_Top15.pdf"))
      plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 2)
      dev.off()
      
      # Collect 10 highest absolute loadings
      high_c1 <- plsda.res$loadings$X %>%
        as.data.frame() %>%
        arrange(-abs(comp1)) %>%
        head(5) %>% purrr::pluck(rownames)
      high_c2 <- plsda.res$loadings$X %>%
        as.data.frame() %>%
        arrange(-abs(comp2)) %>%
        head(5) %>% purrr::pluck(rownames)
      high_loads <- c(high_c1, high_c2)
      
      
      ind.coord = plsda.res$variates$X[, 1:2]
      var.coord = t(cor(ind.coord, plsda.res$X))[high_loads,1:2]
      bip <- ggplot()+
        geom_point(data = ind.coord, aes(x=comp1, y=comp2, col = plsda.res$Y))+
        geom_segment(data = var.coord, aes(x=0,y=0,xend=comp1*10,yend=comp2*10),
                     arrow=arrow(length=unit(0.1,"cm")), color = "#DCDCDC") +
        theme_bw() +
        geom_text(data = var.coord, aes(x=comp1*10, y=comp2*10, label=rownames(var.coord)),color="#006400")
      ggsave(paste0(serum_multi_dir, compname, "_Biplot.pdf"), plot=bip)
    }