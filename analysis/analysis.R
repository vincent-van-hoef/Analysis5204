# Clear global env from earlier session and remove result folder for a clean run when sourcing script - comment out if not necessary
rm(list=ls())
setwd("/home/rstudio/")
unlink("Results", recursive=TRUE)

library("Analysis5204")
data(list=c("plasma_metadata", "plasma_npx", "serum_metadata", "serum_npx", "gsea_sets", "lum_dat"), package = "Analysis5204")
library("ggplot2")
library("dplyr")
library("SummarizedExperiment")
library("pcaExplorer")
library("emmeans")
library("ComplexHeatmap")
library("mixOmics")

################################################################################
# Raw data files are loaded automatically when loading the Analysis5204 package.
################################################################################

# First analyze plasma only

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
  dplyr::rename(SampleID = Sample.ID) %>%
  left_join(plasma_npx %>% dplyr::select(SampleID, Index) %>% distinct()) %>%
  dplyr::select(SampleID, Index, OlinkID,UniProt, Assay,MissingFreq , Panel,Panel_Version, PlateID,QC_Warning, LOD,  NPX, Normalization)

plasma_npx <- rbind(plasma_npx, lum_clean)

# Make a SummarizedExperiment object
npx <- plasma_npx %>% 
  dplyr::select(SampleID, Assay, NPX) %>% 
  tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
  tidyr::unchop(everything()) %>%
  distinct() %>%
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

#One sample sample (VASKA637) is filtered out due to technical concerns (failed sample according to Olink QC rapport) and another one due to relatively high level of NA values (ra1622). Umu63 has a s for GPA indication?
plasma <- plasma[,!colnames(plasma) %in% c("vaska637", "ra1622", "umu63")]

# Save metadata
write.csv2(colData(plasma), "Results/Plasma/plasma_metadata.csv", row.names = TRUE)

# Impute 1 NA assays - better than throwing all sample away
assay(plasma) <- t(apply(assay(plasma), 1, function(x) ifelse(is.na(x), median(x, na.rm=T), x)))

# PCA - centering, no scaling
pl_pca <- pcaplot(plasma, 
        intgroup = c("group"), 
        ellipse = FALSE, 
        text_labels = FALSE, 
        ntop = 100, 
        pcX = 1, 
        pcY = 2)
ggsave(paste0(plasma_qc_dir, "disease_pca_100.png"), plot = pl_pca)

pl_pca_13 <- pcaplot(plasma, 
                  intgroup = "group", 
                  ellipse = FALSE, 
                  text_labels = FALSE, 
                  ntop = 100, 
                  pcX = 1, 
                  pcY = 3)
ggsave(paste0(plasma_qc_dir, "plasma_pca_100_pc13.png"), plot = pl_pca_13)


# PCA only actives - centering, no scaling
plasma_aav_rem <- plasma[,colData(plasma)$group %in% c("aav_remission", "active_disease")]
pl_pca <- pcaplot(plasma_aav_rem, 
                  intgroup = c("group", "diagnosis"), 
                  ellipse = FALSE, 
                  text_labels = FALSE, 
                  ntop = 100, 
                  pcX = 1, 
                  pcY = 2)
ggsave(paste0(plasma_qc_dir, "plasma_aav_pca_100.png"), plot = pl_pca)

pl_pca_13 <- pcaplot(plasma_aav_rem, 
                     intgroup = c("group", "diagnosis"), 
                     ellipse = FALSE, 
                     text_labels = FALSE, 
                     ntop = 100, 
                     pcX = 1, 
                     pcY = 3)
ggsave(paste0(plasma_qc_dir, "plasma_aav_pca_100_pc13.png"), plot = pl_pca_13)


pcaobj_plasma <- prcomp(t(assay(plasma)))

#Screeplot
pl_pca_scree <- pcascree(pcaobj_plasma,type="pev",
         title="Proportion of explained proportion of variance - plasma",
         pc_nr = 20)
ggsave(paste0(plasma_qc_dir, "plasma_pca_scree.png"), plot = pl_pca_scree)

# Correlate PCs
res_pca_plasma <- correlatePCs(pcaobj_plasma,colData(plasma)[,c("age", "sex", "group", "creatinine", "ckd_epi", "time_in_freezer")])
plotPCcorrs(res_pca_plasma)

# hi loadings
png(paste0(plasma_qc_dir, "plasma_pca_pc1_loadings.png"))
hi_loadings(pcaobj_plasma,topN = 10)
dev.off()

# PCA of luminex - olink
df_pca <- t(scale(t(assay(plasma))))
df_prcomp <- prcomp(df_pca, center = FALSE, scale. = FALSE)
png(paste0(plasma_qc_dir, "plasma_olink_lum_pca.png"))
plot(df_prcomp$x, col = ifelse(rownames(df_prcomp$x) %in% c("CCL18_PARC", "TIMP-1", "C5a", "CA15-3_MUC-1"), "red", "grey"), pch=19)
dev.off()


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

png(paste0(plasma_qc_dir, "Heatmap_AllGroups_100_Variable_Prots.png"))
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

#Create olink-function compatible dataset; reduced to the summarizedExperiment samples for the statistical tests
obj <- plasma_npx %>% 
  filter(SampleID %in% colnames(plasma)) %>% 
  left_join(tibble::rownames_to_column(as.data.frame(colData(plasma))), by = c("SampleID" = "rowname")) %>%
  mutate(phenoGroups = paste(group, diagnosis, sep = "_")) %>%
  mutate(phenoGroups = gsub("_NA", "", phenoGroups))

# Add pro and mpo to phenogroups
pro <- subset(obj, pr3.anca == 1)
pro$phenoGroups <- paste(pro$group, "PR3", sep = "_")
mpo <- subset(obj, mpo.anca == 1)
mpo$phenoGroups <- paste(mpo$group, "MPO", sep = "_")

obj <- rbind(obj, mpo, pro)

# Collect samplenames for later comparison selecting of names in heatmap and multivar
sampSel <- obj %>% 
  dplyr::select(SampleID, phenoGroups) %>% 
  distinct() %>% 
  tidyr::pivot_longer(cols = -SampleID)
  

# Group comparisons
grps <- OlinkAnalyze::olink_anova_posthoc(df = obj, 
                                              variable = "phenoGroups", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "phenoGroups",
                                              verbose = FALSE)


# Collect all results
comps <- c(
  "active_disease_GPA - healthy_controls",
  "active_disease_PR3 - healthy_controls",
  "active_disease_MPA - healthy_controls",
  "active_disease_MPO - healthy_controls",
  
  "active_disease_GPA - active_disease_MPA",
  "aav_remission_GPA - aav_remission_MPA",
  "active_disease_MPO - active_disease_PR3",
  "aav_remission_MPO - aav_remission_PR3",
  
  "active_disease_GPA - SLE_nephritis",
  "active_disease_PR3 - SLE_nephritis",
  "active_disease_MPA - SLE_nephritis",
  "active_disease_MPO - SLE_nephritis",

  "active_disease_GPA - RA",
  "active_disease_PR3 - RA",
  "active_disease_MPA - RA",
  "active_disease_MPO - RA",
  
  "aav_remission_GPA - active_disease_GPA",
  "aav_remission_PR3 - active_disease_PR3",
  "aav_remission_MPA - active_disease_MPA",
  "aav_remission_MPO - active_disease_MPO"
)

univariate_results <- grps %>% filter(contrast %in% comps)

for(comp in unique(univariate_results$contrast)){
  
# Create a Results directory in the top directory 
compname <- gsub(" - ", "_vs_", comp)
plasma_uni_dir <- paste0("Results/Plasma/Comparisons/", compname, "/Univariate/")
dir.create(plasma_uni_dir, recursive = TRUE)

df_sub <- univariate_results %>% filter(contrast == comp) %>% mutate(Adjusted_pval = ifelse(Adjusted_pval == 0, 1e-12, Adjusted_pval))
p1 <- ggplot(df_sub, aes(x=estimate, y = -log10(Adjusted_pval))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) + 
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave(paste0(plasma_uni_dir, compname, "_volcano_anova_posthoc.png"), plot = p1)

tab <- df_sub %>%
  as.data.frame() %>%
  dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
  dplyr::slice_head(n=15) %>%
  gt::gt() %>%
  gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
gt::gtsave(tab, paste0(plasma_uni_dir, compname, "_Top15_table.png"))

openxlsx::write.xlsx(df_sub, paste0(plasma_uni_dir, compname, "_anova_posthoc_results.xlsx"))

# GSEA
plasma_gsea_dir <- paste0("Results/Plasma/Comparisons/", compname, "/Univariate/GSEA/")
dir.create(plasma_gsea_dir, recursive = TRUE)
gsea_out <- gsea_olink_run(x=df_sub, gs=gsea_sets, save=TRUE, saveFolder = plasma_gsea_dir, contrastName = compname)
gsea_olink_viz(gsea_res = gsea_out, saveFolder = plasma_gsea_dir, nterms = 10, contrastName = compname)


# Heatmap
# collect samples
hm_samps <- sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% distinct(SampleID) %>% pull(SampleID)
# Collect prots
keep_sig <- df_sub %>% dplyr::filter(Adjusted_pval < 0.05) %>% dplyr::pull(Assay)
# Extract data
hm_data <- t(scale(t(assay(plasma)[keep_sig,hm_samps, drop=FALSE ])))

anots <- sampSel %>% 
  dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% 
  distinct(SampleID, .keep_all = TRUE) %>%
  left_join(data.frame(name=hm_samps), by = c("SampleID" = "name")) %>% 
  pull(value) 
ha <- HeatmapAnnotation(group = anots, col = list(group = c("healthy_controls" = "green",
                                                            "active_disease_GPA" = "red",
                                                            "active_disease_MPA" = "orange",
                                                            "active_disease_MPO" = "red",
                                                            "active_disease_PR3" = "orange",
                                                            "aav_remission_GPA" = "yellow",
                                                            "aav_remission_MPA" = "green",
                                                            "aav_remission_MPO" = "yellow",
                                                            "aav_remission_PR3" = "green",
                                                            "RA" = "blue",
                                                            "SLE_nephritis" = "blue"
                                                            )))

png(paste0(plasma_uni_dir, "Heatmap_Sig_Prots_Prots.png"))
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
  dplyr::rename(protein_A = preferred_name)%>%
  left_join(conv, by = c("protein2" = "protein_external_id")) %>% 
  dplyr::rename(protein_B = preferred_name) %>%
  dplyr::select(protein_A, protein_B) %>%
  as.matrix()

gsea_sets_for_nea <- clusterProfiler::read.gmt("inst/extdata/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_symbol.gmt")
fgs <- gsea_sets_for_nea[grep("MSIGDB", gsea_sets_for_nea$term),]
fgs$term <- droplevels(fgs$term)
fgs <- lapply(split(fgs, f = fgs$term), function(x) x$gene)
names(fgs) <- gsub("%.*$", "", names(fgs))

sigList_up <- lapply(split(univariate_results, univariate_results$contrast), function(x) x %>% dplyr::filter(Adjusted_pval < 0.05 & estimate > 0) %>% mutate(Assay = gsub("-", "", Assay)) %>% pull(Assay)) 
sigList_up <- sigList_up[lengths(sigList_up) > 1]
sigList_down <- lapply(split(univariate_results, univariate_results$contrast), function(x) x %>% dplyr::filter(Adjusted_pval < 0.05 & estimate < 0) %>% mutate(Assay = gsub("-", "", Assay)) %>% pull(Assay)) 
sigList_down <- sigList_down[lengths(sigList_down) > 1]

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
for(marker in c("crp", "sr", "bvas")){
plasma_marker_corr_dir <- paste0("Results/Plasma/Correlation/", marker, "/")
dir.create(plasma_marker_corr_dir, recursive = TRUE)

# Calculate the pearson correlation
corrFun <- function(pheno = "active_disease", desc = "all_active_disease") {

pearson <- obj %>% 
  filter(phenoGroups %in% pheno) %>% 
  group_by(Assay) %>%
  do(broom::tidy(cor.test(.$NPX, .[[marker]]))) %>%
  dplyr::select(Assay, pearson.cor = estimate, conf.low, conf.high, p.value, method, alternative) %>%
  ungroup()  %>%
  mutate(pval.adj = p.adjust (p.value, method='BH')) %>%
  arrange(pval.adj)
openxlsx::write.xlsx(pearson, paste0(plasma_marker_corr_dir, marker, "_", desc, "_pearson_correlation.xlsx"))

}

corrFun("active_disease_GPA", "GPA")
corrFun("active_disease_MPA", "MPA")
corrFun("active_disease_MPO", "MPO")
corrFun("active_disease_PR3", "PR3")

# PLot the linear regression and R squared value per protein
corrPlotFun <- function(pheno = "active_disease", desc = "all_active_disease") {
  
plots <- obj %>% 
  filter(phenoGroups %in% pheno) %>% 
  group_by(Assay) %>%
  group_modify(~tibble(plots=list(
    ggplot(., aes(color = phenoGroups)) +
      aes_string(x = "NPX", y = marker) +
      geom_point() +
      geom_smooth(method = "lm", fullrange = TRUE) +
      theme_bw() +
      ggtitle(.y[[1]]) +
      ggpubr::stat_cor(aes(color = phenoGroups, label = ..r.label..), geom = "label")
    )))
plasma_marker_corr_plot_dir <- paste0(plasma_marker_corr_dir, "regression_plots/")
dir.create(plasma_marker_corr_plot_dir, recursive = TRUE)
 for(fig in 1:nrow(plots)){
   ggsave(plot=plots$plots[[fig]], paste0(plasma_marker_corr_plot_dir, plots$Assay[[fig]], "_", desc, "_regression.png"))
 }
}

corrPlotFun(c("active_disease_GPA", "active_disease_MPA"), "GPA_MPA")
corrPlotFun(c("active_disease_MPO", "active_disease_PR3"), "MP_PR3")


}

#########################
#### Cortisone ##########
#########################

# Create a Results directory in the top directory 
cortisone_res <- obj %>% 
  filter(phenoGroups %in% c("active_disease_GPA", "active_disease_MPA", "active_disease_MPO", "active_disease_PR3")) %>% 
  filter(cortisone %in% c("0", "1")) %>% 
  droplevels() %>%
  mutate(cortisone = recode_factor(cortisone, `1` = "Cortisone", `0` = "No_Cortisone")) %>%
  mutate(pheno_cortisone = paste(phenoGroups, cortisone, sep = "_")) %>%
  OlinkAnalyze::olink_anova_posthoc(variable = "pheno_cortisone", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "pheno_cortisone",
                                              verbose = FALSE)

# Collect all desired results
comps <- c(
  "active_disease_GPA_Cortisone - active_disease_GPA_No_Cortisone",
  "active_disease_MPA_Cortisone - active_disease_MPA_No_Cortisone",
  "active_disease_MPO_Cortisone - active_disease_MPO_No_Cortisone",
  "active_disease_PR3_Cortisone - active_disease_PR3_No_Cortisone"
)

cortisone_results <- cortisone_res %>% filter(contrast %in% comps)

for(comp in unique(cortisone_results$contrast)){
  
  # Create a Results directory in the top directory 
  compname <- gsub(" - ", "_vs_", comp)
  plasma_cort_uni_dir <- paste0("Results/Plasma/Cortisone/",compname, "/Univariate/")
  dir.create(plasma_cort_uni_dir, recursive = TRUE)
  
  df_sub <- cortisone_results %>% filter(contrast == comp) %>% mutate(Adjusted_pval = ifelse(Adjusted_pval == 0, 1e-12, Adjusted_pval)) 
  p1 <- ggplot(df_sub, aes(x=estimate, y = -log10(Adjusted_pval))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) + 
    geom_vline(xintercept = 0) +
    theme_bw()
  ggsave(paste0(plasma_cort_uni_dir, compname, "_volcano_anova_posthoc.png"), plot = p1)
  
  tab <- df_sub %>%
    as.data.frame() %>%
    dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
    dplyr::slice_head(n=15) %>%
    gt::gt() %>%
    gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
  gt::gtsave(tab, paste0(plasma_cort_uni_dir, compname, "_Top15_table.png"))
  
  openxlsx::write.xlsx(cortisone_results, paste0(plasma_cort_uni_dir, compname, "_anova_posthoc_results.xlsx"))
  
  # GSEA
  plasma_cort_gsea_dir <- paste0("Results/Plasma/Cortisone/", compname, "/Univariate/GSEA/")
  dir.create(plasma_cort_gsea_dir, recursive = TRUE)
  gsea_out <- gsea_olink_run(x=df_sub, gs=gsea_sets, save=TRUE, saveFolder = plasma_cort_gsea_dir, contrastName = compname)
  gsea_olink_viz(gsea_res = gsea_out, saveFolder = plasma_cort_gsea_dir, nterms = 10, contrastName = compname)
}

########################
#### Multivariate ######
########################

for(comp in unique(univariate_results$contrast)){
# Create a Results directory in the top directory 
compname <- gsub(" - ", "_vs_", comp)
plasma_multi_dir <- paste0("Results/Plasma/Comparisons/", compname, "/Multivariate/")
dir.create(plasma_multi_dir, recursive = TRUE)

# extract samples and data
use_samps <- sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% pull(SampleID)
pls_obj <- plasma[,use_samps]

# Create outcome list
outs <- obj %>% 
  dplyr::filter(phenoGroups %in% strsplit(comp, " - ")[[1]]) %>% 
  filter(SampleID %in% use_samps) %>% distinct(SampleID, .keep_all = TRUE) %>% 
  dplyr::select(SampleID, phenoGroups) %>%
  tibble::column_to_rownames("SampleID")


plsda.res <-plsda(t(assay(pls_obj)),factor(outs[use_samps,]), ncomp=5, scale=TRUE)
plsda.perf <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, nrepeat = 10)

png(paste0(plasma_multi_dir, compname, "_ClassificationError_5comp.png"))
plot(plsda.perf)
dev.off()

# rerun with 2 components - necessary to get scaling right in biplot
plsda.res <-plsda(t(assay(pls_obj)),factor(outs[use_samps,]), ncomp=2, scale=TRUE)

plotIndiv(plsda.res, comp = c(1:2))
ggsave(paste0(plasma_multi_dir, compname, "_SamplePlot_2comp.png"))

openxlsx::write.xlsx(list(loadings=plsda.res$loadings$X, scores=plsda.res$variates$X),
                     paste0(plasma_multi_dir, compname, "_Loadings_Scores.xlsx"),
                     row.names=TRUE)

png(paste0(plasma_multi_dir, compname, "_Loadings_Comp1_Top15.png"))
plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 1)
dev.off()
png(paste0(plasma_multi_dir, compname, "_Loadings_Comp2_Top15.png"))
plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 2)
dev.off()

# Biplot showing vars above certain absolute correlation, weight in alternative axis
var.coord = t(cor(plsda.res$variates$X[, 1:2], plsda.res$X))[,1:2]
protno9 <- sum(apply(var.coord, 1, function(x) any(abs(x) > 0.9)))
protno8 <- sum(apply(var.coord, 1, function(x) any(abs(x) > 0.8)))

if (protno9 > 4) {
  bip09 <- biplot(plsda.res, ind.names = FALSE, cutoff = 0.9, hline = TRUE, vline = TRUE, pch.size = 1, var.names.size = 3)
  ggsave(paste0(plasma_multi_dir, compname, "_Biplot_cor09.png"), plot=bip09)
} else if (protno8 > 4) { 
  bip08 <- biplot(plsda.res, ind.names = FALSE, cutoff = 0.8, hline = TRUE, vline = TRUE, pch.size = 1, var.names.size = 3)
  ggsave(paste0(plasma_multi_dir, compname, "_Biplot_cor08.png"), plot=bip08)
} 

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
    
    
    # Save metadata
    write.csv2(colData(serum), "Results/Serum/serum_metadata.csv", row.names = TRUE)
    
    # PCA - centering, no scaling
    se_pca <- pcaplot(serum, 
            intgroup = "group", 
            ellipse = FALSE, 
            text_labels = FALSE, 
            ntop = 100, 
            pcX = 1, 
            pcY = 2)
    ggsave(paste0(serum_qc_dir, "serum_pca_100.png"), plot = se_pca)
    
    se_pca_13 <- pcaplot(serum, 
                      intgroup = "group", 
                      ellipse = FALSE, 
                      text_labels = FALSE, 
                      ntop = 100, 
                      pcX = 1, 
                      pcY = 3)
    ggsave(paste0(serum_qc_dir, "serum_pca_100_pc13.png"), plot = se_pca_13)
    
    pcaobj_serum <- prcomp(t(assay(serum)))
    
    #Screeplot
    se_pca_scree <- pcascree(pcaobj_serum,type="pev",
             title="Proportion of explained proportion of variance - serum",
             pc_nr = 20)
    ggsave(paste0(serum_qc_dir, "serum_pca_scree.png"), plot = se_pca_scree)
    
    # Correlate PCs
    res_pca_serum <- correlatePCs(pcaobj_serum,colData(serum)[,c("age", "sex", "group", "creatinine", "ckd_epi", "time_in_freezer")])
    plotPCcorrs(res_pca_serum)
    
    # hi loadings
    png(paste0(serum_qc_dir, "serum_pca_pc1_loadings.png"))
    hi_loadings(pcaobj_serum,topN = 10)
    dev.off()
    
    
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
    
    png(paste0(serum_qc_dir, "Heatmap_AllGroups_100_Variable_Prots.png"))
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
      left_join(tibble::rownames_to_column(as.data.frame(colData(serum))), by = c("SampleID" = "rowname")) %>%
      mutate(phenoGroups = paste(group, diagnosis, sep = "_")) %>%
      mutate(phenoGroups = gsub("_NA", "", phenoGroups))
    
    
    # Add pro and mpo to phenogroups
    pro <- subset(obj, pr3.anca == 1)
    pro$phenoGroups <- paste(pro$group, "PR3", sep = "_")
    mpo <- subset(obj, mpo.anca == 1)
    mpo$phenoGroups <- paste(mpo$group, "MPO", sep = "_")
    
    obj <- rbind(obj, mpo, pro)
    
    # Collect samplenames for later comparison selecting of names in heatmap and multivar
    sampSel <- obj %>% 
      dplyr::select(SampleID, phenoGroups) %>% 
      distinct() %>% 
      tidyr::pivot_longer(cols = -SampleID)
    
    
    # Group comparisons
    grps <- OlinkAnalyze::olink_anova_posthoc(df = obj, 
                                              variable = "phenoGroups", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "phenoGroups",
                                              verbose = FALSE)


    
    # Collect all results
    comps <- c(
      "active_disease_GPA - healthy_controls",
      "active_disease_PR3 - healthy_controls",
      "active_disease_MPA - healthy_controls",
      "active_disease_MPO - healthy_controls",
      
      "active_disease_GPA - active_disease_MPA",
      "aav_remission_GPA - aav_remission_MPA",
      "active_disease_MPO - active_disease_PR3",
      "aav_remission_MPO - aav_remission_PR3",
      
      "active_disease_GPA - SLE_nephritis",
      "active_disease_PR3 - SLE_nephritis",
      "active_disease_MPA - SLE_nephritis",
      "active_disease_MPO - SLE_nephritis",
      
      "aav_remission_GPA - active_disease_GPA",
      "aav_remission_PR3 - active_disease_PR3",
      "aav_remission_MPA - active_disease_MPA",
      "aav_remission_MPO - active_disease_MPO"
    )
    
    univariate_results <- grps %>% filter(contrast %in% comps)
        
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
      ggsave(paste0(serum_uni_dir, compname, "_volcano_anova_posthoc.png"), plot = p1)
      
      tab <- df_sub %>%
        as.data.frame() %>%
        dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
        dplyr::slice_head(n=15) %>%
        gt::gt() %>%
        gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
      gt::gtsave(tab, paste0(serum_uni_dir, compname, "_Top15_table.png"))
      
      openxlsx::write.xlsx(df_sub, paste0(serum_uni_dir, compname, "_anova_posthoc_results.xlsx"))
      
      # GSEA
      serum_gsea_dir <- paste0("Results/Serum/Comparisons/", compname, "/Univariate/GSEA/")
      dir.create(serum_gsea_dir, recursive = TRUE)
      gsea_out <- gsea_olink_run(x=df_sub, gs=gsea_sets, save=TRUE, saveFolder = serum_gsea_dir, contrastName = compname)
      gsea_olink_viz(gsea_res = gsea_out, saveFolder = serum_gsea_dir, nterms = 10, contrastName = compname)
      
      
      # Heatmap
      # collect samples
      hm_samps <- sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% distinct(SampleID) %>% pull(SampleID)
      # Collect prots
      keep_sig <- df_sub %>% dplyr::filter(Adjusted_pval < 0.05) %>% dplyr::pull(Assay)
      # Extract data
      hm_data <- t(scale(t(assay(serum)[keep_sig,hm_samps , drop=FALSE])))
      
      anots <- sampSel %>% 
        dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% 
        distinct(SampleID, .keep_all = TRUE) %>%
        left_join(data.frame(name=hm_samps), by = c("SampleID" = "name")) %>% 
        pull(value)
        
      ha <- HeatmapAnnotation(group = anots, col = list(group = c("healthy_controls" = "green",
                                                                  "active_disease_GPA" = "red",
                                                                  "active_disease_MPA" = "orange",
                                                                  "active_disease_MPO" = "red",
                                                                  "active_disease_PR3" = "orange",
                                                                  "aav_remission_GPA" = "yellow",
                                                                  "aav_remission_MPA" = "green",
                                                                  "aav_remission_MPO" = "yellow",
                                                                  "aav_remission_PR3" = "green",
                                                                  "SLE_nephritis" = "blue")))
      
      png(paste0(serum_uni_dir, "Heatmap_Sig_Prots_Prots.png"))
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
    sigList_up <- sigList_up[lengths(sigList_up) > 1]
    sigList_down <- lapply(split(univariate_results, univariate_results$contrast), function(x) x %>% dplyr::filter(Adjusted_pval < 0.05 & estimate < 0) %>% mutate(Assay = gsub("-", "", Assay)) %>% pull(Assay)) 
    sigList_down <- sigList_down[lengths(sigList_down) > 1]
    
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
    for(marker in c("crp", "sr", "bvas")){
      serum_marker_corr_dir <- paste0("Results/Serum/Correlation/", marker, "/")
      dir.create(serum_marker_corr_dir, recursive = TRUE)
      
      # Calculate the pearson correlation
      corrFun <- function(pheno = "active_disease", desc = "all_active_disease") {
        
        pearson <- obj %>% 
          filter(phenoGroups %in% pheno) %>% 
          group_by(Assay) %>%
          do(broom::tidy(cor.test(.$NPX, .[[marker]]))) %>%
          dplyr::select(Assay, pearson.cor = estimate, conf.low, conf.high, p.value, method, alternative) %>%
          ungroup()  %>%
          mutate(pval.adj = p.adjust (p.value, method='BH')) %>%
          arrange(pval.adj)
        openxlsx::write.xlsx(pearson, paste0(serum_marker_corr_dir, marker, "_", desc, "_pearson_correlation.xlsx"))
        
      }
      corrFun("active_disease_GPA", "GPA")
      corrFun("active_disease_MPA", "MPA")
      corrFun("active_disease_MPO", "MPO")
      corrFun("active_disease_PR3", "PR3")
      # PLot the linear regression and R squared value per protein
      corrPlotFun <- function(pheno = "active_disease", desc = "all_active_disease") {
        
        plots <- obj %>% 
          filter(phenoGroups %in% pheno) %>% 
          group_by(Assay) %>%
          group_modify(~tibble(plots=list(
            ggplot(., aes(color = phenoGroups)) +
              aes_string(x = "NPX", y = marker) +
              geom_point() +
              geom_smooth(method = "lm", fullrange = TRUE) +
              theme_bw() +
              ggtitle(.y[[1]]) +
              ggpubr::stat_cor(aes(color = phenoGroups, label = ..r.label..), geom = "label")
          )))
        serum_marker_corr_plot_dir <- paste0(serum_marker_corr_dir, "regression_plots/")
        dir.create(serum_marker_corr_plot_dir, recursive = TRUE)
        for(fig in 1:nrow(plots)){
          ggsave(plot=plots$plots[[fig]], paste0(serum_marker_corr_plot_dir, plots$Assay[[fig]], "_", desc, "_regression.png"))
        }
      }
      
      corrPlotFun(c("active_disease_GPA", "active_disease_MPA"), "GPA_MPA")
      corrPlotFun(c("active_disease_MPO", "active_disease_PR3"), "MP_PR3")
      
    }
    
    #########################
    #### Cortisone ##########
    #########################
    
    
    # Create a Results directory in the top directory 
    cortisone_res <- obj %>% 
      filter(phenoGroups %in% c("active_disease_GPA", "active_disease_MPA", "active_disease_MPO", "active_disease_PR3")) %>% 
      filter(cortisone %in% c("0", "1")) %>% 
      droplevels() %>%
      mutate(cortisone = recode_factor(cortisone, `1` = "Cortisone", `0` = "No_Cortisone")) %>%
      mutate(pheno_cortisone = paste(phenoGroups, cortisone, sep = "_")) %>%
      OlinkAnalyze::olink_anova_posthoc(variable = "pheno_cortisone", 
                                        covariates = c("age", "ckd_epi"),
                                        effect = "pheno_cortisone",
                                        verbose = FALSE)
    
    # Collect all desired results
    comps <- c(
      "active_disease_GPA_Cortisone - active_disease_GPA_No_Cortisone",
      "active_disease_MPA_Cortisone - active_disease_MPA_No_Cortisone",
      "active_disease_MPO_Cortisone - active_disease_MPO_No_Cortisone",
      "active_disease_PR3_Cortisone - active_disease_PR3_No_Cortisone"
    )
    
    cortisone_results <- cortisone_res %>% filter(contrast %in% comps)

    
    for(comp in unique(cortisone_results$contrast)){
      
      # Create a Results directory in the top directory 
      compname <- gsub(" - ", "_vs_", comp)
      serum_cort_uni_dir <- paste0("Results/Serum/Cortisone/",compname, "/Univariate/")
      dir.create(plasma_cort_uni_dir, recursive = TRUE)
      
      df_sub <- cortisone_results %>% filter(contrast == comp) %>% mutate(Adjusted_pval = ifelse(Adjusted_pval == 0, 1e-12, Adjusted_pval)) 
      p1 <- ggplot(df_sub, aes(x=estimate, y = -log10(Adjusted_pval))) +
        geom_point() +
        geom_hline(yintercept = -log10(0.05)) + 
        geom_vline(xintercept = 0) +
        theme_bw()
      ggsave(paste0(serum_cort_uni_dir, compname, "_volcano_anova_posthoc.png"), plot = p1)
      
      tab <- df_sub %>%
        as.data.frame() %>%
        dplyr::select(-(c(OlinkID, UniProt, Panel, term, contrast))) %>%
        dplyr::slice_head(n=15) %>%
        gt::gt() %>%
        gt::tab_header(paste0("Top 15 Differentially Expressed Proteins in ", compname))
      gt::gtsave(tab, paste0(serum_cort_uni_dir, compname, "_Top15_table.png"))
      
      openxlsx::write.xlsx(cortisone_results, paste0(serum_cort_uni_dir, compname, "_anova_posthoc_results.xlsx"))
      
      # GSEA
      serum_cort_gsea_dir <- paste0("Results/Serum/Cortisone/", compname, "/Univariate/GSEA/")
      dir.create(serum_cort_gsea_dir, recursive = TRUE)
      gsea_out <- gsea_olink_run(x=df_sub, gs=gsea_sets, save=TRUE, saveFolder = serum_cort_gsea_dir, contrastName = compname)
      gsea_olink_viz(gsea_res = gsea_out, saveFolder = serum_cort_gsea_dir, nterms = 10, contrastName = compname)
    }
    
    ########################
    #### Multivariate ######
    ########################
    
    for(comp in unique(univariate_results$contrast)){
     
       # Create a Results directory in the top directory 
      compname <- gsub(" - ", "_vs_", comp)
      serum_multi_dir <- paste0("Results/Serum/Comparisons/", compname, "/Multivariate/")
      dir.create(serum_multi_dir, recursive = TRUE)
      
      # extract samples and data
      use_samps <- sampSel %>% dplyr::filter(value %in% strsplit(comp, " - ")[[1]]) %>% pull(SampleID)
      pls_obj <- serum[,use_samps]
      
      # Create outcome list
      outs <- obj %>% 
        dplyr::filter(phenoGroups %in% strsplit(comp, " - ")[[1]]) %>% 
        filter(SampleID %in% use_samps) %>% distinct(SampleID, .keep_all = TRUE) %>% 
        dplyr::select(SampleID, phenoGroups) %>%
        tibble::column_to_rownames("SampleID")
      
      plsda.res <-plsda(t(assay(pls_obj)),factor(outs[use_samps,]), ncomp=5, scale=TRUE)
      plsda.perf <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, nrepeat = 10)
      
      png(paste0(serum_multi_dir, compname, "_ClassificationError_5comp.png"))
      plot(plsda.perf)
      dev.off()
      
      # rerun with 2 components - necessary to get scaling right in biplot
      plsda.res <-plsda(t(assay(pls_obj)),factor(outs[use_samps,]), ncomp=2, scale=TRUE)
      
      plotIndiv(plsda.res, comp = c(1:2))
      ggsave(paste0(serum_multi_dir, compname, "_SamplePlot_2comp.png"))
      
      openxlsx::write.xlsx(list(loadings=plsda.res$loadings$X, scores=plsda.res$variates$X),
                           paste0(serum_multi_dir, compname, "_Loadings_Scores.xlsx"),
                           row.names=TRUE)
      
      png(paste0(serum_multi_dir, compname, "_Loadings_Comp1_Top15.png"))
      plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 1)
      dev.off()
      png(paste0(serum_multi_dir, compname, "_Loadings_Comp2_Top15.png"))
      plotLoadings(plsda.res, contrib = "max", ndisplay = 15, comp = 2)
      dev.off()
      
      
      # Biplot showing vars above certain absolute correlation, weight in alternative axis
      var.coord = t(cor(plsda.res$variates$X[, 1:2], plsda.res$X))[,1:2]
      protno9 <- sum(apply(var.coord, 1, function(x) any(abs(x) > 0.9)))
      protno8 <- sum(apply(var.coord, 1, function(x) any(abs(x) > 0.8)))
      
      if (protno9 > 4) {
        bip09 <- biplot(plsda.res, ind.names = FALSE, cutoff = 0.9, hline = TRUE, vline = TRUE, pch.size = 1, var.names.size = 3)
        ggsave(paste0(serum_multi_dir, compname, "_Biplot_cor09.png"), plot=bip09)
      } else if (protno8 > 4) { 
        bip08 <- biplot(plsda.res, ind.names = FALSE, cutoff = 0.8, hline = TRUE, vline = TRUE, pch.size = 1, var.names.size = 3)
        ggsave(paste0(serum_multi_dir, compname, "_Biplot_cor08.png"), plot=bip08)
      } 
    }
    