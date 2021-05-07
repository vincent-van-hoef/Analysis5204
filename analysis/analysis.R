# Clear global env from earlier session and remove result folder for a clean run when sourcing script - comment out if not necessary
rm(list=ls())
unlink("Results", recursive=TRUE)

library("Analysis5204")
data(list=c("plasma_metadata", "plasma_npx", "serum_metadata", "serum_npx", "gsea_sets"), package = "Analysis5204")
library("ggplot2")
library("dplyr")
library("SummarizedExperiment")
library("pcaExplorer")
library("emmeans")

#################################################################################
# Raw data files are loaded automatically when loading the Analysis5204 package.#
#################################################################################

# First analyze plasma only

# Three proteins (MCP-1, OPG and uPA) are measured in both the Cardiovascular and the Inflammation panel. As seen in the values correlate well between the panels so only the Inflammation data is kept for these proteins. Has to been done before making the summarizedExperiment because pivot_wider otherwise fails.
doubles <- c("MCP-1", "OPG", "uPA")
plasma_npx <- plasma_npx[!(plasma_npx$Assay %in% doubles & plasma_npx$Panel == "Olink CARDIOVASCULAR III"),]

# Make a SummarizedExperiment object
npx <- plasma_npx %>% 
  dplyr::select(SampleID, Assay, NPX) %>% 
  tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
  tibble::column_to_rownames("SampleID") %>%
  t()
metadata <- plasma_metadata %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("id")
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

#########################
### Univariate Analysis #
#########################

test <- assay(plasma) %>% t() %>% cbind(colData(plasma) %>% as.data.frame() %>% select(group, age, ckd_epi))
mod1 <- aov(IL6 ~ age + ckd_epi + group, data = test)
mod1 %>% emmeans(pairwise ~ "group") %>% purrr::pluck("contrasts")

#Create olink-function compatible dataset; reduced to the summarizedExperiment samples for the statistical tests
obj <- plasma_npx %>% 
  filter(SampleID %in% colnames(plasma)) %>% 
  left_join(tibble::rownames_to_column(as.data.frame(colData(plasma))), by = c("SampleID" = "rowname"))

# Group comparisons
grps <- OlinkAnalyze::olink_anova_posthoc(df = obj, 
                                              variable = "group", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "group",
                                              verbose = FALSE)
# GPA - MPA
gpa_mpa <- obj %>% 
  filter(diagnosis %in% c("0", "1")) %>%
  mutate(diagnosis = recode_factor(diagnosis, `0` = "MPA", `1` = "GPA")) %>%
  OlinkAnalyze::olink_anova_posthoc(variable = "diagnosis", 
                                    covariates = c("age", "ckd_epi"),
                                    effect = "diagnosis",
                                    verbose = FALSE)

# pr3 - mpo
pr3_mpo <- obj %>%
 mutate(abs = dplyr::case_when(
    pr3.anca == 1 ~ "PR3_pos", 
    mpo.anca == 1 ~ "MPO_pos")) %>%
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
  select(Assay, pearson.cor = estimate, conf.low, conf.high, p.value, method, alternative) %>%
  arrange(p.value)
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
      stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")
    )))
plasma_marker_corr_plot_dir <- paste0(plasma_marker_corr_dir, "regression_plots/")
dir.create(plasma_marker_corr_plot_dir, recursive = TRUE)
 for(fig in 1:nrow(plots)){
   ggsave(plot=plots$plots[[fig]], paste0(plasma_marker_corr_plot_dir, plots$Assay[[fig]], "_regression.pdf"))
   }
}


# Cortisone effect
# Create a Results directory in the top directory 
plasma_cortisone_dir <- paste0("Results/Plasma/Cortisone/")
dir.create(plasma_cortisone_dir, recursive = TRUE)
obj %>% 
  filter(group == "active_disease") %>% 
  filter(cortisone %in% c("0", "1")) %>% 
  droplevels() %>%
  OlinkAnalyze::olink_anova_posthoc(variable = "cortisone", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "cortisone",
                                              verbose = FALSE) %>%
  openxlsx::write.xlsx(paste0(plasma_cortisone_dir, "Cortisone_Effect_Plasma_ActiveAAV.xlsx"))


########################
#### Multivariate ######
########################

opls_obj <- plasma[,plasma$group %in% c("active_disease", "healthy_controls")]

MyResult.splsda <-splsda(t(assay(opls_obj)), factor(colData(opls_obj)$group))
plotIndiv(MyResult.splsda)
plotVar(MyResult.splsda, var.names=TRUE, cutoff = 0.8, cex = 2)
vars <- selectVar(MyResult.splsda)
plotLoadings(MyResult.splsda, contrib = "max", ndisplay = 15)

testExp <- plasma


colData(testExp) <- colData(testExp) %>% 
  as.data.frame() %>%
  mutate(abs = dplyr::case_when(
    pr3.anca == 1 ~ "PR3_pos", 
    mpo.anca == 1 ~ "MPO_pos")) %>%
  DataFrame()
opls_obj <- testExp[,testExp$abs %in% c("PR3_pos", "MPO_pos")] 
plsda.res <- plsda(t(assay(opls_obj)), factor(colData(opls_obj)$abs), ncomp = 5)
perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 10) 
plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
MyResult.splsda <-splsda(t(assay(opls_obj)), factor(colData(opls_obj)$abs))
plotIndiv(MyResult.splsda)
plotVar(MyResult.splsda, var.names=TRUE, cutoff = 0.5, cex = 2)
vars <- selectVar(MyResult.splsda)
plotLoadings(MyResult.splsda, contrib = "max", ndisplay = 15, comp = 2)

