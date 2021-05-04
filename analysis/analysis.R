# Clear global env from earlier session
rm(list=ls())

library("Analysis5204")
data(list=c("plasma_metadata", "plasma_npx", "serum_metadata", "serum_npx"), package = "Analysis5204")
library("ggplot2")
library("SummarizedExperiment")
library("pcaExplorer")
library("emmeans")

# Raw data files are loaded automatically when loading the Analysis5204 package.

# First analyse plasma only

# Three proteins (MCP-1, OPG and uPA) are measured in both the Cardiovascular and the Inflammation panel. As seen in the values correlate well between the panels so only the Inflammation data is kept for these proteins. Has to been done before making the summarizedExperiment because pivot_wider otherwise fails.
doubles <- c("MCP-1", "OPG", "uPA")
plasma_npx <- plasma_npx[!(plasma_npx$Assay %in% doubles & plasma_npx$Panel == "Olink CARDIOVASCULAR III"),]

# Make a data=metadata object
npx <- plasma_npx %>% 
  dplyr::select(SampleID, Assay, NPX) %>% 
  tidyr::pivot_wider(names_from = Assay, values_from = NPX) %>%
  tibble::column_to_rownames("SampleID") %>%
  t()
metadata <- plasma_metadata %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("id")
# Make sure order of assays columns and colData ros is the same!!
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

olink_dist(plasma_npx, "Olink CARDIOVASCULAR III", paste0(plasma_qc_dir, "plasma_CVIII_dist.png"))
olink_dist(plasma_npx, "Olink INFLAMMATION", paste0(plasma_qc_dir, "plasma_INF_dist.png"))

plasma_qc <- OlinkAnalyze::olink_qc_plot(plasma_npx)
ggsave(paste0(plasma_qc_dir, "plasma_qc.png"), plot = plasma_qc)

#One sample sample (VASKA637) is filtered out due to technical concerns (failed sample according to Olink QC rapport) and another one due to relatively high level of NA values (ra1622.
plasma <- plasma[,!colnames(plasma) %in% c("vaska637", "ra1622")]

# Impute 1 NA assays - better than throwing all sample away
assay(plasma) <-  apply(assay(plasma), 2, function(x) ifelse(is.na(x), median(x, na.rm=T), x))

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

# Anova
# Create a Results directory in the top directory 
plasma_uni_dir <- "Results/Plasma/Univariate/"
dir.create(plasma_uni_dir, recursive = TRUE)
#Create olink-function compatible dataset; reduced to the summarizedExperiment samples
obj <- plasma_npx %>% 
  filter(SampleID %in% colnames(plasma)) %>% 
  left_join(tibble::rownames_to_column(as.data.frame(colData(plasma))), by = c("SampleID" = "rowname"))

ph_anova <- OlinkAnalyze::olink_anova_posthoc(df = obj, 
                                              variable = "group", 
                                              covariates = c("age", "ckd_epi"),
                                              effect = "group",
                                              verbose = FALSE)

for(comp in unique(ph_anova$contrast)){
  
compname <- gsub(" - ", "_vs_", comp)
dir.create(paste0(plasma_uni_dir, compname, "/"), recursive = TRUE)

df_sub <- ph_anova %>% filter(contrast == comp)
p1 <- ggplot(df_sub, aes(x=estimate, y = -log10(Adjusted_pval))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) + 
  geom_vline(xintercept = 0) +
  #xlim(c(-3,3)) +
  theme_bw()
ggsave(paste0(plasma_uni_dir, compname, "/", compname, "_volcano_anova_posthoc.pdf"), plot = p1)

}
