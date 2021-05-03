library("Analysis5204")
library("ggplot2")
library("SummarizedExperiment")
library("pcaExplorer")

# Raw data files are loaded automatically when loading the Analysis5204 package.

# First analyse plasma only

# Three proteins (MCP-1, OPG and uPA) are measured in both the Cardiovascular and the Inflammation panel. As seen in the values correlate well between the panels so only the Inflammation data is kept for these proteins.
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
plasma <- SummarizedExperiment(assays = list(npx = npx), colData = metadata)

# correct classes
colData(plasma)$age <- as.numeric(colData(plasma)$age)
colData(plasma)$crp <- as.numeric(colData(plasma)$crp)
colData(plasma)$sr <- as.numeric(colData(plasma)$sr)
colData(plasma)$esr <- as.numeric(colData(plasma)$esr)
colData(plasma)$das28 <- as.numeric(colData(plasma)$das28)
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

# PCA
pcaplot(plasma, 
        intgroup = "group", 
        ellipse = FALSE, 
        text_labels = FALSE, 
        ntop = 50, 
        pcX = 1, 
        pcY = 2)
