library("Analysis5204")
library("ggplot2")

# Raw data files are loaded automatically when loading the Analysis5204 package.

# First analyse plasma only

# Create a Results directory in the top directory 
plasma_qc_dir <- "Results/Plasma/QC/"
dir.create(plasma_qc_dir, recursive = TRUE)

olink_dist(plasma_npx, "Olink CARDIOVASCULAR III", paste0(plasma_qc_dir, "plasma_CVIII_dist.png"))
olink_dist(plasma_npx, "Olink INFLAMMATION", paste0(plasma_qc_dir, "plasma_INF_dist.png"))

plasma_qc <- OlinkAnalyze::olink_qc_plot(plasma_npx)
ggsave(paste0(plasma_qc_dir, "plasma_qc.png"), plot = plasma_qc)

# Three proteins (MCP-1, OPG and uPA) are measured in both the Cardiovascular and the Inflammation panel. As seen in the values correlate well between the panels so only the Inflammation data is kept for these proteins.
doubles <- c("MCP-1", "OPG", "uPA")
plasma_npx <- plasma_npx[!(plasma_npx$Assay %in% doubles & plasma_npx$Panel == "Olink CARDIOVASCULAR III"),]

#One sample sample (VASKA637) is filtered out due to technical concerns (failed sample according to Olink QC rapport) and another one due to relatively high level of NA values (ra1622.
plasma_npx <- plasma_npx[!plasma_npx$SampleID %in% c("vaska637", "ra1622"),]

# Impute 1 NA assays - better than throwing all sample away
plasma_npx <- plasma_npx %>%
  group_by(Assay) %>%
  mutate(NPX = ifelse(is.na(NPX), median(NPX, na.rm=TRUE), NPX)) %>%
  ungroup(Assay)

