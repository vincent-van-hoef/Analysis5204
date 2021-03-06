## code to prepare `plasma_npx` dataset goes here

plasma_npx <- OlinkAnalyze::read_NPX("./inst/extdata/Biomarkorer_vid_ANCA_plasma_NPX.xlsx")

# Match id to metdata
plasma_npx$SampleID <- tolower(plasma_npx$SampleID)
plasma_npx$SampleID <- gsub(" ", "", plasma_npx$SampleID)

# Remove double sample
plasma_npx <- plasma_npx[!plasma_npx$SampleID == "lu1082",]
# Remove no counterpart in metadata
plasma_npx <- plasma_npx[!plasma_npx$SampleID %in% c("lun1059", "iga2156", "umu78"),]

usethis::use_data(plasma_npx, overwrite = TRUE)
