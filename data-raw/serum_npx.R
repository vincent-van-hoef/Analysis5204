## code to prepare `serum_npx` dataset goes here

serum_npx <- OlinkAnalyze::read_NPX("inst/extdata/Biomarkörer_vid_ANCA_serum_NPX.xlsx")

# Match id to metdata
serum_npx$SampleID <- tolower(serum_npx$SampleID)
serum_npx$SampleID <- gsub(" ", "", serum_npx$SampleID)
serum_npx$SampleID <- gsub("\\(u81\\)", "", serum_npx$SampleID)
serum_npx$SampleID <- gsub("\\(u82\\)", "", serum_npx$SampleID)

usethis::use_data(serum_npx, overwrite = TRUE)
