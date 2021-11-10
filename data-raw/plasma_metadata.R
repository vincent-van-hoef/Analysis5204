## code to prepare `plasma_metadata` dataset goes here
library("dplyr")

# Put original data file in data-raw
raw <- "inst/extdata/210429_Vasculitis.xlsx"

# Create a list of the separate "blocks" of metadata in the raw data excel file
if (!file.exists(raw)) {
  stop("Raw data seems to be missing from the extdata folder!")
}

active_disease <- readxl::read_excel(raw,
                          sheet = "Plasma",
                          range = "A3:M72",
                          col_names = TRUE)

healthy_controls <- readxl::read_excel(raw,
                                     sheet = "Plasma",
                                     range = "O3:V141",
                                     col_names = TRUE)

aav_remission <- readxl::read_excel(raw,
                                       sheet = "Plasma",
                                       range = "X3:AG48",
                                       col_names = TRUE)

disease_controls <- readxl::read_excel(raw,
                                     sheet = "Plasma",
                                     range = "AI3:AS50",
                                     col_names = TRUE)
disease_controls <- disease_controls[!disease_controls$ID %in% c("UMU78", "IgA 2156"),]
RA <- split(disease_controls, disease_controls$Disease)$RA
SLE_nephritis <- split(disease_controls, disease_controls$Disease)$`SLE nephritis`

plasmaList <- list("active_disease" = active_disease,
                   "healthy_controls" = healthy_controls,
                   "aav_remission" = aav_remission,
                   "RA" = RA,
                   "SLE_nephritis" = SLE_nephritis)

## Rename all columns
# Create a lookup table for the original column names
old_colnames <- unique(unlist(sapply(plasmaList, function(x) colnames(x))))
new_colnames <- dplyr::recode(old_colnames,
              `Date of sampling` = "sampling_date",
              `Date of birth` = "birth_date",
              `Age at sampling` = "age",
              `Sex (1=male, 0=female, )` = "sex",
              `Diagnosis (1=GPA, 0=MPA)` = "diagnosis",
              `...13` = "cortisone",
              `...11` = "das28",
              `CKD_EPI` = "ckd_epi",
              `CKD-EPI` = "ckd_epi",
              `CKD-EPI real crea/from table` = "ckd_epi",
              `Matched for case:` = "match_case",
              `Matched set number` = "match_set",
              `creatinine umol/l` = "creatinine") %>% tolower()
names(new_colnames) <- old_colnames
# Rename
plasmaList <- lapply(plasmaList, function(x) {colnames(x) <- new_colnames[colnames(x)]; x})
# Add group name
plasmaList <- lapply(names(plasmaList), function(x) {plasmaList[[x]]$group <- x; plasmaList[[x]]})

# Merge list dfs into one df
plasma_metadata <- data.table::rbindlist(plasmaList, fill=TRUE)

# Fix times; conversions to numeric and back to character are necessary
plasma_metadata$sampling_date[grep("^[0-9]{5}$", plasma_metadata$sampling_date)] <- as.character(lubridate::as_date(as.numeric(plasma_metadata$sampling_date[grep("^[0-9]{5}$", plasma_metadata$sampling_date)]), origin = "1899-12-30"))
plasma_metadata$sampling_date[grep("^[0-9]{4}$", plasma_metadata$sampling_date)] <- as.character(lubridate::ymd(as.numeric(plasma_metadata$sampling_date[grep("^[0-9]{4}$", plasma_metadata$sampling_date)]), truncated = 2L))
plasma_metadata$sampling_date <- lubridate::parse_date_time(plasma_metadata$sampling_date, orders = c("ymd"))

plasma_metadata$birth_date[grep("^[0-9]{5}$", plasma_metadata$birth_date)] <- as.character(lubridate::as_date(as.numeric(plasma_metadata$birth_date[grep("^[0-9]{5}$", plasma_metadata$birth_date)]), origin = "1899-12-30"))
plasma_metadata$birth_date[grep("^[0-9]{4}$", plasma_metadata$birth_date)] <- as.character(lubridate::ymd(as.numeric(plasma_metadata$birth_date[grep("^[0-9]{4}$", plasma_metadata$birth_date)]), truncated = 2L))
plasma_metadata$birth_date <- lubridate::parse_date_time(plasma_metadata$birth_date, orders = c("ymd", "y"))

# Clean id
plasma_metadata$id <- tolower(plasma_metadata$id)
plasma_metadata$id <- gsub(" ", "", plasma_metadata$id)
plasma_metadata$id <- gsub("!", "1", plasma_metadata$id)


# Add the BVAS scores
bvas <- "inst/extdata/210603 All active total BVAS score.xlsx"
plasma_bvas <- readxl::read_excel(bvas,
                                  sheet = "Sheet1",
                                  range = "A3:B72",
                                  col_names = TRUE) %>%
  dplyr::rename(bvas = `BVAS score`) %>%
  mutate(ID = tolower(gsub(" ", "", ID))) %>%
  mutate(ID = gsub("vaska122", "vaska122_1", ID),
         ID = gsub("vaska089", "vaska089_1", ID),
         ID = gsub("vaska648", "vaska648_2", ID),
         ID = gsub("vaska114", "vaska114_2", ID),
         ID = gsub("daniel56", "danie56", ID))
plasma_metadata <- plasma_metadata %>% left_join(plasma_bvas, by = c("id" = "ID"))

# Add time in freezer
plasma_metadata$time_in_freezer <- lubridate::time_length(lubridate::as.duration(lubridate::ymd("2020-04-01") - lubridate::ymd(plasma_metadata$sampling_date)), "months")

# Add organ info
organ_info <- "inst/extdata/211109_BVAS_all_active_AAV_patients.xlsx"
plasma_orgs <- readxl::read_excel(organ_info,
                                  sheet = "Plasma",
                                  range = "A4:K65",
                                  col_names = TRUE) %>%
  mutate(ID = tolower(gsub(" ", "", ID))) %>%
  mutate(ID = gsub("vaska089_!", "vaska089_1", ID),
        ID = gsub("daniel56", "danie56", ID)) %>%
  rename(`MUC MEMB /EYES` = "muc_memb_eyes",
         `NERV SYSTEM` = "nerv_system",
         `FINAL SCORE` = "final_score") %>%
  rename_with(tolower)
plasma_metadata <- plasma_metadata %>% left_join(plasma_orgs, by = c("id" = "id"))
  
# Remove no counterpart in data
plasma_metadata <- plasma_metadata[!plasma_metadata$id == "lun1065",]

usethis::use_data(plasma_metadata, overwrite = TRUE)
