## code to prepare `plasma_metadata` dataset goes here

# Create a list of the separate "blocks" of metadata in the raw data excel file
if (!file.exists(system.file("extdata", "210303_Vasculitis.xlsx", package = "Analysis5204"))) {
  stop("Raw data seems to be missing from the data-raw folder!")
}

# Put original data file in data-raw
raw <- system.file("extdata", "210303_Vasculitis.xlsx", package = "Analysis5204")

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
RA <- split(disease_controls, disease_controls$Disease)$RA
IgA_nephritis <- split(disease_controls, disease_controls$Disease)$`IgA nephritis`
SLE_nephritis <- split(disease_controls, disease_controls$Disease)$`SLE nephritis`

plasmaList <- list("active_disease" = active_disease,
                   "healthy_controls" = healthy_controls,
                   "aav_remission" = aav_remission,
                   "RA" = RA,
                   "IgA_nephritis" = IgA_nephritis,
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


usethis::use_data(plasma_metadata, overwrite = TRUE)
