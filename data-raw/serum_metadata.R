## code to prepare `serum_metadata` dataset goes here
library("dplyr")

# Put original data file in data-raw
raw <- "inst/extdata/210303_Vasculitis.xlsx"

# Create a list of the separate "blocks" of metadata in the raw data excel file
if (!file.exists(raw)) {
  stop("Raw data seems to be missing from the extdata folder!")
}

active_disease <- readxl::read_excel(raw,
                                     sheet = "Serum",
                                     range = "A3:M30",
                                     col_names = TRUE)

healthy_controls <- readxl::read_excel(raw,
                                       sheet = "Serum",
                                       range = "O3:W57",
                                       col_names = TRUE)

aav_remission <- readxl::read_excel(raw,
                                    sheet = "Serum",
                                    range = "Y3:AH10",
                                    col_names = TRUE)

disease_controls <- readxl::read_excel(raw,
                                       sheet = "Serum",
                                       range = "AJ3:AT10",
                                       col_names = TRUE)
SLE_nephritis <- split(disease_controls, disease_controls$Disease)$SLE

serumList <- list("active_disease" = active_disease,
                   "healthy_controls" = healthy_controls,
                   "aav_remission" = aav_remission,
                   "SLE_nephritis" = SLE_nephritis)

## Rename all columns
# Create a lookup table for the original column names
old_colnames <- unique(unlist(sapply(serumList, function(x) colnames(x))))
new_colnames <- dplyr::recode(old_colnames,
                              `Date of sampling` = "sampling_date",
                              `Date of birth` = "birth_date",
                              `Age at sampling` = "age",
                              `Sex (1=male, 0=female, )` = "sex",
                              `Diagnosis (1=GPA, 0=MPA)` = "diagnosis",
                              `...13` = "cortisone",
                              `...11` = "das28",
                              `Year of sampling` = "sampling_date",
                              `Year of birth` = "birth_date",
                              `Matched for case:` = "match_case",
                              `Matched set number` = "match_set") %>% tolower()
names(new_colnames) <- old_colnames
# Rename
serumList <- lapply(serumList, function(x) {colnames(x) <- new_colnames[colnames(x)]; x})
# Add group name
serumList <- lapply(names(serumList), function(x) {serumList[[x]]$group <- x; serumList[[x]]})

# Merge list dfs into one df; first change dates to characters
serumList <- lapply(serumList, function(x) {x$sampling_date <- as.character(x$sampling_date); x})
serumList <- lapply(serumList, function(x) {x$birth_date <- as.character(x$birth_date); x})
serum_metadata <- data.table::rbindlist(serumList, fill=TRUE)

# Fix times; conversions to numeric and back to character are necessary
serum_metadata$sampling_date[grep("^[0-9]{5}$", serum_metadata$sampling_date)] <- as.character(lubridate::as_date(as.numeric(serum_metadata$sampling_date[grep("^[0-9]{5}$", serum_metadata$sampling_date)]), origin = "1899-12-30"))
serum_metadata$sampling_date[grep("^[0-9]{4}$", serum_metadata$sampling_date)] <- as.character(lubridate::ymd(as.numeric(serum_metadata$sampling_date[grep("^[0-9]{4}$", serum_metadata$sampling_date)]), truncated = 2L))
serum_metadata$sampling_date <- lubridate::parse_date_time(serum_metadata$sampling_date, orders = c("ymd"))

serum_metadata$birth_date[grep("^[0-9]{5}$", serum_metadata$birth_date)] <- as.character(lubridate::as_date(as.numeric(serum_metadata$birth_date[grep("^[0-9]{5}$", serum_metadata$birth_date)]), origin = "1899-12-30"))
serum_metadata$birth_date[grep("^[0-9]{4}$", serum_metadata$birth_date)] <- as.character(lubridate::ymd(as.numeric(serum_metadata$birth_date[grep("^[0-9]{4}$", serum_metadata$birth_date)]), truncated = 2L))
serum_metadata$birth_date <- lubridate::parse_date_time(serum_metadata$birth_date, orders = c("ymd", "y"))

# Clean id
serum_metadata$id <- tolower(serum_metadata$id)
serum_metadata$id <- gsub(" ", "", serum_metadata$id)

# Remove no counterpart in data
serum_metadata <- serum_metadata[!serum_metadata$id %in% c("lin113", "sle_121/a", "sle_156/a"),]

usethis::use_data(serum_metadata, overwrite = TRUE)

