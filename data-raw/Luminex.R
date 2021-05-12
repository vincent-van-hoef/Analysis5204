## code to prepare `Luminex` dataset goes here

fileList <- list.files(path = "inst/extdata/", pattern = "^ANCA.*xlsx", full.names = TRUE)
fileNames <- list.files(path = "inst/extdata/", pattern = "^ANCA.*xlsx", full.names = FALSE)
fileNames <- gsub("ANCA_JohannaDahlqvist_", "", fileNames)
fileNames <- gsub(".xlsx", "", fileNames)
conc <- lapply(fileList, function(x) openxlsx::read.xlsx(x, 
                                               sheet = 1, 
                                               rows = c(16:103), 
                                               colNames = TRUE, 
                                               rowNames = FALSE)) 
names(conc) <- fileNames
conc <- lapply(conc, function(x) tidyr::pivot_longer(x, cols = !(Plate:Sample.ID), names_to = "Analyte", values_to = "Conc"))
conc <- purrr::map_df(conc, ~as.data.frame(.x), .id = "id")

lum_dat <- conc %>% 
  dplyr::filter(!Sample.ID %in% c("Blank", "S1", "S2", "S3", "S4", "S5", "S6")) %>%
  dplyr::select(Sample.ID, Analyte, Conc) %>%
  mutate(Analyte = gsub("Complement.Component.C5a|Complement.Componenet.C5a", "C5a", Analyte)) %>%
  tidyr::pivot_wider(id_cols = Sample.ID, names_from = Analyte, values_from = Conc) %>%
  mutate(Sample.ID = tolower(Sample.ID),
         Sample.ID = gsub(" ", "", Sample.ID),
         Sample.ID = gsub("vaska648a_2", "vaska648_2", Sample.ID)) %>%
  dplyr::filter(!grepl('_rep', Sample.ID)) %>%
  dplyr::filter(!Sample.ID %in% c("lu1082", "umu78","lun1059","iga2156"))

lum_dat$`CCL18/PARC` <- as.numeric(lum_dat$`CCL18/PARC`)
lum_dat$`CA15-3/MUC-1` <- as.numeric(lum_dat$`CA15-3/MUC-1`)
lum_dat$`TIMP-1` <- as.numeric(lum_dat$`TIMP-1`)
lum_dat$C5a <- as.numeric(lum_dat$C5a)


usethis::use_data(lum_dat, overwrite = TRUE)
