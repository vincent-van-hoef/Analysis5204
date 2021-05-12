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
  dplyr::filter(!Sample.ID %in% c("lu1082","vaska637","umu78","lun1059","iga2156","ra1622"))


usethis::use_data(lum_dat, overwrite = TRUE)
