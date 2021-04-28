## code to prepare `plasma_metadata` dataset goes here
library("dplyr")
library("readxl")

# Create a list of the separate "blocks" of metadata in the raw data excel file
if (!file.exists("data-raw/210303_Vasculitis.xlsx")) {
  stop("Raw data seems to be missing from the data-raw folder!")
}

raw <- 'data-raw/210303_Vasculitis.xlsx'

active_disease <- readxl::read_excel('data-raw/210303_Vasculitis.xlsx',
                          sheet = "Plasma",
                          range = "A3:M72",
                          col_names = TRUE)

healthy_controls <- readxl::read_excel('data-raw/210303_Vasculitis.xlsx',
                                     sheet = "Plasma",
                                     range = "O3:V141",
                                     col_names = TRUE)

aav_remission <- readxl::read_excel('data-raw/210303_Vasculitis.xlsx',
                                       sheet = "Plasma",
                                       range = "X3:AG48",
                                       col_names = TRUE)

disease_controls <- readxl::read_excel('data-raw/210303_Vasculitis.xlsx',
                                     sheet = "Plasma",
                                     range = "AI3:AS50",
                                     col_names = TRUE)


usethis::use_data(plasma_metadata, overwrite = TRUE)
