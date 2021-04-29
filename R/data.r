#' Clinical data of all the plasma samples
#'
#' A dataset containing the metadata and clinical data of the plasma samples.
#' All samples have been put together in one data frame and the dates have been made uniform.
#' Sample id have been lower cased and spaces removed. One sample did not have a counterpart in the data and was removed ("lun1065")
#'
#'
#' @format A data frame with 298 rows and 19 variables:
#' \describe{
#'   \item{id}{ID code; unique for each sample}
#'   \item{sampling_date}{date of sampling}
#'   \item{birth_date}{date of birth}
#'   \item{age}{age of patient at sampling date}
#'   \item{sex}{gender; 1 for male, 0 for female}
#'   \item{diagnosis}{diagnosis; 1 for GPA, 0 for MPA}
#'   \item{pr3-anca}{Presence of pr3 auto-antibody}
#'   \item{mpo-anca}{Presence of mpo auto-antibody}
#'   \item{group}{Group membership}
#'   \item{creatinine}{LEvels of creatinine}
#'   \item{ckd_epi}{Measure for kidney activity}
#'   \item{esr}{unknown}
#'   \item{match_case}{matching case}
#'   \item{match_set}{matching set}
#'   \item{sr}{unknown}
#'   \item{das28}{unknown}
#'   \item{cortisone}{unknown}
#'   \item{crp}{unknown}
#'   \item{disease}{disease type}
#' }
#' @source Delivered by Johanna
"plasma_metadata"

#' Clinical data of all the serum samples
#'
#' A dataset containing the metadata and clinical data of the serum samples.
#' All samples have been put together in one data frame and the dates have been made uniform.
#' Sample id have been lower cased and spaces removed. Three samples did not have a counterpart in the data and were removed.
#'
#' @format A data frame with 92 rows and 19 variables:
#' \describe{
#'   \item{id}{ID code; unique for each sample}
#'   \item{sampling_date}{date of sampling}
#'   \item{birth_date}{date of birth}
#'   \item{age}{age of patient at sampling date}
#'   \item{sex}{gender; 1 for male, 0 for female}
#'   \item{diagnosis}{diagnosis; 1 for GPA, 0 for MPA}
#'   \item{pr3-anca}{Presence of pr3 auto-antibody}
#'   \item{mpo-anca}{Presence of mpo auto-antibody}
#'   \item{group}{Group membership}
#'   \item{creatinine}{Levels of creatinine}
#'   \item{ckd-epi}{Measure for kidney activity}
#'   \item{esr}{unknown}
#'   \item{match_case}{matching case}
#'   \item{match_set}{matching set}
#'   \item{sr}{unknown}
#'   \item{das28}{unknown}
#'   \item{cortisone}{unknown}
#'   \item{crp}{unknown}
#'   \item{disease}{disease type}
#' }
#' @source Delivered by Johanna
"serum_metadata"

#' @title plasma_npx
#' @description NPX data of the plasma samples
#' @format A data frame with 54832 rows and 13 variables:
#' \describe{
#'   \item{\code{SampleID}}{character Sample ID; matches metadata of plasma}
#'   \item{\code{Index}}{integer Sample numbering}
#'   \item{\code{OlinkID}}{character OlinkID}
#'   \item{\code{UniProt}}{character UniProc code for the protein}
#'   \item{\code{Assay}}{character Protein Symbol}
#'   \item{\code{MissingFreq}}{character Frequency of missingness}
#'   \item{\code{Panel}}{character Panel Name}
#'   \item{\code{Panel_Version}}{character Paner Version}
#'   \item{\code{PlateID}}{character Plate ID}
#'   \item{\code{QC_Warning}}{character Is there a QC waringn?}
#'   \item{\code{LOD}}{double Limit of Detection}
#'   \item{\code{NPX}}{double NPX value; abundance measurement}
#'   \item{\code{Normalization}}{character Normalization method}
#'}
#' @details DETAILS
"plasma_npx"

#' @title serum_npx
#' @description NPX data of the serum samples
#' @format A data frame with 16928 rows and 13 variables:
#' \describe{
#'   \item{\code{SampleID}}{character Sample ID; matches metadata of serum}
#'   \item{\code{Index}}{integer Sample numbering}
#'   \item{\code{OlinkID}}{character OlinkID}
#'   \item{\code{UniProt}}{character UniProc code for the protein}
#'   \item{\code{Assay}}{character Protein Symbol}
#'   \item{\code{MissingFreq}}{character Frequency of missingness}
#'   \item{\code{Panel}}{character Panel Name}
#'   \item{\code{Panel_Version}}{character Paner Version}
#'   \item{\code{PlateID}}{character Plate ID}
#'   \item{\code{QC_Warning}}{character Is there a QC waringn?}
#'   \item{\code{LOD}}{double Limit of Detection}
#'   \item{\code{NPX}}{double NPX value; abundance measurement}
#'   \item{\code{Normalization}}{character Normalization method}
#'}
#' @details DETAILS
"serum_npx"
