#' Clinical data of all the plasma samples
#'
#' A dataset containing the metadata and clinical data of the plasma samples.
#' All samples have been put together in one data frame and the dates have been made uniform.
#'
#' @format A data frame with 299 rows and 19 variables:
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
#'
#' @format A data frame with 95 rows and 19 variables:
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
