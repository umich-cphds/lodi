#' Synthetic toy data for clmi
#'
#' @format A data.frame with 100 observations on 6 variables:
#' \describe{
#'   \item{id}{Patient ID number.}
#'   \item{case_cntrl}{Patient's case-control status. Either 1 or 0.}
#'   \item{poll}{Concentration of pollutant in patient's blood sample.}
#'   \item{smoking}{Smoking status. Either 1 or 0.}
#'   \item{gender}{Gender. 1 for male, 0 for female.}
#'   \item{batch1}{Batch status. Integer}
#'   \item{lod}{batch's limit of detection for patient.}
#' }
"toy_data"
