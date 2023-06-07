#' Simulated diagnosis data
#'
#' Diagnosis data simulated using a sinusoidal curve
#'
#' @format A data frame with columns:
#' \describe{
#'  \item{ID}{ID of the individual}
#'  \item{ICD}{ICD code of the disease endpoint}
#'  \item{EVENT_DATE}{Date of diagnosis}
#' }
"diagnosis_data"


#' ICD codes for seasonality replication study
#'
#' ICD code table with codes used for the seasonality replication study
#'
#' @format A data frame with columns:
#' \describe{
#'  \item{disease_endpoint}{ID of the disease endpoint}
#'  \item{ICD_version}{Version of ICD}
#'  \item{ICD}{ICD code}
#' }
"ICD_codes"
