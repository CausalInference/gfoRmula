#' Example Dataset for a Survival Outcome without Censoring
#'
#' A dataset consisting of 13,170 observations on 2,500 individuals over 7 time points. Each row in the dataset corresponds to the record of one individual at one time point.
#'
#' @docType data
#'
#' @format A data table with 13,170 rows and 7 variables:
#' \describe{
#'   \item{t0}{Time index.}
#'   \item{id}{Unique identifier for each individual.}
#'   \item{L1}{Binary covariate.}
#'   \item{L2}{Continuous covariate.}
#'   \item{L3}{Continuous baseline covariate. For each individual, the baseline values are repeated at each time point.}
#'   \item{A}{Binary treatment variable.}
#'   \item{Y}{Outcome of interest; time-varying indicator of failure.}
#' }
"basicdata_nocomp"

#' Example Dataset for a Survival Outcome with Censoring
#'
#' A dataset consisting of 11,332 observations on 2,500 individuals over 7 time points. Each row in the dataset corresponds to the record of one individual at one time point. Individuals who are censored at time \eqn{k+1} only have a total of \eqn{k+1} records, which correspond to time indices \eqn{0,..., k}.
#'
#' @docType data
#'
#' @format A data table with 11,332 rows and 8 variables:
#' \describe{
#'   \item{t0}{Time index.}
#'   \item{id}{Unique identifier for each individual.}
#'   \item{L1}{Binary time-varying covariate.}
#'   \item{L2}{Continuous time-varying covariate.}
#'   \item{L3}{Continuous baseline covariate. For each individual, the baseline values are repeated at each time point.}
#'   \item{A}{Binary treatment variable.}
#'   \item{D}{Competing event; time-varying indicator of failure.}
#'   \item{Y}{Outcome of interest; time-varying indicator of failure.}
#' }
"basicdata"


#' Example Dataset for a Binary Outcome at End of Follow-Up
#'
#' A dataset consisting of 17,500 observations on 2,500 individuals over 7 time points. Each row in the dataset corresponds to the record of one individual at one time point.
#'
#' @docType data
#'
#' @format A data table with 17,500 rows and 7 variables:
#' \describe{
#'   \item{time}{Time index.}
#'   \item{id_num}{Unique identifier for each individual.}
#'   \item{cov1}{Binary time-varying covariate.}
#'   \item{cov2}{Continuous time-varying covariate.}
#'   \item{cov3}{Continuous baseline covariate. For each individual, the baseline values are repeated at each time point.}
#'   \item{treat}{Binary treatment variable.}
#'   \item{outcome}{Binary outcome of interest. Because this outcome is only defined at the end of follow-up, values of \code{NA} are given in all other time points.}
#' }
"binary_eofdata"


#' Example Dataset for a Continuous Outcome at End of Follow-Up
#'
#' A dataset consisting of 17,500 observations on 2,500 individuals over 7 time points. Each row in the dataset corresponds to the record of one individual at one time point.
#'
#' @docType data
#'
#' @format A data table with 17,500 rows and 7 variables:
#' \describe{
#'   \item{t0}{Time index.}
#'   \item{id}{Unique identifier for each individual.}
#'   \item{L1}{Categorical time-varying covariate.}
#'   \item{L2}{Continuous time-varying covariate.}
#'   \item{L3}{Continuous baseline covariate. For each individual, the baseline values are repeated at each time point.}
#'   \item{A}{Binary treatment variable.}
#'   \item{Y}{Continuous outcome of interest. Because this outcome is only defined at the end of follow-up, values of \code{NA} are given in all other time points.}
#' }
"continuous_eofdata"
