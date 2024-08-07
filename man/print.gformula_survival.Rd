% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3methods.R
\name{print.gformula_survival}
\alias{print.gformula_survival}
\alias{print.gformula_continuous_eof}
\alias{print.gformula_binary_eof}
\alias{summary.gformula}
\alias{print.summary.gformula}
\title{Print and summary methods for "gformula" objects}
\usage{
\method{print}{gformula_survival}(
  x,
  all_times = FALSE,
  coefficients = FALSE,
  stderrs = FALSE,
  rmses = FALSE,
  hazardratio = FALSE,
  fits = FALSE,
  ...
)

\method{print}{gformula_continuous_eof}(
  x,
  coefficients = FALSE,
  stderrs = FALSE,
  rmses = FALSE,
  fits = FALSE,
  ...
)

\method{print}{gformula_binary_eof}(
  x,
  coefficients = FALSE,
  stderrs = FALSE,
  rmses = FALSE,
  fits = FALSE,
  ...
)

\method{summary}{gformula}(object, ...)

\method{print}{summary.gformula}(
  x,
  all_times = TRUE,
  coefficients = FALSE,
  stderrs = FALSE,
  rmses = FALSE,
  hazardratio = FALSE,
  fits = TRUE,
  ...
)
}
\arguments{
\item{x}{Object of class "gformula_survival", "gformula_continuous_eof", "gformula_binary_eof", or "summary.gformula" (for \code{print}).}

\item{all_times}{Logical scalar indicating whether to print the results for all time points. This argument is only applicable to objects of class "gformula_survival". If this argument is set to \code{FALSE}, the results are only shown for the final time point. The default is \code{FALSE} for \code{print} and \code{TRUE} for \code{summary}.}

\item{coefficients}{Logical scalar indicating whether to print the model coefficients. The default is \code{FALSE}.}

\item{stderrs}{Logical scalar indicating whether to print the standard error of the model coefficients. The default is \code{FALSE}.}

\item{rmses}{Logical scalar indicating whether to print the model root mean square errors (RMSEs). The default is \code{FALSE}.}

\item{hazardratio}{Logical scalar indicating whether to print the hazard ratio between two interventions (if computed). If bootstrapping was used, 95\% confidence intervals will be given. This argument is only applicable to objects of class "gformula_survival". The default is \code{FALSE}.}

\item{fits}{Logical scalar indicating whether to print summaries of the fitted models for the time-varying covariates, outcome, and competing event (if applicable). This argument is only effective if the argument \code{model_fits} was set to \code{TRUE} in \code{gformula}. The default is \code{FALSE} for \code{print} and \code{TRUE} for \code{summary}.}

\item{...}{Other arguments.}

\item{object}{Object of class "gformula" (for \code{summary}).}
}
\value{
No value is returned for the \code{print} functions. The \code{summary} function returns the object passed to it and adds the class "summary.gformula" to it.
}
\description{
Print and summary method for objects of class "gformula_survival", "gformula_continuous_eof", or "gformula_binary_eof".
}
\examples{
## Estimating the effect of static treatment strategies on risk of a
## failure event
\donttest{
id <- 'id'
time_points <- 7
time_name <- 't0'
covnames <- c('L1', 'L2', 'A')
outcome_name <- 'Y'
outcome_type <- 'survival'
covtypes <- c('binary', 'bounded normal', 'binary')
histories <- c(lagged, lagavg)
histvars <- list(c('A', 'L1', 'L2'), c('L1', 'L2'))
covparams <- list(covmodels = c(L1 ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
                                  L3 + t0,
                                L2 ~ lag1_A + L1 + lag_cumavg1_L1 +
                                  lag_cumavg1_L2 + L3 + t0,
                                A ~ lag1_A + L1 + L2 + lag_cumavg1_L1 +
                                  lag_cumavg1_L2 + L3 + t0))
ymodel <- Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0
intervention1.A <- list(static, rep(0, time_points))
intervention2.A <- list(static, rep(1, time_points))
int_descript <- c('Never treat', 'Always treat')
nsimul <- 10000

gform_basic <- gformula(obs_data = basicdata_nocomp, id = id,
                        time_points = time_points,
                        time_name = time_name, covnames = covnames,
                        outcome_name = outcome_name,
                        outcome_type = outcome_type, covtypes = covtypes,
                        covparams = covparams, ymodel = ymodel,
                        intervention1.A = intervention1.A,
                        intervention2.A = intervention2.A,
                        int_descript = int_descript,
                        histories = histories, histvars = histvars,
                        basecovs = c('L3'), nsimul = nsimul,
                        seed = 1234)
summary(gform_basic)
}

}
\seealso{
\code{\link{gformula}}
}
