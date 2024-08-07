% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interventions.R
\name{static}
\alias{static}
\title{Static Intervention}
\usage{
static(newdf, pool, intvar, intvals, time_name, t)
}
\arguments{
\item{newdf}{Data table containing the simulated data at time \eqn{t}.}

\item{pool}{Data table containing the simulated data at times before \eqn{t}.}

\item{intvar}{Character string specifying the name of the variable to be intervened on in each round of the simulation.}

\item{intvals}{A list of length 1. The entry is the value of static treatment to be assigned to \code{intvar}.}

\item{time_name}{Character string specifying the name of the time variable in \code{pool} and \code{newdf}.}

\item{t}{Integer specifying the current time index.}
}
\value{
No value is returned. The data table \code{newdf} is modified in place.
}
\description{
This function implements a static intervention (i.e., either constant treatment or no treatment over
all time points) for the specified intervention variable in the data table \code{newdf}.
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
gform_basic
}

}
