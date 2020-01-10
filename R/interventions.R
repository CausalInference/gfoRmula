#' Natural Course Intervention
#'
#' This function implements the natural course (i.e., model-based simulation) for the specified
#' intervention variable in the data table newdf. Because \code{newdf} is initiated with the
#' natural course, this function does nothing.
#'
#' @param time_name  Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf      Data table containing the simulated data at time \eqn{t}.
#' @param pool       Data table containing the simulated data at times before \eqn{t}.
#' @param intvar     Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param intvals    A list of user-specified values for the intervention.
#' @param t          Integer specifying the current time index.
#' @return No value is returned.
#' @keywords internal
natural <- function(newdf, pool, intvar, intvals, time_name, t){
}

#' Static Intervention
#'
#' This function implements a static intervention (i.e., either constant treatment or no treatment over
#' all time points) for the specified intervention variable in the data table \code{newdf}.
#'
#' @param time_name  Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf      Data table containing the simulated data at time \eqn{t}.
#' @param pool       Data table containing the simulated data at times before \eqn{t}.
#' @param intvar     Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param intvals    A list of length 1. The entry is the value of static treatment to be assigned to \code{intvar}.
#' @param t          Integer specifying the current time index.
#' @return No value is returned. The data table \code{newdf} is modified in place.
#' @import data.table
#'
#' @examples
#' ## Estimating the effect of static treatment strategies on risk of a
#' ## failure event
#' \donttest{
#' id <- 'id'
#' time_points <- 7
#' time_name <- 't0'
#' covnames <- c('L1', 'L2', 'A')
#' outcome_name <- 'Y'
#' covtypes <- c('binary', 'bounded normal', 'binary')
#' histories <- c(lagged, lagavg)
#' histvars <- list(c('A', 'L1', 'L2'), c('L1', 'L2'))
#' covparams <- list(covmodels = c(L1 ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
#'                                   L3 + t0,
#'                                 L2 ~ lag1_A + L1 + lag_cumavg1_L1 +
#'                                   lag_cumavg1_L2 + L3 + t0,
#'                                 A ~ lag1_A + L1 + L2 + lag_cumavg1_L1 +
#'                                   lag_cumavg1_L2 + L3 + t0))
#' ymodel <- Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0
#' intvars <- list('A', 'A')
#' interventions <- list(list(c(static, rep(0, time_points))),
#'                       list(c(static, rep(1, time_points))))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_basic <- gformula_survival(obs_data = basicdata_nocomp, id = id,
#'                                  time_points = time_points,
#'                                  time_name = time_name, covnames = covnames,
#'                                  outcome_name = outcome_name,
#'                                  covtypes = covtypes,
#'                                  covparams = covparams, ymodel = ymodel,
#'                                  intvars = intvars,
#'                                  interventions = interventions,
#'                                  int_descript = int_descript,
#'                                  histories = histories, histvars = histvars,
#'                                  basecovs = c('L3'), nsimul = nsimul,
#'                                  seed = 1234)
#' gform_basic
#' }
#'
#' @export
static <- function(newdf, pool, intvar, intvals, time_name, t){
  set(newdf, j = intvar, value = intvals[[t+1]])
}

#' Threshold Intervention
#'
#' This function implements a threshold intervention (i.e., once treatment bypasses a certain
#' threshold, it remains at that threshold until end of follow-up) for the specified
#' intervention variable in the data table \code{newdf}.
#'
#' @param time_name  Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf      Data table containing the simulated data at time \eqn{t}.
#' @param pool       Data table containing the simulated data at times before \eqn{t}.
#' @param intvar     Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param intvals    A list of length 2. The first entry is lower bound of the threshold, and the second entry is the upper bound.
#' @param t          Integer specifying the current time index.
#' @return No value is returned. The data table \code{newdf} is modified in place.
#'
#' @examples
#'
#' ## Estimating the effect of threshold interventions on the mean of a binary
#' ## end of follow-up outcome
#' \donttest{
#' id <- 'id_num'
#' time_name <- 'time'
#' covnames <- c('cov1', 'cov2', 'treat')
#' outcome_name <- 'outcome'
#' histories <- c(lagged, cumavg)
#' histvars <- list(c('treat', 'cov1', 'cov2'), c('cov1', 'cov2'))
#' covtypes <- c('binary', 'zero-inflated normal', 'normal')
#' covparams <- list(covmodels = c(cov1 ~ lag1_treat + lag1_cov1 + lag1_cov2 + cov3 +
#'                                   time,
#'                                 cov2 ~ lag1_treat + cov1 + lag1_cov1 + lag1_cov2 +
#'                                   cov3 + time,
#'                                 treat ~ lag1_treat + cumavg_cov1 +
#'                                   cumavg_cov2 + cov3 + time))
#' ymodel <- outcome ~  treat + cov1 + cov2 + lag1_cov1 + lag1_cov2 + cov3
#' intvars <- list('treat', 'treat')
#' interventions <- list(list(c(static, rep(0, 7))),
#'                       list(c(threshold, 1, Inf)))
#' int_descript <- c('Never treat', 'Threshold - lower bound 1')
#' nsimul <- 10000
#' ncores <- 2
#'
#' gform_bin_eof <- gformula_binary_eof(obs_data = binary_eofdata, id = id,
#'                                      time_name = time_name,
#'                                      covnames = covnames,
#'                                      outcome_name = outcome_name,
#'                                      covtypes = covtypes,
#'                                      covparams = covparams,
#'                                      ymodel = ymodel,
#'                                      intvars = intvars,
#'                                      interventions = interventions,
#'                                      int_descript = int_descript,
#'                                      histories = histories, histvars = histvars,
#'                                      basecovs = c("cov3"), seed = 1234,
#'                                      parallel = TRUE, nsamples = 5,
#'                                      nsimul = nsimul, ncores = ncores)
#' gform_bin_eof
#' }
#'
#' @import data.table
#' @export
threshold <- function(newdf, pool, intvar, intvals, time_name, t){
  if (nrow(newdf[newdf[[intvar]] < intvals[[1]]]) != 0){
    classtmp <- class(newdf[[intvar]])
    myclass <- paste('as.', classtmp, sep = "")
    newdf[newdf[[intvar]] < intvals[[1]], (intvar) := get(myclass)(intvals[[1]])]
  }
  if (nrow(newdf[newdf[[intvar]] > intvals[[2]]]) != 0){
    classtmp <- class(newdf[[intvar]])
    myclass <- paste('as.', classtmp, sep = "")
    newdf[newdf[[intvar]] > intvals[[2]], (intvar) := get(myclass)(intvals[[2]])]
  }
}

#' Execute Intervention
#'
#' This internal function executes the intervention of interest on the specified intervention variable
#' in the data table \code{newdf}.
#'
#' @param time_name     Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf         Data table containing the simulated data at time \eqn{t}.
#' @param pool          Data table containing the simulated data at times before \eqn{t}.
#' @param intervention  List, whose elements are lists of vectors. Each vector contains a function
#'                      implementing a particular intervention on a single variable, optionally
#'                      followed by one or more "intervention values" (i.e.,
#'                      integers used to specify the treatment regime).
#' @param intvar        Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param int_time      Vector specifying the time points in which the intervention is applied.
#' @param t             Integer specifying the current time index.
#' @return No value is returned. The data table \code{newdf} is modified in place.
#' @keywords internal
intfunc <- function(newdf, pool, intervention, intvar, int_time, time_name, t){
  if (t %in% int_time){
    lapply(seq_along(intervention), FUN = function(i) {
      if (length(intervention[[i]]) == 1){ # Check if intervention contains just function
        intervention[[i]][[1]](newdf, pool, intvar[i], intvals = NULL, time_name, t)
      } else {
        # If intervention contains intervention values, pass those values to intervention
        # function
        intervention[[i]][[1]](newdf, pool, intvar[i], intvals = intervention[[i]][2:length(intervention[[i]])],
                               time_name, t)
      }
    })
  }
}
