#' Simple Restriction
#'
#' This function assists the implementation of a restriction on a covariate in the data
#' table \code{newdf} by setting lines where the covariate is restricted to a user-specified value.
#'
#' @param newdf       Data table containing the simulated data at time \eqn{t}.
#' @param pool        Data table containing the simulated data at times before \eqn{t}.
#' @param restriction List of vectors. Each vector contains as its first entry
#'                    the covariate affected by the restriction; its second entry
#'                    the condition that must be \code{TRUE} for the covariate to be
#'                    modeled; its third entry a function that executes other specific actions
#'                    based on the condition (in this case, this function); and its fourth
#'                    entry some value used by the function (in this case, the value
#'                    the user desires to assign to the covariate when it is not modeled).
#' @param time_name   Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param t           Integer specifying the current time index.
#' @return No value is returned. The data table \code{newdf} is modified in place.
#' @examples
#' ## Estimating the effect of static treatment strategies on risk of a
#' ## failure event with a simple restriction
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
#' # At t0 == 5, assume we have deterministic knowledge that L1 equals 0
#' restrictions <- list(c('L1', 't0 != 5', simple_restriction, 0))
#'
#' gform_basic <- gformula_survival(obs_data = basicdata_nocomp, id = id,
#'                                  time_points = time_points,
#'                                  time_name = time_name, covnames = covnames,
#'                                  outcome_name = outcome_name,
#'                                  covtypes = covtypes,
#'                                  covparams = covparams, ymodel = ymodel,
#'                                  intvars = intvars,
#'                                  interventions = interventions,
#'                                  restrictions = restrictions,
#'                                  int_descript = int_descript,
#'                                  histories = histories, histvars = histvars,
#'                                  basecovs = c('L3'), nsimul = nsimul,
#'                                  seed = 1234)
#' gform_basic
#' }
#'
#' @import data.table
#' @export
simple_restriction <- function(newdf, pool, restriction, time_name, t){
  classtmp <- class(newdf[!eval(parse(text = restriction[[2]]))][[restriction[[1]]]])
  myclass <- paste('as.', classtmp, sep = "")
  newdf[!eval(parse(text = restriction[[2]])), (restriction[[1]]) :=
          get(myclass)(restriction[4])]
}

#' Carry Forward
#'
#' This function assists the implemention of a restriction on a covariate in the date
#' table \code{newdf}. A particular covariate is simulated only when some condition
#' (usually a covariate representing whether a doctor's visit occurred or not) is \code{TRUE}.
#' If the condition is \code{FALSE}, the covariate value is not simulated for that time point
#' and the value is instead carried over from the previous time point.
#'
#' @param newdf       Data table containing the simulated data at time \eqn{t}.
#' @param pool        Data table containing the simulated data at times before \eqn{t}.
#' @param restriction List of vectors. Each vector contains as its first entry
#'                    the covariate affected by the restriction; its second entry
#'                    the condition that must be \code{TRUE} for the covariate to be
#'                    modeled; its third entry a function that executes other specific actions
#'                    based on the condition (in this case, this function); and its fourth
#'                    entry some value used by the function (in this case, this entry is not used).
#' @param time_name   Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param t           Integer specifying the current time index.
#' @return No value is returned. The data table \code{newdf} is modified in place.
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
#' # At t0 == 5, assign L1 its value at the previous time point
#' restrictions <- list(c('L2', 't0 != 5', carry_forward))
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
#'                                  restrictions = restrictions,
#'                                  histories = histories, histvars = histvars,
#'                                  basecovs = c('L3'), nsimul = nsimul,
#'                                  seed = 1234)
#' gform_basic
#' }
#'
#' @import data.table
#' @export
carry_forward <- function(newdf, pool, restriction, time_name, t){
  restrict_ids <- newdf[!eval(parse(text = restriction[[2]]))]$id
  # For restricted individuals, carry covariate value over from prior visit
  classtmp <- class(newdf[newdf$id %in% restrict_ids][[restriction[[1]]]])
  myclass <- paste('as.', classtmp, sep = "")
  newdf[newdf$id %in% restrict_ids, (restriction[[1]]) :=
          get(myclass)(pool[pool[[time_name]] == (t - 1)][newdf$id %in% restrict_ids][[restriction[[1]]]])]
}
