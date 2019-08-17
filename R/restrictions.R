#' Simple Restriction
#'
#' This function assists the implementation of a restriction on a covariate in the data
#' table \code{newdf} by setting lines where the covariate is restricted to a user-specified value.
#' The data table \code{newdf} is modified in place.
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
#' @import data.table
#' @export
simple_restriction <- function(newdf, pool, restriction, time_name, t){
  classtmp <- class(newdf[!eval(parse(text = restriction[[2]]))][[restriction[[1]]]])
  myclass <- paste('as.', classtmp, sep = "")
  newdf[!eval(parse(text = restriction[[2]])), (restriction[[1]]) :=
          get(myclass)(restriction[4])]
}

#' Visit Process
#'
#' This function assists the implemention of a visit process, where a particular covariate is simulated
#' only when some condition (usually a covariate representing whether a doctor's visit
#' occurred or not) is \code{TRUE}. If the condition is \code{FALSE}, the covariate value is not
#' simulated for that time point and the value is instead carried over from the previous
#' time point. If the condition is \code{FALSE} for more than a certain number of times
#' consecutively (i.e., if more than a certain number of visits are missed in a row),
#' the individual is censored.
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
