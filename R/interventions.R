# Copyright (c) 2019 The President and Fellows of Harvard College

#' Example Intervention 1
#'
#' This function implements an example intervention for the specified intervention variable in
#' the data table \code{newdf}. It allows the user to define a cutoff for the time variable, where a
#' different intervention is assigned depending on whether the time is before or after the cutoff.
#' The data table \code{newdf} is modified in place.
#'
#' @param time_name  Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf      Data table containing the simulated data at time \eqn{t}.
#' @param pool       Data table containing the simulated data at times before \eqn{t}.
#' @param intvar     Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param intvals    A list. The first entry is user-specified cutoff for the time variable, above
#'                   which \code{intvar} will be assigned some user-specified
#'                   value and below which the intervention variable will be assigned some
#'                   other user-specified value. The second and third entries are the values
#'                   assigned to \code{intvar} if the time value is above and below, respectively, the cutoff.
#' @param t          Integer specifying the current time index.
#' @import data.table
#' @export
example_intervention_1 <- function(newdf, pool, intvar, intvals, time_name, t){
  newdf[newdf[[time_name]] < intvals[[1]], (intvar) := intvals[[2]]]
  newdf[newdf[[time_name]] >= intvals[[1]], (intvar) := intvals[[3]]]
}

#' Example Intervention 2
#'
#' This function implements an example intervention for the specified intervention variable in
#' the data table \code{newdf}. It allows the user to override the natural course, such that if the
#' simulated natural course treatment is equal to some user-defined value, the treatment
#' variable is set to 1, and otherwise the treatment variable is set to 0. The data table
#' \code{newdf} is modified in place.
#'
#' @param time_name  Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf      Data table containing the simulated data at time \eqn{t}.
#' @param pool       Data table containing the simulated data at times before \eqn{t}.
#' @param intvar     Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param intvals    A list of length 1. The entry in the list is a user-specified value for the variable \code{A}, such that
#'                   if the simulated value of \code{A} is equal to some value, \code{intvar} is set equal
#'                   to 1, while otherwise the intervention variable is set equal to 0.
#' @param t          Integer specifying the current time index.
#' @import data.table
#' @export
example_intervention_2 <- function(newdf, pool, intvar, intvals, time_name, t){
  newdf[newdf$A == intvals[[1]], (intvar) := 1]
  newdf[newdf$A != intvals[[1]], (intvar) := 0]
}

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
#' @export
natural <- function(newdf, pool, intvar, intvals, time_name, t){
}

#' Static Intervention
#'
#' This function implements a static intervention (i.e., either constant treatment or no treatment over
#' all time points) for the specified intervention variable in the data table \code{newdf}. The data table \code{newdf} is modified in place.
#'
#' @param time_name  Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf      Data table containing the simulated data at time \eqn{t}.
#' @param pool       Data table containing the simulated data at times before \eqn{t}.
#' @param intvar     Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param intvals    A list of length 1. The entry is the value of static treatment to be assigned to \code{intvar}.
#' @param t          Integer specifying the current time index.
#' @import data.table
#' @export
static <- function(newdf, pool, intvar, intvals, time_name, t){
  set(newdf, j = intvar, value = intvals[[t+1]])
}

#' Threshold Intervention
#'
#' This function implements a threshold intervention (i.e., once treatment bypasses a certain
#' threshold, it remains at that threshold until end of follow-up) for the specified
#' intervention variable in the data table \code{newdf}. The data table \code{newdf} is modified in place.
#'
#' @param time_name  Character string specifying the name of the time variable in \code{pool} and \code{newdf}.
#' @param newdf      Data table containing the simulated data at time \eqn{t}.
#' @param pool       Data table containing the simulated data at times before \eqn{t}.
#' @param intvar     Character string specifying the name of the variable to be intervened on in each round of the simulation.
#' @param intvals    A list of length 2. The first entry is lower bound of the threshold, and the second entry is the upper bound.
#' @param t          Integer specifying the current time index.
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
#' in the data table \code{newdf}. The data table \code{newdf} is modified in place.
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
intfunc <- function(newdf, pool, intervention, intvar, int_time, time_name, t){
  if (t %in% int_time){
    lapply(1:length(intervention), FUN = function(i) {
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
