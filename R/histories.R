#' History functions
#'
#' These functions create new columns in an input data table for covariate histories. Users must specify which covariates are to be used in the history functions.
#'
#' \code{lagged} creates new columns for lagged versions of existing
#' variables in the dataset. The user must specify which variables are to be lagged.
#'
#' \code{cumavg} creates new columns for the cumulative average up until
#' time \eqn{t} of existing variables in the dataset.
#'
#' \code{lagavg} creates new columns for the "lagged cumulative average"
#' (cumulative average up until time t, then lagged by one time unit) up until time \eqn{t} of existing
#' variables in the dataset.
#'
#' @param t          Integer specifying the current time index.
#' @param time_name  Character string specifying the name of the time variable in \code{pool}.
#' @param pool       Data table containing all information prior to time \eqn{t} (\eqn{t} noninclusive).
#' @param histvars   Vector of character strings specifying the names of the variables for which history functions are to be applied.
#' @param histvals   For \code{lagged}, this argument is a vector specifying the lags used in the model statements (e.g., if \code{lag1_varname} and \code{lag2_varname} were included in the model statements, this vector would be \code{c(1,2)}). For \code{lagavg}, this argument is a numeric vector specifying the lag averages used in the model statements.
#' @param id_name    Character string specifying the name of the ID variable in \code{pool}.
#' @param baselags   Logical scalar for specifying the convention used for lagi and lag_cumavgi terms in the model statements when pre-baseline times are not
#'                   included in \code{obs_data} and when the current time index, \eqn{t}, is such that \eqn{t < i}. If this argument is set to \code{FALSE}, the value
#'                   of all lagi and lag_cumavgi terms in this context are set to 0 (for non-categorical covariates) or the reference
#'                   level (for categorical covariates). If this argument is set to \code{TRUE}, the value of lagi and lag_cumavgi terms
#'                   are set to their values at time 0. The default is \code{FALSE}.
#' @param below_zero_indicator Logical scalar indicating whether the observed data set contains rows for time \eqn{t < 0}.
#' @return No value is returned. The data table \code{pool} is modified in place.
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
#' @import data.table
#' @export
lagged <- function(pool, histvars, histvals, time_name, t, id_name, baselags,
                   below_zero_indicator){
  for (i in histvals){
    if (t < i & baselags & !below_zero_indicator){
      current_ids <- unique(pool[pool[[time_name]] == t][[id_name]])
      lapply(histvars, FUN = function(histvar){
        classtmp <- class(pool[pool[[time_name]] == t][[histvar]])
        myclass <- paste('as.', classtmp, sep = "")
        pool[pool[[time_name]] == t, (paste("lag", i, "_", histvar, sep = "")) :=
               get(myclass)(pool[pool[[time_name]] == 0 & get(id_name) %in% current_ids][[histvar]])]
      })
    } else if (t < i & !below_zero_indicator){
      lapply(histvars, FUN = function(histvar){
        classtmp <- class(pool[pool[[time_name]] == t][[histvar]])
        myclass <- paste('as.', classtmp, sep = "")
        if (is.factor(pool[pool[[time_name]] == t][[histvar]])){
          reflevel <- levels(pool[pool[[time_name]] == t][[histvar]])[1]
          pool[pool[[time_name]] == t, (paste("lag", i, "_", histvar, sep = "")) := reflevel]
        } else {
          pool[pool[[time_name]] == t, (paste("lag", i, "_", histvar, sep = "")) := get(myclass)(0)]
        }
      })
    } else {
      current_ids <- unique(pool[pool[[time_name]] == t][[id_name]])
      lapply(histvars, FUN = function(histvar){
        classtmp <- class(pool[pool[[time_name]] == t][[histvar]])
        myclass <- paste('as.', classtmp, sep = "")
        pool[pool[[time_name]] == t, (paste("lag", i, "_", histvar, sep = "")) :=
               get(myclass)(pool[pool[[time_name]] == t - i & get(id_name) %in% current_ids][[histvar]])]
      })
    }
  }
}

#' @rdname lagged
#' @export
cumavg <- function(pool, histvars, time_name, t, id_name, below_zero_indicator){
  if (below_zero_indicator){
    denom <- t + 1 - min(pool[[time_name]])
  } else {
    denom <- t + 1
  }
  if (t == 0 & !below_zero_indicator){
    # At first time point, set all cumulative averages equal to the actual value of the
    # variable
    lapply(histvars, FUN = function(histvar){
      pool[get(time_name) == t, (paste("cumavg_", histvar, sep = "")) :=
             as.double(pool[get(time_name) == t][[histvar]])]
    })
  } else {
    # At subsequent time points, create new column containing calculated cumulative
    # average until that time point
    current_ids <- unique(pool[get(time_name) == t][[id_name]])
    colnam <- colnames(pool)
    if (!(paste("cumavg_", "_", histvars[1], sep = "") %in% colnam)){
      # The cumulative average variable was not created yet
      # Therefore, cannot use recursive formula
      id_factor <- is.factor(pool[[id_name]])
      if (id_factor){
        lapply(histvars, FUN = function(histvar){
          pool[get(time_name) == t, (paste("cumavg_", histvar, sep = "")) :=
                 as.double(tapply(pool[get(id_name) %in% current_ids &
                                         get(time_name) <= t][[histvar]],
                                  droplevels(pool[get(id_name) %in% current_ids &
                                                    get(time_name) <= t][[id_name]]),
                                  FUN = mean))]
        })
      } else {
        lapply(histvars, FUN = function(histvar){
          pool[get(time_name) == t, (paste("cumavg_", histvar, sep = "")) :=
                 as.double(tapply(pool[get(id_name) %in% current_ids &
                                         get(time_name) <= t][[histvar]],
                                  pool[get(id_name) %in% current_ids &
                                         get(time_name) <= t][[id_name]], FUN = mean))]
        })
      }
    } else {
      # The cumulative average variable was already created
      # Therefore, can use recursive formula
      for (histvar in histvars){
        pool[get(time_name) == t, (paste("cumavg_", histvar, sep = "")) :=
               as.double(
                 (pool[get(id_name) %in% current_ids & get(time_name) == (t-1)][[paste("cumavg_", histvar, sep = "")]] * (denom - 1)
                  + pool[get(id_name) %in% current_ids & get(time_name) == t][[histvar]]) / denom
               )
        ]
      }
    }
  }
}

#' @rdname lagged
#' @export
lagavg <- function(pool, histvars, histvals, time_name, t, id_name, baselags,
                   below_zero_indicator){
  for (i in histvals){
    if (t < i & !below_zero_indicator){
      if (baselags){
        current_ids <- unique(pool[pool[[time_name]] == t][[id_name]])
        lapply(histvars, FUN = function(histvar){
          pool[pool[[time_name]] == t, (paste("lag_cumavg", i, "_", histvar, sep = "")) :=
                 as.double(pool[pool[[time_name]] == 0 & get(id_name) %in% current_ids, ][[histvar]])]
        })
      } else {
        # At time point less than i, set all lagged cumulative averages equal to 0
        lapply(histvars, FUN = function(histvar){
          pool[pool[[time_name]] == t, (paste("lag_cumavg", i, "_", histvar, sep = "")) := as.double(0)]

        })
      }

    } else if (t == (i + min(pool[[time_name]]))){
      # At time point i, set all lagged cumulative averages equal to the actual value
      # of the variable
      current_ids <- unique(pool[pool[[time_name]] == t][[id_name]])
      lapply(histvars, FUN = function(histvar){
        pool[pool[[time_name]] == t, (paste("lag_cumavg", i, "_", histvar, sep = "")) :=
               as.double(pool[pool[[time_name]] == t - i & get(id_name) %in% current_ids, ][[histvar]])]
      })
    } else if (t > (i + min(pool[[time_name]]))) {

      # At time points after i, create new column containing calculated lagged cumulative
      # average until that time point
      current_ids <- unique(pool[pool[[time_name]] == t][[id_name]])
      colnam <- colnames(pool)
      if (!(paste("lag_cumavg", i, "_", histvars[1], sep = "") %in% colnam)){
        # The lagged cumulative average variable was not created yet
        # Therefore, cannot use recursive formula
        id_factor <- is.factor(pool[[id_name]])
        if (id_factor){
          lapply(histvars, FUN = function(histvar){
            pool[pool[[time_name]] == t, (paste("lag_cumavg", i, "_", histvar, sep = "")) :=
                   as.double(tapply(pool[pool[[time_name]] < t-i+1 & get(id_name) %in% current_ids, ][[histvar]],
                                    droplevels(pool[pool[[time_name]] < t-i+1 & get(id_name) %in% current_ids][[id_name]]), FUN = mean))]
          })
        } else {
          lapply(histvars, FUN = function(histvar){
            pool[pool[[time_name]] == t, (paste("lag_cumavg", i, "_", histvar, sep = "")) :=
                   as.double(tapply(pool[pool[[time_name]] < t-i+1 & get(id_name) %in% current_ids, ][[histvar]],
                                    pool[pool[[time_name]] < t-i+1 & get(id_name) %in% current_ids][[id_name]], FUN = mean))]
          })
        }
      } else {
        if (below_zero_indicator){
          denom <- t - i + 1 - min(pool[[time_name]])
        } else {
          denom <- t - i + 1
        }
        # The lagged cumulative average variable was already created
        # Therefore, can use recursive formula
        for (histvar in histvars){
          pool[get(time_name) == t, (paste("lag_cumavg", i, "_", histvar, sep = "")) :=
                 as.double((pool[get(id_name) %in% current_ids & get(time_name) == (t-1)][[paste("lag_cumavg", i, "_", histvar, sep = "")]] * (denom - 1)
                            + pool[get(id_name) %in% current_ids & get(time_name) == (t - i)][[histvar]]) / denom)
          ]
        }
      }
    }
  }
}


#' Create Visit Sum Covariate
#'
#' This internal function assists in the implementation of a visit process by creating a covariate,
#' \code{ts_visit}, that counts the number of visits in the past \code{max_visits} time points. If this
#' number is greater than 0, then the individual has not missed more than the maximum number
#' of visits.
#'
#' @param t          Integer specifying the current time index.
#' @param time_name  Character string specifying the name of the time variable in \code{pool}.
#' @param pool       Data table containing all information prior to time \eqn{t} (\eqn{t} noninclusive).
#' @param histvars   Vector of character strings specifying the names of the variables for which lagged cumulative averages are to
#'                   be created.
#' @param id_name    Character string specifying the name of the ID variable in \code{pool}.
#' @param max_visits A vector of one or more values denoting the maximum number of times
#'                   a binary covariate representing a visit process may be missed before
#'                   the individual is censored from the data (in the observed data) or
#'                   a visit is forced (in the simulated data). Multiple values exist in the
#'                   vector when the modeling of more than covariate is attached to a visit
#'                   process.
#' @return No value is returned. The data table \code{pool} is modified in place.
#' @import data.table
#' @keywords internal
visit_sum <- function(pool, histvars, time_name, t, id_name, max_visits){

  current_ids <- unique(pool[pool[[time_name]] ==t ,get(id_name)])
  if(t == 0){
      lapply(histvars,FUN=function(histvar){
      pool[ pool[[time_name]] == t   , c( paste("ts_", histvar, sep = "") ,  paste("lag1_ts_", histvar, sep = "") ):= list(as.integer(0) , as.integer(0))]
  ##    pool[ pool[[time_name]] == t , paste("lag1_ts_", histvar, sep = "") := as.integer(0)]
      })

  } else if (t > 0) {
    lapply(histvars,FUN=function(histvar) {
    pool[pool[[time_name]] == t   ,
         paste('lag1_ts_',histvar,sep="") :=  as.integer( pool[ pool[[time_name]] == (t-1) & get(id_name) %in% current_ids , eval(parse(text=paste('ts_',histvar,sep=""))) ] )]
    pool[( pool[[time_name]] == t & pool[[histvar]] == 1 )  , paste('ts_',histvar,sep="") := as.integer(0)]
    pool[( pool[[time_name]] == t & pool[[histvar]] == 0 )  , paste('ts_',histvar,sep="") := as.integer(eval(parse(text=paste('lag1_ts_',histvar,sep=""))) + 1)]
    }
    )
  }
}





#' Generates Functions of History of Existing Covariates
#'
#' This internal function applies the history functions to create new columns in an input data table containing new variables that are functions
#' of the histories of existing variables in the dataset.
#'
#' @param t          Integer specifying the current time index.
#' @param time_name  Character string specifying the name of the time variable in \code{pool}.
#' @param pool       Data table containing all information prior to time \eqn{t} (\eqn{t} noninclusive).
#' @param histvars   List of vectors. The kth vector specifies the names of the variables for which the kth history function
#'                   in \code{histories} is to be applied.
#' @param histvals   List of length two. The first element is a numeric vector specifying the lags used in the model statements (e.g., if \code{lag1_varname} and \code{lag2_varname} were included in the model statements, this vector would be \code{c(1,2)}). The second element is a numeric vector specifying the lag averages used in the model statements.
#' @param histories  Vector of history functions to apply to the variables specified in \code{histvars}.
#' @param id         Character string specifying the name of the ID variable in \code{pool}.
#' @param max_visits A vector of one or more values denoting the maximum number of times
#'                   a binary covariate representing a visit process may be missed before
#'                   the individual is censored from the data (in the observed data) or
#'                   a visit is forced (in the simulated data). Multiple values exist in the
#'                   vector when the modeling of more than covariate is attached to a visit
#'                   process. A value of \code{NA} should be provided when there is no visit process.
#' @param baselags   Logical scalar for specifying the convention used for lagi and lag_cumavgi terms in the model statements when pre-baseline times are not
#'                   included in \code{obs_data} and when the current time index, \eqn{t}, is such that \eqn{t < i}. If this argument is set to \code{FALSE}, the value
#'                   of all lagi and lag_cumavgi terms in this context are set to 0 (for non-categorical covariates) or the reference
#'                   level (for categorical covariates). If this argument is set to \code{TRUE}, the value of lagi and lag_cumavgi terms
#'                   are set to their values at time 0. The default is \code{FALSE}.
#' @param below_zero_indicator Logical scalar indicating whether the observed data set contains rows for time \eqn{t < 0}.
#' @return No value is returned. The data table \code{pool} is modified in place.
#' @keywords internal

make_histories <- function(pool, histvars, histvals, histories, time_name, t, id,
                           max_visits, baselags, below_zero_indicator){
  if (!is.na(histvars[1]) && !is.na(histories[1])){
    lapply(seq_along(histories), FUN = function(i){
      if (isTRUE(all.equal(histories[[i]], visit_sum))){
        visit_sum(pool = pool, histvars = histvars[[i]], time_name = time_name,
                       t = t, id_name = id, max_visits = max_visits[i])
      } else if (isTRUE(all.equal(histories[[i]], lagged))){
        lagged(pool = pool, histvars = histvars[[i]],
               histvals = histvals$lag_indicator,
               time_name = time_name, t = t, id_name = id,
               baselags = baselags, below_zero_indicator = below_zero_indicator)
      } else if (isTRUE(all.equal(histories[[i]], lagavg))){
        lagavg(pool = pool, histvars = histvars[[i]],
               histvals = histvals$lagavg_indicator,
               time_name = time_name, t = t, id_name = id,
               baselags = baselags, below_zero_indicator = below_zero_indicator)
      } else if (isTRUE(all.equal(histories[[i]], cumavg))){
        cumavg(pool = pool, histvars = histvars[[i]],
               time_name = time_name, t = t, id_name = id,
               below_zero_indicator = below_zero_indicator)
      } else {
        histories[[i]](pool = pool, histvars = histvars[[i]], time_name = time_name,
                       t = t, id_name = id)
      }
    })
  }
}
