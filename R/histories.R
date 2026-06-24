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
#' outcome_type <- 'survival'
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
#' intervention1.A <- list(static, rep(0, time_points))
#' intervention2.A <- list(static, rep(1, time_points))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_basic <- gformula(obs_data = basicdata_nocomp, id = id,
#'                         time_points = time_points,
#'                         time_name = time_name, covnames = covnames,
#'                         outcome_name = outcome_name,
#'                         outcome_type = outcome_type, covtypes = covtypes,
#'                         covparams = covparams, ymodel = ymodel,
#'                         intervention1.A = intervention1.A,
#'                         intervention2.A = intervention2.A,
#'                         int_descript = int_descript,
#'                         histories = histories, histvars = histvars,
#'                         basecovs = c('L3'), nsimul = nsimul,
#'                         seed = 1234)
#' gform_basic
#' }
#'
#' @import data.table
#' @export
lagged <- function(pool, histvars, histvals, time_name, t, id_name, baselags,
                   below_zero_indicator){
  # Cache the time-t rows once and align source slices to them by id. This
  # assumes one row per id per time, as gformula's pool always has.
  tcol <- pool[[time_name]]
  t_idx <- which(tcol == t)
  if (length(t_idx) == 0) return(invisible(pool))
  t_ids <- pool[[id_name]][t_idx]
  for (i in histvals){
    if (t < i & baselags & !below_zero_indicator){
      src_idx <- which(tcol == 0)
      src_pos <- src_idx[match(t_ids, pool[[id_name]][src_idx])]
      for (histvar in histvars){
        myclass <- paste0('as.', class(pool[[histvar]][t_idx])[1])
        set(pool, i = t_idx, j = paste0("lag", i, "_", histvar),
            value = get(myclass)(pool[[histvar]][src_pos]))
      }
    } else if (t < i & !below_zero_indicator){
      for (histvar in histvars){
        col_t <- pool[[histvar]][t_idx]
        if (is.factor(col_t)){
          reflevel <- levels(col_t)[1]
          set(pool, i = t_idx, j = paste0("lag", i, "_", histvar), value = reflevel)
        } else {
          myclass <- paste0('as.', class(col_t)[1])
          set(pool, i = t_idx, j = paste0("lag", i, "_", histvar), value = get(myclass)(0))
        }
      }
    } else {
      src_idx <- which(tcol == t - i)
      src_pos <- src_idx[match(t_ids, pool[[id_name]][src_idx])]
      for (histvar in histvars){
        myclass <- paste0('as.', class(pool[[histvar]][t_idx])[1])
        set(pool, i = t_idx, j = paste0("lag", i, "_", histvar),
            value = get(myclass)(pool[[histvar]][src_pos]))
      }
    }
  }
  invisible(pool)
}

#' @rdname lagged
#' @export
cumavg <- function(pool, histvars, time_name, t, id_name, below_zero_indicator){
  if (below_zero_indicator){
    denom <- t + 1 - min(pool[[time_name]])
  } else {
    denom <- t + 1
  }
  tcol <- pool[[time_name]]
  t_idx <- which(tcol == t)
  if (length(t_idx) == 0) return(invisible(pool))
  if (t == 0 & !below_zero_indicator){
    # At first time point, set all cumulative averages equal to the actual value of the
    # variable
    for (histvar in histvars){
      set(pool, i = t_idx, j = paste0("cumavg_", histvar),
          value = as.double(pool[[histvar]][t_idx]))
    }
  } else {
    # At subsequent time points, create new column containing calculated cumulative
    # average until that time point
    t_ids <- pool[[id_name]][t_idx]
    colnam <- colnames(pool)
    if (!(paste("cumavg_", "_", histvars[1], sep = "") %in% colnam)){
      # The cumulative average variable was not created yet
      # Therefore, cannot use recursive formula
      keep <- which(tcol <= t)
      keep_ids <- pool[[id_name]][keep]
      for (histvar in histvars){
        gm <- vapply(split(pool[[histvar]][keep], keep_ids), mean, numeric(1))
        set(pool, i = t_idx, j = paste0("cumavg_", histvar),
            value = as.double(gm[as.character(t_ids)]))
      }
    } else {
      # The cumulative average variable was already created
      # Therefore, can use recursive formula
      prev_idx <- which(tcol == (t - 1))
      prev_pos <- prev_idx[match(t_ids, pool[[id_name]][prev_idx])]
      for (histvar in histvars){
        cumcol <- paste0("cumavg_", histvar)
        set(pool, i = t_idx, j = cumcol,
            value = as.double(
              (pool[[cumcol]][prev_pos] * (denom - 1) + pool[[histvar]][t_idx]) / denom
            ))
      }
    }
  }
  invisible(pool)
}

#' @rdname lagged
#' @export
lagavg <- function(pool, histvars, histvals, time_name, t, id_name, baselags,
                   below_zero_indicator){
  tcol <- pool[[time_name]]
  t_idx <- which(tcol == t)
  if (length(t_idx) == 0) return(invisible(pool))
  t_ids <- pool[[id_name]][t_idx]
  min_time <- min(tcol)
  for (i in histvals){
    if (t < i & !below_zero_indicator){
      if (baselags){
        src_idx <- which(tcol == 0)
        src_pos <- src_idx[match(t_ids, pool[[id_name]][src_idx])]
        for (histvar in histvars){
          set(pool, i = t_idx, j = paste0("lag_cumavg", i, "_", histvar),
              value = as.double(pool[[histvar]][src_pos]))
        }
      } else {
        # At time point less than i, set all lagged cumulative averages equal to 0
        for (histvar in histvars){
          set(pool, i = t_idx, j = paste0("lag_cumavg", i, "_", histvar), value = as.double(0))
        }
      }

    } else if (t == (i + min_time)){
      # At time point i, set all lagged cumulative averages equal to the actual value
      # of the variable
      src_idx <- which(tcol == t - i)
      src_pos <- src_idx[match(t_ids, pool[[id_name]][src_idx])]
      for (histvar in histvars){
        set(pool, i = t_idx, j = paste0("lag_cumavg", i, "_", histvar),
            value = as.double(pool[[histvar]][src_pos]))
      }
    } else if (t > (i + min_time)) {

      # At time points after i, create new column containing calculated lagged cumulative
      # average until that time point
      colnam <- colnames(pool)
      if (!(paste("lag_cumavg", i, "_", histvars[1], sep = "") %in% colnam)){
        # The lagged cumulative average variable was not created yet
        # Therefore, cannot use recursive formula
        keep <- which(tcol < t - i + 1)
        keep_ids <- pool[[id_name]][keep]
        for (histvar in histvars){
          gm <- vapply(split(pool[[histvar]][keep], keep_ids), mean, numeric(1))
          set(pool, i = t_idx, j = paste0("lag_cumavg", i, "_", histvar),
              value = as.double(gm[as.character(t_ids)]))
        }
      } else {
        if (below_zero_indicator){
          denom <- t - i + 1 - min_time
        } else {
          denom <- t - i + 1
        }
        # The lagged cumulative average variable was already created
        # Therefore, can use recursive formula
        prev_idx <- which(tcol == (t - 1))
        prev_pos <- prev_idx[match(t_ids, pool[[id_name]][prev_idx])]
        ti_idx <- which(tcol == (t - i))
        ti_pos <- ti_idx[match(t_ids, pool[[id_name]][ti_idx])]
        for (histvar in histvars){
          lacol <- paste0("lag_cumavg", i, "_", histvar)
          set(pool, i = t_idx, j = lacol,
              value = as.double(
                (pool[[lacol]][prev_pos] * (denom - 1) + pool[[histvar]][ti_pos]) / denom
              ))
        }
      }
    }
  }
  invisible(pool)
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

  tcol <- pool[[time_name]]
  t_idx <- which(tcol == t)
  if (length(t_idx) == 0) return(invisible(pool))
  if(t == 0){
    for (histvar in histvars){
      set(pool, i = t_idx, j = paste0("ts_", histvar), value = as.integer(0))
      set(pool, i = t_idx, j = paste0("lag1_ts_", histvar), value = as.integer(0))
    }
  } else if (t > 0) {
    t_ids <- pool[[id_name]][t_idx]
    prev_idx <- which(tcol == (t - 1))
    prev_pos <- prev_idx[match(t_ids, pool[[id_name]][prev_idx])]
    for (histvar in histvars){
      ts_col   <- paste0("ts_", histvar)
      lag_col  <- paste0("lag1_ts_", histvar)
      lag_vals <- as.integer(pool[[ts_col]][prev_pos])
      set(pool, i = t_idx, j = lag_col, value = lag_vals)
      # ts = 0 where the visit occurred this period, lag1_ts + 1 where it did not.
      hv <- pool[[histvar]][t_idx]
      touch1 <- which(hv == 1L)
      touch0 <- which(hv == 0L)
      if (length(touch1)) set(pool, i = t_idx[touch1], j = ts_col, value = 0L)
      if (length(touch0)) set(pool, i = t_idx[touch0], j = ts_col,
                              value = lag_vals[touch0] + 1L)
    }
  }
  invisible(pool)
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
#'                   vector when the modeling of more than one covariate is attached to a visit
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
    # Skip the costlier all.equal() when the functions are identical.
    same_fun <- function(f, g) identical(f, g) || isTRUE(all.equal(f, g))
    for (i in seq_along(histories)){
      hf <- histories[[i]]
      if (same_fun(hf, visit_sum)){
        visit_sum(pool = pool, histvars = histvars[[i]], time_name = time_name,
                       t = t, id_name = id, max_visits = max_visits[i])
      } else if (same_fun(hf, lagged)){
        lagged(pool = pool, histvars = histvars[[i]],
               histvals = histvals$lag_indicator,
               time_name = time_name, t = t, id_name = id,
               baselags = baselags, below_zero_indicator = below_zero_indicator)
      } else if (same_fun(hf, lagavg)){
        lagavg(pool = pool, histvars = histvars[[i]],
               histvals = histvals$lagavg_indicator,
               time_name = time_name, t = t, id_name = id,
               baselags = baselags, below_zero_indicator = below_zero_indicator)
      } else if (same_fun(hf, cumavg)){
        cumavg(pool = pool, histvars = histvars[[i]],
               time_name = time_name, t = t, id_name = id,
               below_zero_indicator = below_zero_indicator)
      } else {
        hf(pool = pool, histvars = histvars[[i]], time_name = time_name,
                       t = t, id_name = id)
      }
    }
  }
}
