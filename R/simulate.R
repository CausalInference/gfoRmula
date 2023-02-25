#' Simulate Binary Values
#'
#' This internal function simulates covariate values from a binomial distribution.
#'
#' @param x     Integer specifying the number of observations to be simulated.
#' @param size  Integer specifying the number of trials.
#' @param prob  Numeric vector specifying the probabilities.
#' @return      Numeric vector of simulated covariate values under the binomial distribution.
#' @keywords internal
#'
predict_binomial <- function(x, size, prob){
  return (stats::rbinom(x, size, prob))
}

#' Simulate Normal Values
#'
#' This internal function simulates covariate values from a normal distribution.
#'
#' @param x       Integer specifying the number of observations to be simulated.
#' @param mean    Numeric scalar specifying the mean of the distribution.
#' @param est_sd  Numeric scalar specifying the standard deviation of the distribution.
#' @return        Numeric vector of simulated covariate values under the normal distribution.
#' @keywords internal
#'
predict_normal <- function(x, mean, est_sd = NA){
  return (stats::rnorm(x, mean, est_sd))
}


#' Simulate Truncated Normal Values
#'
#' This internal function simulates covariate values from a normal distribution truncated
#' on one side.
#'
#' @param x         Integer specifying the number of observations to be simulated.
#' @param mean      Numeric scalar specifying the mean of the distribution.
#' @param est_sd    Numeric scalar specifying the standard deviation of the distribution.
#' @param a         Numeric scalar specifying the lower bound of truncation.
#' @param b         Numeric scalar specifying the upper bound of truncation.
#' @return          Numeric vector of simulated covariate values under the truncated normal
#'                  distribution.
#' @keywords internal
#'
predict_trunc_normal <- function(x, mean, est_sd, a, b){
  return (truncnorm::rtruncnorm(n = x, mean = mean, sd = est_sd, a = a, b = b))
}

#' Simulate Counterfactual Outcomes Under Intervention
#'
#' This internal function simulates a new dataset containing covariates, outcome probabilities, competing event
#' probabilities (if any), outcomes, and competing events (if any) based on an observed
#' dataset and a user-specified intervention.
#'
#' @param o                       Integer specifying the index of the current intervention.
#' @param fitcov                  List of model fits for the time-varying covariates.
#' @param fitY                    Model fit for the outcome variable.
#' @param fitD                    Model fit for the competing event variable, if any.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#' @param compevent_restrictions  List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the competing event variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the competing event variable takes on the value in the
#'                                second entry.
#' @param restrictions            List of vectors. Each vector contains as its first entry a covariate for which
#'                                \emph{a priori} knowledge of its distribution is available; its second entry a condition
#'                                under which no knowledge of its distribution is available and that must be \code{TRUE}
#'                                for the distribution of that covariate given that condition to be estimated via a parametric
#'                                model or other fitting procedure; its third entry a function for estimating the distribution
#'                                of that covariate given the condition in the second entry is false such that \emph{a priori} knowledge
#'                                of the covariate distribution is available; and its fourth entry a value used by the function in the
#'                                third entry. The default is \code{NA}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param compevent_name          Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}.
#' @param intvars                 List, whose elements are vectors of character strings. The kth vector in \code{intvars} specifies the name(s) of the variable(s) to be intervened
#'                                on in each round of the simulation under the kth intervention in \code{interventions}.
#' @param interventions           List, whose elements are lists of vectors. Each list in \code{interventions} specifies a unique intervention on the relevant variable(s) in \code{intvars}. Each vector contains a function
#'                                implementing a particular intervention on a single variable, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#' @param int_times               List, whose elements are lists of vectors. The kth list in \code{int_times} corresponds to the kth intervention in \code{interventions}. Each vector specifies the time points in which the relevant intervention is applied on the corresponding variable in \code{intvars}.
#'                                When an intervention is not applied, the simulated natural course value is used. By default, this argument is set so that all interventions are applied in all time points.
#' @param histvars                List of vectors. The kth vector specifies the names of the variables for which the kth history function
#'                                in \code{histories} is to be applied.
#' @param histvals                List of length two. The first element is a numeric vector specifying the lags used in the model statements (e.g., if \code{lag1_varname} and \code{lag2_varname} were included in the model statements, this vector would be \code{c(1,2)}). The second element is a numeric vector specifying the lag averages used in the model statements.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}.
#' @param comprisk                Logical scalar indicating the presence of a competing event.
#' @param ranges                  List of vectors. Each vector contains the minimum and
#'                                maximum values of one of the covariates in \code{covnames}.
#' @param outcome_type            Character string specifying the "type" of the outcome. The possible "types" are: \code{"survival"}, \code{"continuous_eof"}, and \code{"binary_eof"}.
#' @param subseed                 Integer specifying the seed for this simulation.
#' @param obs_data                Data table containing the observed data.
#' @param time_points             Number of time points to simulate.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param covnames                Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}.
#' @param max_visits              A vector of one or more values denoting the maximum number of times
#'                                a binary covariate representing a visit process may be missed before
#'                                the individual is censored from the data (in the observed data) or
#'                                a visit is forced (in the simulated data). Multiple values exist in the
#'                                vector when the modeling of more than covariate is attached to a visit
#'                                process.
#' @param baselags                Logical scalar for specifying the convention used for lagi and lag_cumavgi terms in the model statements when pre-baseline times are not
#'                                included in \code{obs_data} and when the current time index, \eqn{t}, is such that \eqn{t < i}. If this argument is set to \code{FALSE}, the value
#'                                of all lagi and lag_cumavgi terms in this context are set to 0 (for non-categorical covariates) or the reference
#'                                level (for categorical covariates). If this argument is set to \code{TRUE}, the value of lagi and lag_cumavgi terms
#'                                are set to their values at time 0. The default is \code{FALSE}.
#' @param below_zero_indicator    Logical scalar indicating whether the observed data set contains rows for time \eqn{t < 0}.
#' @param min_time                Numeric scalar specifying lowest value of time \eqn{t} in the observed data set.
#' @param show_progress           Logical scalar indicating whether to print a progress bar for the number of bootstrap samples completed in the R console. This argument is only applicable when \code{parallel} is set to \code{FALSE} and bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0). The default is \code{TRUE}.
#' @param pb                      Progress bar R6 object. See \code{\link[progress]{progress_bar}} for further details.
#' @param int_visit_type          Vector of logicals. The kth element is a logical specifying whether to carry forward the intervened value (rather than the natural value) of the treatment variables(s) when performing a carry forward restriction type for the kth intervention in \code{interventions}.
#'                                When the kth element is set to \code{FALSE}, the natural value of the treatment variable(s) in the kth intervention in \code{interventions} will be carried forward.
#'                                By default, this argument is set so that the intervened value of the treatment variable(s) is carried forward for all interventions.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#' @return                        A data table containing simulated data under the specified intervention.
#' @keywords internal
#' @import data.table
simulate <- function(o, fitcov, fitY, fitD,
                     yrestrictions, compevent_restrictions, restrictions,
                     outcome_name, compevent_name, time_name,
                     intvars, interventions, int_times, histvars, histvals, histories,
                     comprisk, ranges,
                     outcome_type, subseed, obs_data, time_points, parallel,
                     covnames, covtypes, covparams, covpredict_custom,
                     basecovs, max_visits, baselags, below_zero_indicator,
                     min_time, show_progress, pb, int_visit_type, ...){
  set.seed(subseed)

  # Mechanism of passing intervention variable and intervention is different for parallel
  # and non-parallel versions
  if (parallel){
    intvar <- intvars[[o]]
    intervention <- interventions[[o]]
    int_time <- int_times[[o]]
    int_visit_type <- int_visit_type[o]
  } else {
    intvar <- intvars
    intervention <- interventions
    int_time <- int_times
  }

  if (!is.null(fitcov)){
    rmses <- lapply(seq_along(fitcov), FUN = rmse_calculate, fits = fitcov, covnames = covnames,
                    covtypes = covtypes)
  }

  # Initialize
  ids_unique <- unique(obs_data$newid)
  data_len <- length(ids_unique)
  restrict_ids <- rep(0, data_len)
  restrict_counts <- rep(list(rep(0, data_len)), length(restrictions))
  if (!is.na(restrictions[[1]][[1]])){
    restrict_covs <- lapply(restrictions, FUN = function(restriction){restriction[[1]]})
  }

  # Create histories_int and histvars_int, which are the necessary histories to create after the intervention
  nat_course <- length(intvar == 1) && intvar == 'none'
  if (!nat_course) {
    intvar_vec <- unique(unlist(intvar))
    histvars_int <- histories_int <- rep(list(NA), length(histvars))
    for (l in seq_along(histvars)){
      histvars_temp <- histvars[[l]][histvars[[l]] %in% intvar_vec]
      if (length(histvars_temp) > 0){
        histvars_int[[l]] <- histvars_temp
        histories_int[[l]] <- histories[[l]]
      }
    }
    histvars_int <- histvars_int[!is.na(histvars_int)]
    histories_int <- histories_int[!is.na(histories_int)]
  }

  for (t in ((1:time_points) - 1)){
    if (t == 0){
      # Set simulated covariate values at time t = 0 equal to observed covariate values
      if (!is.na(basecovs[[1]])){
        pool <- obs_data[obs_data[[time_name]] <= t, ][, .SD, .SDcols = c(covnames, basecovs, time_name)]
      } else {
        pool <- obs_data[obs_data[[time_name]] <= t, ][, .SD, .SDcols = c(covnames, time_name)]
      }
      set(pool, j = 'id', value = rep(ids_unique, each = 1 - min_time))
      set(pool, j = 'eligible_pt', value = TRUE)
      if (!is.na(basecovs[[1]])){
        setcolorder(pool, c('id', time_name, covnames, basecovs))
      } else {
        setcolorder(pool, c('id', time_name, covnames))
      }
      newdf <- pool[pool[[time_name]] == 0]
      # Update datatable with specified treatment regime / intervention for this
      # simulation
      if (!nat_course){
        mycols <- match(intvar, names(newdf))
        temp_intvar <- newdf[, ..mycols]

        if (!int_visit_type){
          for (var in intvar){
            newdf[, eval(paste0(var, '_natural')) := newdf[[var]]]
          }
        }
      }
      intfunc(newdf, pool = pool, intervention, intvar, unlist(int_time), time_name, t)
      # Check if intervened
      intervened <- rep(0, times = nrow(newdf))
      if (!nat_course){
        for (var in intvar){
          # Check if the natural value of the intervention variable equals the intervened value
          intervened <- intervened + (abs(temp_intvar[[var]] - newdf[[var]]) > 1e-6)
        }
        intervened <- ifelse(newdf$eligible_pt, intervened >= 1, NA)
      }
      set(newdf, j = 'intervened', value = intervened)

      if (ncol(newdf) > ncol(pool)){
        pool <- rbind(pool[pool[[time_name]] < t], newdf, fill = TRUE)
        pool <- pool[order(id, get(time_name))]
      } else {
        pool[pool[[time_name]] == t] <- newdf
      }
      # Update datatable with new covariates that are functions of history of existing
      # covariates
      make_histories(pool = pool, histvars = histvars, histvals = histvals,
                     histories = histories, time_name = time_name, t = t, id = 'id',
                     max_visits = max_visits, baselags = baselags, below_zero_indicator = below_zero_indicator)
      newdf <- pool[pool[[time_name]] == t]
      # Generate outcome probabilities
      if (outcome_type == 'survival'){
        set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
      } else if (outcome_type == 'continuous_eof'){
        if (t < (time_points - 1)){
          set(newdf, j = 'Ey', value = as.double(NA))
        } else if (t == (time_points - 1)){
          set(newdf, j = 'Ey', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      } else if (outcome_type == 'binary_eof'){
        if (t < (time_points - 1)){
          set(newdf, j = 'Py', value = as.double(NA))
        } else if (t == (time_points - 1)){
          set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      }
      if (!is.na(yrestrictions[[1]][[1]])){ # Check if there are restrictions on outcome
        # variable simulation
        for (yrestriction in yrestrictions){
          # Set non-modeled outcome variable values equal to user-specified value
          if (outcome_type == 'survival'){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))),
                  "Py" := as.double(yrestriction[2])]
          }
        }
      }
      # Simulate outcome variable
      if (outcome_type == 'survival'){
        set(newdf, j = 'Y', value = stats::rbinom(data_len, 1, newdf$Py))
      }
      if (outcome_type == 'survival')
      {
        if (comprisk){
          # Predict competing event probabilities
          set(newdf, j = 'Pd', value = stats::predict(fitD, type = 'response', newdata = newdf))
          if (!is.na(compevent_restrictions[[1]][[1]])){ # Check if there are restrictions
            # on competing event variable simulation
            for (compevent_restriction in compevent_restrictions){
              # Set non-modeled competing event values equal to user-specified value
              newdf[!eval(parse(text = compevent_restriction[1])),
                    "Pd" := as.double(compevent_restriction[2])]
            }
          }
          # Simulate competing event variable
          set(newdf, j = 'D', value = stats::rbinom(data_len, 1, newdf$Pd))
          # Calculate probability of death by main event rather than competing event at
          # time t
          set(newdf, j = 'prodp1', value = newdf$Py * (1 - newdf$Pd))
        } else {
          set(newdf, j = 'D', value = 0)
          # Calculate probability of death by main event without competing event
          if (outcome_type == 'survival'){
            set(newdf, j = 'prodp1', value = newdf$Py)
          }
        }
        set(newdf[newdf$D == 1], j = 'Y', value = NA)
        set(newdf, j = 'prodp0', value = 1 - newdf$Py)
      }
      # If competing event occurs, outcome cannot also occur because
      # both presumably lead to death
      # Calculate probability of survival or death by competing event at time t
      pool <- rbind(pool[pool[[time_name]] < t], newdf, fill = TRUE)
      pool <- pool[order(id, get(time_name))]
      col_types <- sapply(pool, class)
    } else {
      # Set initial simulated values at time t to simulated values at time t - 1, to be
      # updated later
      newdf <- pool[pool[[time_name]] == t - 1]
      set(newdf, j = time_name, value = rep(t, data_len))
      if ('categorical time' %in% covtypes){
        time_name_f <- paste(time_name, "_f", sep = "")
        newdf[, (time_name_f) :=
                obs_data[get(time_name) == t, get(time_name_f)][1]]
      }
      pool <- rbind(newdf, pool)
      make_histories(pool = pool, histvars = histvars, histvals = histvals, histories = histories,
                     time_name = time_name, t = t, id = 'id', max_visits = max_visits,
                     baselags = baselags, below_zero_indicator = below_zero_indicator)
      newdf <- pool[pool[[time_name]] == t]
      for (i in seq_along(covnames)){
        cast <- get(paste0('as.',unname(col_types[covnames[i]])))
        if (covtypes[i] == 'binary'){
          set(newdf, j = covnames[i],
              value = cast(predict_binomial(data_len, 1, stats::predict(fitcov[[i]], type = 'response',
                                                                        newdata = newdf))))
        } else if (covtypes[i] == 'normal'){
          set(newdf, j = covnames[i],
              value = cast(predict_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                       newdata = newdf),
                                     est_sd = rmses[[i]])))
        } else if (covtypes[i] == 'categorical'){
          set(newdf, j = covnames[i],
              value = cast(stats::predict(fitcov[[i]], type = 'class', newdata = newdf)))
        } else if (covtypes[i] == 'zero-inflated normal'){
          set(newdf, j = paste("I_", covnames[i], sep = ""),
              value = predict_binomial(data_len, 1, stats::predict(fitcov[[i]][[1]], type = 'response',
                                                            newdata = newdf)))
          # Where indicator is nonzero, simulate from Gaussian model
          set(newdf, j = covnames[i],
              value = cast(exp(predict_normal(data_len,
                                         stats::predict(fitcov[[i]][[2]], type = 'response',
                                                 newdata = newdf),
                                         est_sd =  rmses[[i]]))))
          set(newdf, j = covnames[i],
              value = cast(newdf[[paste("I_", covnames[i], sep = "")]] * newdf[[covnames[i]]]))
          # Remove indicator
          set(newdf, j = paste("I_", covnames[i], sep = ""), value = NULL)
        } else if (covtypes[i] == 'bounded normal'){
          if (!is.na(restrictions[[1]][[1]])){
            restrictnames <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[1]]})
            # Create list of conditions where covariates are modeled
            conditions <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[2]]
            })
            if (covnames[i] %in% restrictnames){
              j <- which(restrictnames %in% covnames[i])
              condition <- ""
              if (length(j) > 1){
                for (k in j){
                  condition_var <- sub("<.*$", "", restrictions[[k]][[2]])
                  condition_var <- sub(">.*$", "", condition_var)
                  condition_var <- sub("=.*$", "", condition_var)
                  condition_var <- sub("!.*$", "", condition_var)
                  condition_var <- sub("%in%.*$", "", condition_var)
                  condition_var <- sub(" ", "", condition_var)
                  if (condition_var %in% names(obs_data)){
                    if (condition[1] == ""){
                      condition <- conditions[j]
                    } else {
                      condition <- paste(condition, conditions[j], sep = "||")
                    }
                  }
                }
              } else {
                condition_var <- sub("<.*$", "", restrictions[[j]][[2]])
                condition_var <- sub(">.*$", "", condition_var)
                condition_var <- sub("=.*$", "", condition_var)
                condition_var <- sub("!.*$", "", condition_var)
                condition_var <- sub("%in%.*$", "", condition_var)
                condition_var <- sub(" ", "", condition_var)
                if (condition_var %in% names(obs_data)){
                  condition <- conditions[j]
                }
              }
              if (condition[1] != ""){
                sub_obs_data <- subset(obs_data[obs_data[[time_name]] >= 0], eval(parse(text = condition)))
              } else {
                sub_obs_data <- obs_data[obs_data[[time_name]] >= 0]
              }
            } else {
              sub_obs_data <- obs_data[obs_data[[time_name]] >= 0]
            }
          } else {
            sub_obs_data <- obs_data[obs_data[[time_name]] >= 0]
          }
          set(newdf, j = paste("norm_", covnames[i], sep = ""),
              value = predict_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                       newdata = newdf), est_sd = rmses[[i]]))
          set(newdf, j = covnames[i],
              value = cast((newdf[[paste("norm_", covnames[i], sep = "")]] *
                         (max(sub_obs_data[[covnames[i]]]) - min(sub_obs_data[[covnames[i]]]))) +
                min(sub_obs_data[[covnames[i]]])))
          set(newdf, j = paste("norm_", covnames[i], sep = ""), value = NULL)
        } else if (covtypes[i] == 'truncated normal'){
          if (covparams$direction[i] == 'left'){
            set(newdf, j = covnames[i],
                value = cast(predict_trunc_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                               newdata = newdf),
                                             est_sd = rmses[[i]], a = covparams$point[i], b = Inf)))
          } else if (covparams$direction[i] == 'right'){
            set(newdf, j = covnames[i],
                value = cast(predict_trunc_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                               newdata = newdf),
                                             est_sd = rmses[[i]], a = - Inf, b = covparams$point[i])))
          }
        } else if (covtypes[i] == 'custom'){
          if (!is.na(restrictions[[1]][[1]])){
            restrictnames <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[1]]})
            # Create list of conditions where covariates are modeled
            conditions <- lapply(seq_along(restrictions), FUN = function(r){
              restrictions[[r]][[2]]
            })
            if (covnames[i] %in% restrictnames){
              j <- which(restrictnames %in% covnames[i])
              condition <- ""
              if (length(j) > 1){
                for (k in j){
                  condition_var <- sub("<.*$", "", restrictions[[k]][[2]])
                  condition_var <- sub(">.*$", "", condition_var)
                  condition_var <- sub("=.*$", "", condition_var)
                  condition_var <- sub("!.*$", "", condition_var)
                  condition_var <- sub("%in%.*$", "", condition_var)
                  condition_var <- sub(" ", "", condition_var)
                  if (condition_var %in% names(obs_data)){
                    if (condition[1] == ""){
                      condition <- conditions[j]
                    } else {
                      condition <- paste(condition, conditions[j], sep = "||")
                    }
                  }
                }
              } else {
                condition_var <- sub("<.*$", "", restrictions[[j]][[2]])
                condition_var <- sub(">.*$", "", condition_var)
                condition_var <- sub("=.*$", "", condition_var)
                condition_var <- sub("!.*$", "", condition_var)
                condition_var <- sub("%in%.*$", "", condition_var)
                condition_var <- sub(" ", "", condition_var)
                if (condition_var %in% names(obs_data)){
                  condition <- conditions[j]
                }
              }
            }
          }
          set(newdf, j = covnames[i],
              value = cast(covpredict_custom[[i]](obs_data, newdf, fitcov[[i]],time_name, t,
                                             condition, covnames[i], ...)))
        }
        if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
           covtypes[i] == 'truncated normal'){

             if (length(newdf[newdf[[covnames[i]]] < ranges[[i]][1]][[covnames[i]]]) != 0){
                newdf[newdf[[covnames[i]]] < ranges[[i]][1], (covnames[i]) := cast(ranges[[i]][1])]
             }
             if (length(newdf[newdf[[covnames[i]]] > ranges[[i]][2]][[covnames[i]]]) != 0){
                newdf[newdf[[covnames[i]]] > ranges[[i]][2], (covnames[i]) := cast(ranges[[i]][2])]
             }
        } else if (covtypes[i] == 'zero-inflated normal') {
             if (length(newdf[newdf[[covnames[i]]] < ranges[[i]][1] & newdf[[covnames[i]]] > 0][[covnames[i]]]) != 0){
                newdf[newdf[[covnames[i]]] < ranges[[i]][1] & newdf[[covnames[i]]] > 0, (covnames[i]) := cast(ranges[[i]][1])]
             }
             if (length(newdf[newdf[[covnames[i]]] > ranges[[i]][2]][[covnames[i]]]) != 0){
                newdf[newdf[[covnames[i]]] > ranges[[i]][2], (covnames[i]) := cast(ranges[[i]][2])]
             }
        }
        # Check if there are restrictions on covariate simulation
        if (!is.na(restrictions[[1]][[1]])){
          lapply(seq_along(restrictions), FUN = function(r){
            if (restrictions[[r]][[1]] == covnames[i]){
              restrict_ids <- newdf[!eval(parse(text = restrictions[[r]][[2]]))]$id
              if (length(restrict_ids) != 0){
                restrictions[[r]][[3]](newdf, pool[pool[[time_name]] < t & pool[[time_name]] >= 0], restrictions[[r]], time_name, t, int_visit_type, intvar)
              }
            }
          })
        }
        pool[pool[[time_name]] == t] <- newdf
        if (covnames[i] %in% unlist(histvars)){
          ind <- unlist(lapply(histvars, FUN = function(x) {
            covnames[i] %in% x
            }))
          make_histories(pool = pool, histvars = rep(list(covnames[i]), sum(ind)),
                         histvals = histvals, histories = histories[ind],
                         time_name = time_name, t = t, id = 'id', max_visits = max_visits,
                         baselags = baselags, below_zero_indicator = below_zero_indicator)
          newdf <- pool[pool[[time_name]] == t]
        }
      }
      # Update datatable with specified treatment regime / intervention for this
      # simulation
      newdf <- pool[pool[[time_name]] == t]
      if (!nat_course){
        mycols <- match(intvar, names(newdf))
        temp_intvar <- newdf[, ..mycols]

        if (!int_visit_type){
          for (var in intvar){
            newdf[, eval(paste0(var, '_natural')) := newdf[[var]]]
          }
        }
      }
      intfunc(newdf, pool, intervention, intvar, unlist(int_time), time_name, t)
      # Check if intervened
      intervened <- rep(0, times = nrow(newdf))
      if (!nat_course){
        for (var in intvar){
          # Check if the natural value of the intervention variable equals the intervened value
          intervened <- intervened + (abs(temp_intvar[[var]] - newdf[[var]]) > 1e-6)
        }
        intervened <- ifelse(newdf$eligible_pt, intervened >= 1, NA)
      }
      set(newdf, j = 'intervened', value = intervened)

      # Update datatable with new covariates that are functions of history of existing
      # covariates
      pool[pool[[time_name]] == t] <- newdf

      if (!(length(intvar) == 1 && intvar == 'none') && length(histvars_int) > 0){
        make_histories(pool = pool, histvars = histvars_int, histvals = histvals,
                       histories = histories_int, time_name = time_name, t = t, id = 'id',
                       max_visits = max_visits, baselags = baselags, below_zero_indicator = below_zero_indicator)
      }

      newdf <- pool[pool[[time_name]] == t]

      # Predict outcome probabilities
      if (outcome_type == 'survival'){
        if (comprisk){
          # Predict competing event probabilities
          set(newdf, j = 'Pd', value = stats::predict(fitD, type = 'response', newdata = newdf))
          if (!is.na(compevent_restrictions[[1]][[1]])){ # Check if there are restrictions
            # on competing event variable
            # simulation
            for (compevent_restriction in compevent_restrictions){
              # Set non-modeled competing event values equal to user-specified value
              newdf[!eval(parse(text = compevent_restriction[1])), "Pd" := as.double(compevent_restriction[2])]
            }
          }
          # Simulate competing event variable
          set(newdf, j = 'D', value = stats::rbinom(data_len, 1, newdf$Pd))
        } else {
          set(newdf, j = 'D', value = 0)
        }
        set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
      } else if (outcome_type == 'continuous_eof'){
        if (t < (time_points - 1)){
          set(newdf, j = 'Ey', value = as.double(NA))
        } else if (t == (time_points - 1)){
          set(newdf, j = 'Ey', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      } else if (outcome_type == 'binary_eof'){
        if (t < (time_points - 1)){
          set(newdf, j = 'Py', value = as.double(NA))
        } else if (t == (time_points - 1)){
          set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      }
      if (!is.na(yrestrictions[[1]][[1]])){ # Check if there are restrictions on outcome
        # variable simulation
        for (yrestriction in yrestrictions){
          # Set non-modeled outcome variable values equal to user-specified value
          if (outcome_type == 'survival'){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Py" := as.double(yrestriction[2])]
          } else if (outcome_type == 'continuous_eof' && t == (time_points - 1)){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Ey" := as.double(yrestriction[2])]
          } else if (outcome_type == 'binary_eof' && t == (time_points - 1)){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Py" := as.double(yrestriction[2])]
          }
        }
      }
      # Calculate probability of survival or death from competing event (if any) at time t
      if (outcome_type == 'survival'){
        set(newdf, j = 'prodp0', value = 1 - newdf$Py)
      }
      # Simulate outcome variable
      if (outcome_type == 'survival'){
        set(newdf, j = 'Y', value = stats::rbinom(data_len, 1, newdf$Py))
      }
      if (outcome_type == 'survival')
        newdf[newdf$D == 1, 'Y' := NA]
      # If competing event occurs, outcome cannot also occur because
      # both presumably lead to death
      # Calculate probability of death from main event at time t
      if (comprisk){
        set(newdf, j = 'prodp1',
            value = newdf$Py * tapply(pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$prodp0,
                                      pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$id, FUN = prod) *
              tapply(1 - pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$Pd,
                     pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$id,
                     FUN = prod) * (1 - newdf$Pd))
      } else if (outcome_type == 'survival'){
        set(newdf, j = 'prodp1', value = newdf$Py * tapply(pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$prodp0,
                                                           pool[pool[[time_name]] < t & pool[[time_name]] >= 0]$id,
                                                           FUN = prod))
      }
      # Add simulated data for time t to aggregate simulated data over time
      pool[pool[[time_name]] == t] <- newdf
    }
  }
  colnames(pool)[colnames(pool) == time_name] <- 't0'
  setorder(pool, id, t0)
  colnames(pool)[colnames(pool) == 't0'] <- time_name
  # Calculate probabiity of death from main event at or before time t for each individual
  # at each time point
  pool <- pool[pool[[time_name]] >= 0]
  if (outcome_type == 'survival'){
    pool[, 'poprisk' := stats::ave(pool$prodp1, by = pool$id, FUN = cumsum)]
    pool[, 'survival' := 1 - pool$poprisk]
  }
  pool2 <- copy(pool)
  if (show_progress){
    pb$tick()
  }
  return (pool2)
}
