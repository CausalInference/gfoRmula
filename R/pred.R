#' Fit GLM on Covariate
#'
#' This internal function fits a generalized linear model (GLM) for a single covariate using the observed data.
#'
#' @param covparams   List of vectors, where each vector contains information for
#'                    one parameter used in the modeling of the time-varying covariates (e.g.,
#'                    model statement, family, link function, etc.).
#' @param covlink     Vector of link functions.
#' @param covfam      Vector of character strings specifying the names of the family functions to be
#'                    used for fitting the GLM.
#' @param obs_data    Data on which the model is fit.
#' @param j           Integer specifying the index of the covariate.
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal

fit_glm <- function(covparams, covlink = NA, covfam, obs_data, j){
  # Get model parameters
  covmodels <- covparams$covmodels
  if (!is.null(covparams$covlink)){
    covlink <- covparams$covlink
  }

  if (is.na(covlink[j])){
    famtext <- paste(covfam, "()", sep="")
  } else {
    famtext <- paste(covfam, "(link = ", covlink[j], ")", sep="")
  }

  # Fit GLM for covariate using user-specified formula
  if (!is.null(covparams$control) && !is.na(covparams$control[j])){
    fit <- stats::glm(stats::as.formula(paste(covmodels[j])), family = eval(parse(text = famtext)),
                      data = obs_data, y = TRUE, control = covparams$control[j])
  } else {
    fit <- stats::glm(stats::as.formula(paste(covmodels[j])), family = eval(parse(text = famtext)),
                      data = obs_data, y = TRUE)
  }

  fit$rmse <- add_rmse(fit)
  fit$stderrs <- add_stderr(fit)
  fit <- trim_glm(fit)
  return (fit)
}

#' Fit Multinomial Model on Covariate
#'
#' This internal function fits a multinomial regression model for a categorical covariate using the observed data.
#'
#' @param covparams   List of vectors, where each vector contains information for
#'                    one parameter used in the modeling of the time-varying covariates (e.g.,
#'                    model statement, family, link function, etc.).
#' @param obs_data    Data on which the model is fit.
#' @param j           Integer specifying the index of the covariate.
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal

fit_multinomial <- function(covparams, obs_data, j){
  covmodels <- covparams$covmodels
  if (!is.null(covparams$control) && !is.na(covparams$control[j])){
    args <- c(list(formula = stats::as.formula(paste(covmodels[j])),
                   data = obs_data), trace = FALSE, covparams$control[j])
    fit <- do.call(nnet::multinom, args = args)
  } else {
    fit <- nnet::multinom(stats::as.formula(paste(covmodels[j])),
                          data = obs_data, trace = FALSE)
  }

  fit$stderr <- add_stderr(fit)
  fit <- trim_multinom(fit)
  return (fit)
}

#' Fit Zero-Inflated Normal Model on Covariate
#'
#' This internal function models a zero-inflated normal distribution through the combined
#' use of a generalized linear model (GLM) fit on a zero vs. non-zero indicator
#' and a GLM fit on all non-zero values.
#'
#' @param covparams   List of vectors, where each vector contains information for
#'                    one parameter used in the modeling of the time-varying covariates (e.g.,
#'                    model statement, family, link function, etc.).
#' @param covlink     Vector of link functions.
#' @param covname     Name of the covariate at index \eqn{j}.
#' @param obs_data    Data on which the model is fit.
#' @param j           Integer specifying the index of the covariate.
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table

fit_zeroinfl_normal <- function(covparams, covlink = NA, covname, obs_data, j){
  covmodels <- covparams$covmodels
  if (!is.null(covparams$covlink)){
    covlink <- covparams$covlink
  }

  if (is.na(covlink[j])){
    famtext <- "gaussian()"
  } else {
    famtext <- paste("gaussian(link = ", covlink[j], ")", sep = "")
  }

  obs_data[, paste("I_", covname, sep = "")] <- obs_data[[covname]]
  obs_data[[paste("I_", covname, sep = "")]][obs_data[[covname]] != 0] <- 1

  # Take log to ensure that no negative values are predicted
  obs_data[, paste("log_", covname, sep = "")] <- 0
  obs_data[obs_data[[covname]] != 0][, paste("log_", covname, sep = "")] <-
    log(obs_data[obs_data[[covname]] != 0][[covname]])
  # Fit binomial model on indicator of whether covariate is 0 or not
  if (!is.null(covparams$control) && !is.na(covparams$control[j]) &&
      !is.null(covparams$control[j][[1]])){
    fit1 <- stats::glm(stats::as.formula(paste("I_", covmodels[j], sep = "")), family = stats::binomial(),
                       control = covparams$control[j][[1]], data = obs_data)
  } else {
    fit1 <- stats::glm(stats::as.formula(paste("I_", covmodels[j], sep = "")), family = stats::binomial(),
                       data = obs_data)
  }


  # Fit Gaussian model on data points for which covariate does not equal 0
  if (!is.null(covparams$control) && !is.na(covparams$control[j]) &&
      !is.null(covparams$control[j][[2]])){
    fit2 <- stats::glm(stats::as.formula(paste("log_", covmodels[j], sep = "")),
                       family = eval(parse(text = famtext)),
                       control = covparams$control[j][[2]],
                       data = obs_data[obs_data[[covname]] != 0])
  } else {
    fit2 <- stats::glm(stats::as.formula(paste("log_", covmodels[j], sep = "")),
                       family = eval(parse(text = famtext)),
                       data = obs_data[obs_data[[covname]] != 0])
  }

  fit1$rmse <- add_rmse(fit1)
  fit2$rmse <- add_rmse(fit2)
  fit1$stderr <- add_stderr(fit1)
  fit2$stderr <- add_stderr(fit2)
  fit1 <- trim_glm(fit1)
  fit2 <- trim_glm(fit2)

  return (list(fit1, fit2))
}

#' Fit Bounded Normal Model on Covariate
#'
#' This internal function models a covariate using a "bounded normal" distribution
#' by first standardizing the covariate values to the range [0, 1], noninclusive,
#' then fitting a generalized linear model (GLM) under the Gaussian family function.
#'
#' @param covparams   List of vectors, where each vector contains information for
#'                    one parameter used in the modeling of the time-varying covariates (e.g.,
#'                    model statement, family, link function, etc.).
#' @param covlink     Vector of link functions.
#' @param covname     Name of the covariate at index \eqn{j}.
#' @param obs_data    Data on which the model is fit.
#' @param j           Integer specifying the index of the covariate.
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table

fit_bounded_continuous <- function(covparams, covlink = NA, covname, obs_data, j){

  covmodels <- covparams$covmodels
  if (!is.null(covparams$covlink)){
    covlink <- covparams$covlink
  }

  obs_data[, paste("norm_", covname, sep = "")] <-
    (obs_data[[covname]] - min(obs_data[[covname]]))/(max(obs_data[[covname]]) - min(obs_data[[covname]]))
  if (!is.na(covlink[j])){
    if (!is.null(covparams$control) && !is.na(covparams$control[j])){
      fit <- stats::glm(stats::as.formula(paste("norm_", covmodels[j], sep = "")),
                        family = stats::gaussian(link = covlink[j]), data = obs_data, y = TRUE,
                        control = covparams$control[j])
    } else {
      fit <- stats::glm(stats::as.formula(paste("norm_", covmodels[j], sep = "")),
                        family = stats::gaussian(link = covlink[j]), data = obs_data, y = TRUE)
    }
  } else {
    if (!is.null(covparams$control) && !is.na(covparams$control[j])){
      fit <- stats::glm(stats::as.formula(paste("norm_", covmodels[j], sep = "")), family = stats::gaussian(),
                        data = obs_data, y = TRUE, control = covparams$control[j])
    } else {
      fit <- stats::glm(stats::as.formula(paste("norm_", covmodels[j], sep = "")), family = stats::gaussian(),
                        data = obs_data, y = TRUE)
    }
  }
  fit$rmse <- add_rmse(fit)
  fit$stderr <- add_stderr(fit)
  fit <- trim_glm(fit)
  return (fit)
}

#' Fit Truncated Normal Model on Covariate
#'
#' This internal function models a covariate using a normal distribution truncated on one
#' side at a user-specified cutoff.
#'
#' @param covparams   List of vectors, where each vector contains information for
#'                    one parameter used in the modeling of the time-varying covariates (e.g.,
#'                    model statement, family, link function, etc.).
#' @param obs_data    Data on which the model is fit.
#' @param j           Integer specifying the index of the covariate.
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal

fit_trunc_normal <- function(covparams, obs_data, j){
  covmodels <- covparams$covmodels
  point <- covparams$point[j]
  direction <- covparams$direction[j]
  if (!is.null(covparams$control) && !is.na(covparams$control[j])){
    args <- c(list(formula = stats::as.formula(paste(covmodels[j])),
                   data = obs_data, point = point, direction = direction,
                   y = TRUE), covparams$control[j])
    fit <- do.call(truncreg::truncreg, args = args)
  } else {
    fit <- truncreg::truncreg(stats::as.formula(paste(covmodels[j])), data = obs_data, point = point,
                              direction = direction, y = TRUE)
  }

  fit$rmse <- add_rmse(fit)
  fit$stderr <- add_stderr(fit)
  fit <- trim_truncreg(fit)
  return (fit)
}

#' Fit Covariate Models
#'
#' This internal function fits a model for each covariate using the observed data.
#'
#' @param covparams       List of vectors, where each vector contains information for
#'                        one parameter used in the modeling of the time-varying covariates (e.g.,
#'                        model statement, family, link function, etc.). Each vector
#'                        must be the same length as \code{covnames} and in the same order.
#' @param covnames        Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes        Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param covfits_custom  Vector containing custom fit functions for time-varying covariates that
#'                        do not fall within the pre-defined covariate types. It should be in
#'                        the same order \code{covnames}. If a custom fit function is not
#'                        required for a particular covariate (e.g., if the first
#'                        covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                        index should be set to \code{NA}.
#' @param restrictions    List of vectors. Each vector contains as its first entry a covariate for which
#'                        \emph{a priori} knowledge of its distribution is available; its second entry a condition
#'                        under which no knowledge of its distribution is available and that must be \code{TRUE}
#'                        for the distribution of that covariate given that condition to be estimated via a parametric
#'                        model or other fitting procedure; its third entry a function for estimating the distribution
#'                        of that covariate given the condition in the second entry is false such that \emph{a priori} knowledge
#'                        of the covariate distribution is available; and its fourth entry a value used by the function in the
#'                        third entry. The default is \code{NA}.
#' @param time_name       Character string specifying the name of the time variable in \code{obs_data}.
#' @param obs_data        Data on which the models are fit.
#' @return                A list of fitted models, one for each covariate in \code{covnames}.
#' @keywords internal
#' @import data.table

pred_fun_cov <- function(covparams, covnames, covtypes, covfits_custom,
                         restrictions, time_name, obs_data){
  if (!is.na(restrictions[[1]][[1]])){ # Check for restrictions
    # Create list of covariates whose modeling is affected by restrictions
    restrictnames <- lapply(seq_along(restrictions), FUN = function(r){
      restrictions[[r]][[1]]
    })
    # Create list of conditions where covariates are modeled
    conditions <- lapply(seq_along(restrictions), FUN = function(r){
      restrictions[[r]][[2]]
    })
  } else {
    restrictnames <- conditions <- NA
  }

  # Create list of models for covariates

  subdata <- obs_data[obs_data[[time_name]] > 0]

  fits <- lapply(seq_along(covnames), FUN = function(j){
    if (!is.na(restrictions[[1]][[1]])){
      if (covnames[j] %in% restrictnames){
        i <- which(restrictnames %in% covnames[j])
        condition <- ""
        if (length(i) > 1){
          for (k in i){
            condition_var <- sub("<.*$", "", restrictions[[k]][[2]])
            condition_var <- sub(">.*$", "", condition_var)
            condition_var <- sub("=.*$", "", condition_var)
            condition_var <- sub("!.*$", "", condition_var)
            condition_var <- sub("%in%.*$", "", condition_var)
            condition_var <- sub(" ", "", condition_var)
            if (condition_var %in% names(obs_data)){
              if (condition[1] == ""){
                condition <- conditions[k]
              } else {
                condition <- paste(condition, conditions[k], sep = "||")
              }
            }
          }
        } else {
          condition_var <- sub("<.*$", "", restrictions[[i]][[2]])
          condition_var <- sub(">.*$", "", condition_var)
          condition_var <- sub("=.*$", "", condition_var)
          condition_var <- sub("!.*$", "", condition_var)
          condition_var <- sub("%in%.*$", "", condition_var)
          condition_var <- sub(" ", "", condition_var)
          if (condition_var %in% names(obs_data)){
            condition <- conditions[i]
          }
        }
        if (condition != ""){
          subdata <- subset(subdata, eval(parse(text = condition)))
        }
      }
    }
    if (covtypes[j] == 'binary'){
      fit_glm(covparams = covparams, covfam = 'binomial', obs_data = subdata, j = j)
    } else if (covtypes[j] == 'normal'){
      fit_glm(covparams = covparams, covfam = 'gaussian', obs_data = subdata, j = j)
    } else if (covtypes[j] == 'categorical'){
      fit_multinomial(covparams, obs_data = subdata, j = j)
    } else if (covtypes[j] == 'zero-inflated normal'){
      fit_zeroinfl_normal(covparams, covname = covnames[j], obs_data = subdata, j = j)
    } else if (covtypes[j] == 'bounded normal'){
      fit_bounded_continuous(covparams, covname = covnames[j], obs_data = subdata, j = j)
    } else if (covtypes[j] == 'truncated normal'){
      fit_trunc_normal(covparams = covparams, obs_data = subdata, j = j)
    } else if (covtypes[j] == 'custom'){
      covfits_custom[[j]](covparams, covname = covnames[j], obs_data = subdata, j = j)
    }
  })
  return (fits)
}

#' Fit Outcome Model
#'
#' This internal function fits a generalized linear model (GLM) for the outcome variable using the observed data.
#'
#' @param model          Model statement for the outcome variable.
#' @param yrestrictions  List of vectors. Each vector containins as its first entry
#'                       a condition and its second entry an integer. When the
#'                       condition is \code{TRUE}, the outcome variable is simulated
#'                       according to the fitted model; when the condition is \code{FALSE},
#'                       the outcome variable takes on the value in the second entry.
#' @param outcome_type   Character string specifying the "type" of the outcome. The possible "types" are: \code{"survival"}, \code{"continuous_eof"}, and \code{"binary_eof"}.
#' @param outcome_name   Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param time_name      Character string specifying the name of the time variable in \code{obs_data}.
#' @param obs_data       Data on which the model is fit.
#' @return               Fitted model for the outcome variable.
#' @keywords internal
#' @import data.table

pred_fun_Y <- function(model, yrestrictions, outcome_type, outcome_name, time_name, obs_data){
  if (outcome_type == 'continuous' || outcome_type == 'continuous_eof'){
    outcome_fam <- stats::gaussian()
  } else if (outcome_type == 'survival' || outcome_type == 'binary_eof'){
    outcome_fam <- stats::binomial()
  }
  if (grepl('eof', outcome_type)){
    obs_data <- obs_data[obs_data[[time_name]] == max(unique(obs_data[[time_name]]))]
  }
  if (!is.na(yrestrictions[[1]][[1]])){ # Check for restrictions on outcome variable modeling
    # Set condition where outcome variable is modeled
    condition <- yrestrictions[[1]][1]
    if (length(yrestrictions) > 1){
      # If more than one restriction on outcome variable, set condition such that modeling
      # occurs when any one of the restriction conditions is fulfilled
      for (yrestriction in yrestrictions[-1]){
        condition <- paste(condition, yrestriction[1], sep = "||")
      }
    }
    # Fit GLM for outcome variable using user-specified model
    if (grepl('continuous', outcome_type)){
      obs_data[, paste("norm_", outcome_name, sep = "")] <-
        (obs_data[[outcome_name]] - min(obs_data[[outcome_name]]))/(max(obs_data[[outcome_name]]) - min(obs_data[[outcome_name]]))
      fitY <- stats::glm(stats::as.formula(paste("norm_", model, sep = "")), family = outcome_fam,
                  data = subset(obs_data, eval(parse(text = condition))), y = TRUE)
    } else {
      fitY <- stats::glm(model, family = outcome_fam,
                  data = subset(obs_data, eval(parse(text = condition))), y = TRUE)
    }
  } else { # Case where there are no restrictions on outcome variable
    # Fit GLM for outcome variable using user-specified model and entire dataset
    fitY <- stats::glm(model, family = outcome_fam, data = obs_data, y = TRUE)
  }
  fitY$rmse <- add_rmse(fitY)
  fitY$stderr <- add_stderr(fitY)
  fitY <- trim_glm(fitY)
  return (fitY)
}

#' Fit Competing Event Model
#'
#' This internal function fits a generalized linear model (GLM) for the competing event variable using the observed data.
#'
#' @param model                  Model statement for competing event variable.
#' @param compevent_restrictions List of vectors. Each vector containins as its first entry
#'                               a condition and its second entry an integer. When the
#'                               condition is \code{TRUE}, the competing event variable is simulated
#'                               according to the fitted model; when the condition is \code{FALSE},
#'                               the competing event variable takes on the value in the
#'                               second entry.
#' @param obs_data               Data on which the model is fit.
#' @return                       Fitted model for the competing event variable.
#' @keywords internal
#' @import data.table

pred_fun_D <- function(model, compevent_restrictions, obs_data){
  if (!is.na(compevent_restrictions[[1]][[1]])){ # Check for restrictions on compevent event variable modeling
    # Set condition where competing event variable is modeled
    condition <- compevent_restrictions[[1]][1]
    if (length(compevent_restrictions) > 1){
      # If more than one restriction on compeeting event variable, set condition such
      # that modeling occurs when any one of the restriction conditions is fulfilled
      for (compevent_restriction in compevent_restrictions[-1]){
        condition <- paste(condition, compevent_restriction[1], sep = "||")
      }
    }
    # Fit GLM assuming binomial distribution
    # Restrict data used for modeling to rows where condition is true
    fitD <- stats::glm(model, family = stats::binomial(),
                data = subset(obs_data, eval(parse(text = condition))), y = TRUE)
  } else { # Case where there are no restrictions on competing event variable
    #Fit GLM assuming binomial distribution
    fitD <- stats::glm(model, family = stats::binomial(), data = obs_data, y = TRUE)
  }
  fitD$rmse <- add_rmse(fitD)
  fitD$stderr <- add_stderr(fitD)
  fitD <- trim_glm(fitD)
  return (fitD)
}
