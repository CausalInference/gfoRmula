#' Calculate Observed Covariate Means and Risk
#'
#' This internal function calculates the mean observed values of covariates at each time point, as well as mean
#' observed risk.
#'
#' @param outcome_name    Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param compevent_name  Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param time_name       Character string specifying the name of the time variable in \code{obs_data}.
#' @param covnames        Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param comprisk        Logical scalar indicating the presence of a competing event.
#' @param outcome_type    Character string specifying the "type" of the outcome. The possible "types" are: \code{"survival"}, \code{"continuous_eof"}, and \code{"binary_eof"}.
#' @param obs_data        Data table containing the observed data.
#' @return                A list. Its first entry is a list of mean covariate values at each time point;
#'                        its second entry is a vector of the mean observed risk (for \code{"survival"}
#'                        outcome types) or the mean observed outcome (for \code{"continuous_eof"} and
#'                        \code{"binary_eof"} outcome types); for \code{"survival"} outcome types, its
#'                        third entry is a vector of mean observed survival.
#' @keywords internal
#' @import data.table
obs_calculate <- function(outcome_name, compevent_name, time_name, covnames, comprisk,
                          outcome_type, obs_data){
  obs_means <- lapply(covnames, FUN = function(covname){
    tapply(as.numeric(obs_data[[covname]]), obs_data[[time_name]], FUN = mean,
           na.rm = TRUE)
  })
  names(obs_means) <- covnames
  if (outcome_type == 'survival'){
    # Calculate mean observed outcome probability at each time point
    meanPy <- tapply(obs_data[[outcome_name]], obs_data[[time_name]], FUN = mean,
                     na.rm = TRUE)
    if (comprisk){
      # Calculate mean observed competing event probability at each time point
      meanPd <- tapply(obs_data[[compevent_name]], obs_data[[time_name]], FUN = mean,
                       na.rm = TRUE)
    }
    # Initialize
    obs_prodp0 <- rep(NA, length(meanPy))
    obs_prodp1 <- rep(NA, length(meanPy))
    # Calculate mean observed probability of death from main event at first time point
    if (comprisk){
      obs_prodp1[1] <- meanPy[1] * (1 - meanPd[1])
    } else {
      obs_prodp1[1] <- meanPy[1]
    }
    # Calculate mean observed probability of survival at first time point
    obs_prodp0[1] <- 1 - obs_prodp1[1]
    if (length(meanPy) > 1){
      for (k in 2:length(meanPy)){
        # Calculate mean observed probability of survival at time k - 1
        obs_prodp0[k] <- 1 - meanPy[k]
        # Calculated mean observed probability of death from main event at time k - 1
        if (comprisk){
          obs_prodp1[k] <- meanPy[k] * (cumprod(obs_prodp0[1:(k - 1)] * (1-meanPd[k]))[k - 1])
        }
        else {
          obs_prodp1[k] <- meanPy[k] * (cumprod(obs_prodp0[1:(k - 1)])[k - 1])
        }
      }
    }
    # Calculate observed probability of death at or before each time point
    obs_risk <- cumsum(obs_prodp1)
    obs_survival <- cumprod(obs_prodp0)
    return (list(obs_means, obs_risk, obs_survival))
  } else if (outcome_type == 'continuous_eof' || outcome_type == 'binary_eof'){
    meanEy <- tapply(obs_data[[outcome_name]], obs_data[[time_name]], FUN = mean,
                     na.rm = TRUE)
    return (list(obs_means, meanEy))
  }
}

#' Calculate RMSE for Covariate, Outcome, and Competing Risk Models
#'
#' This internal function calculates root mean square error (RMSE) for the covariate (and the outcome and competing event, if applicable) models fit on
#' the observed data.
#'
#' @param i                         Integer specifying the index of \code{fits}.
#' @param fits                      List of fitted models.
#' @param covnames                  Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                  Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param obs_data                  Data table containing the observed data.
#' @param outcome_name              Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param time_name                 Character string specifying the name of the time variable in \code{obs_data}.
#' @param restrictions              List of vectors. Each vector contains as its first entry a covariate for which
#'                                  \emph{a priori} knowledge of its distribution is available; its second entry a condition
#'                                  under which no knowledge of its distribution is available and that must be \code{TRUE}
#'                                  for the distribution of that covariate given that condition to be estimated via a parametric
#'                                  model or other fitting procedure; its third entry a function for estimating the distribution
#'                                  of that covariate given the condition in the second entry is false such that \emph{a priori} knowledge
#'                                  of the covariate distribution is available; and its fourth entry a value used by the function in the
#'                                  third entry. The default is \code{NA}.
#' @param yrestrictions             List of vectors. Each vector containins as its first entry
#'                                  a condition and its second entry an integer. When the
#'                                  condition is \code{TRUE}, the outcome variable is simulated
#'                                  according to the fitted model; when the condition is \code{FALSE},
#'                                  the outcome variable takes on the value in the second entry.
#' @param compevent_restrictions    List of vectors. Each vector containins as its first entry
#'                                  a condition and its second entry an integer. When the
#'                                  condition is \code{TRUE}, the competing event variable is simulated
#'                                  according to the fitted model; when the condition is \code{FALSE},
#'                                  the competing event variable takes on the value in the
#'                                  second entry.
#'
#' @return                          The RMSE for the model.
#' @keywords internal
#' @import data.table

rmse_calculate <- function(i, fits, covnames, covtypes, obs_data, outcome_name, time_name,
                           restrictions, yrestrictions, compevent_restrictions){
  fit <- fits[[i]]
  obs_data <- obs_data[obs_data[[time_name]] > 0]
  if (i <= length(covtypes)){
    if (covtypes[i] == 'normal' || covtypes[i] == 'binary' ||
        covtypes[i] == 'bounded normal' || covtypes[i] == 'truncated normal'){
      fit$rmse
    } else if (covtypes[i] == 'zero-inflated normal') {
      fit <- fits[[i]][[2]]
      return (fit$rmse)
    }
  } else {
    if (i == length(covtypes) + 1) {
      if (!is.na(yrestrictions[[1]][[1]])){ # Check for restrictions on outcome variable modeling
        # Set condition where outcome variable is modeled
        ycondition <- yrestrictions[[1]][1]
        if (length(yrestrictions) > 1){
          # If more than one restriction on outcome variable, set condition such that modeling
          # occurs when any one of the restriction conditions is fulfilled
          for (yrestriction in yrestrictions[-1]){
            ycondition <- paste(ycondition, yrestriction[1], sep = "||")
          }
        }
        obs_data <- subset(obs_data, eval(parse(text = ycondition)))
      }
    }
    if (i == length(covtypes) + 2) {
      if (!is.na(compevent_restrictions[[1]][[1]])){ # Check for restrictions on compevent event variable modeling
        # Set condition where competing event variable is modeled
        dcondition <- compevent_restrictions[[1]][1]
        if (length(compevent_restrictions) > 1){
          # If more than one restriction on compeeting event variable, set condition such
          # that modeling occurs when any one of the restriction conditions is fulfilled
          for (compevent_restriction in compevent_restrictions[-1]){
            dcondition <- paste(dcondition, compevent_restriction[1], sep = "||")
          }
        }
        obs_data <- subset(obs_data, eval(parse(text = dcondition)))
      }
    }
    return (sqrt(mean(obs_data[[outcome_name]] - stats::predict(fit, newdata = obs_data), na.rm = TRUE)^2))
  }
}



#' Get Plotting Information
#'
#' This internal function obtains the data tables necessary for plotting. For continuous and binary covariates, the mean observed and simulated values are obtained for each time point. For categorical covariates, the observed and simulated counts of the levels of the factors are obtained for each time point.  When the outcome is of type \code{"survival"}, the observed and simulated risk and survival are obtained for each time point.
#'
#' @param outcome_name    Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param compevent_name  Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param time_name       Character string specifying the name of the time variable in \code{obs_data}.
#' @param time_points     Number of time points to simulate.
#' @param covnames        Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes        Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param nat_pool        Pooled-over-time data table containing simulated data under the natural course.
#' @param nat_result      Vector containing the mean outcome over all subjects at each time for natural course.
#' @param comprisk        Logical scalar indicating the presence of a competing event.
#' @param outcome_type    Character string specifying the "type" of the outcome. The possible "types" are: \code{"survival"}, \code{"continuous_eof"}, and \code{"binary_eof"}.
#' @param obs_data        Data table containing observed data.
#' @return                A list with the following components:
#' \item{obs_results}{A list of the mean observed values at each time point for covariates and - if the outcome is of type \code{"survival"} - the risk and survival.}
#' \item{dt_cov_plot}{A list of data tables The data tables contain the observed and simulated mean values of the covariates under each time point.}
#' \item{dt_obs_plot}{For outcomes of type \code{"survival"}, a list of data tables. The data tables contain the observed and simulated risks and survival under each time point. For other outcomes, a value of \code{NA} is given.}
#' @keywords internal
#' @import data.table
#' @import ggplot2
get_plot_info <- function(outcome_name, compevent_name, time_name, time_points,
                          covnames, covtypes, nat_pool, nat_result, comprisk,
                          outcome_type, obs_data){

  # Calculate mean observed values at each time point for covariates, risk, and survival
  obs_results <- obs_calculate(outcome_name, compevent_name, time_name, covnames, comprisk,
                               outcome_type, obs_data[obs_data[[time_name]] < time_points])

  # Calculate mean simulated values at each time point for covariates
  sim_results_cov <- lapply(covnames, FUN = function(covname){
    tapply(as.numeric(nat_pool[[covname]]), nat_pool[[time_name]], FUN = mean,
           na.rm = TRUE)
  })
  names(sim_results_cov) <- covnames
  if (outcome_type == 'survival'){
    sim_results_surv <- tapply(nat_pool$survival, nat_pool[[time_name]], FUN = mean)
    sim_results <- list(sim_results_cov, nat_result, sim_results_surv)
  } else {
    sim_results <- list(sim_results_cov)
  }

  # Generate data tables for plotting each covariate
  dt_cov_plot <- lapply(1:length(covnames), FUN = function(i){
    covname <- covnames[i]
    if (covtypes[i] == 'categorical'){
      sub_obs_data <- obs_data[obs_data[[time_name]] < time_points,
                               summary(eval(parse(text = covname))), by = time_name]
      sub_obs_data[, (covname) := rep(levels(obs_data[[covname]]), time_points)]
      sub_obs_data[, 'legend' := 'nonparametric estimates']
      sub_est_data <- nat_pool[, summary(eval(parse(text = covname))), by = time_name]
      sub_est_data[, (covname) := rep(levels(nat_pool[[covname]]), time_points)]
      sub_est_data[, 'legend' := 'parametric g-formula estimates']
      comb_data <- rbind(sub_est_data, sub_obs_data)
      names(comb_data)[names(comb_data) == time_name] <- "t0"
      comb_data
    } else if (covtypes[i] == 'categorical time'){
      NA
    } else {
      obs_cov_data <- data.table(t0 = 0:(time_points - 1), cov = obs_results[[1]][[covname]],
                                 legend = "nonparametric estimates")
      est_cov_data <- data.table(t0 = 0:(time_points - 1), cov = sim_results[[1]][[covname]],
                                 legend = "parametric g-formula estimates")
      comb_data <- rbind(obs_cov_data, est_cov_data)
      colnames(comb_data)[colnames(comb_data) == 't0'] <- time_name
      comb_data
    }
  })
  names(dt_cov_plot) <- covnames

  # Generate data tables for plotting outcomes for survival objects
  if (outcome_type == 'survival'){
    obs_risk_data <- data.table(t0 = 0:(time_points - 1),
                                risk = obs_results[[2]],
                                survival = obs_results[[3]],
                                legend = "nonparametric estimates")
    est_risk_data <- data.table(t0 = 0:(time_points - 1),
                                risk = sim_results[[2]],
                                survival = sim_results[[3]],
                                legend = "parametric g-formula estimates")
    dt_out_plot <- rbind(obs_risk_data, est_risk_data)
    colnames(dt_out_plot)[colnames(dt_out_plot) == 't0'] <- time_name
  } else {
    dt_out_plot <- NA
  }

  return(list(obs_results = obs_results, dt_cov_plot = dt_cov_plot,
              dt_out_plot = dt_out_plot))
}

#' Get Covariate Plots
#'
#' This internal function obtains the covariate plots for \code{\link{gformula_survival}}, \code{\link{gformula_continuous_eof}}, and \code{\link{gformula_binary_eof}}.
#'
#' @param x              Object of class "gformula_survival", "gformula_continuous_eof", or "gformula_binary_eof".
#' @param covnames       Vector of character strings specifying the names of the time-varying covariates to be plotted. The ordering of covariates given here is used in the plot grid.
#' @param covtypes       Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param xlab           Character string for the x axes of all plots.
#' @param ylab_cov       Vector of character strings for the y axes of the plots for the covariates. This argument must be the same length as \code{covnames}. The i-th element of this argument corresponds to the plot for the i-th element of \code{covnames}.
#' @return               A list of covariate plots.
#' @keywords internal
#' @import data.table
#' @import ggplot2
get_cvgrphs <- function(x, covnames, covtypes, xlab, ylab_cov){

  cvgrphs <- lapply(1:length(covnames), FUN = function(i){
    covname <- covnames[i]
    covtype <- covtypes[i]
    comb_cov_data <- x$dt_cov_plot[[covname]]

    if (covtype == "categorical"){
      ggplot(data = comb_cov_data, aes_string(x = covname, y = "V1", fill = "legend")) +
        geom_bar(stat = 'identity', position = 'dodge') + facet_grid(~ t0) + xlab(xlab) +
        ylab(ylab_cov[i])
    } else {
      ggplot(data = comb_cov_data, aes_string(x = x$time_name, y = "cov", color = "legend")) +
        geom_point() + geom_line() + xlab(xlab) + ylab(ylab_cov[i])
    }
    })
  names(cvgrphs) <- covnames
  return(cvgrphs)
}


#' Get Risk and Survival Plots
#'
#' This internal function obtains the risk and survival plots for \code{\link{gformula_survival}}.
#'
#' @param x              Object of class "gformula_survival".
#' @param risk           Logical scalar indicating whether to include a plot for the risk.
#' @param survival       Logical scalar indicating whether to include a plot for the survival.
#' @param xlab           Character string for the x axes of all plots.
#' @param ylab_risk      Character string for the y axis of the plot for the risk (if applicable).
#' @param ylab_surv      Character string for the y axis of the plot for the survival (if applicable).
#' @param ci_risk        Logical scalar specifying whether to include error bars for the 95\% confidence intervals of the estimated risk under the natural course. This argument is only effective if the argument \code{nsamples} was set to a positive value in \code{\link{gformula_survival}}.
#' @return               A list of plots for the risk and survival.
#' @keywords internal
#' @import data.table
#' @import ggplot2
get_outgrphs <- function(x, risk, survival, xlab, ylab_risk, ylab_surv, ci_risk){

  risk_outgrph <-
    ggplot(data = x$dt_out_plot, aes_string(x = x$time_name, y = "risk", color = "legend")) +
    geom_point() + geom_line() + xlab(xlab) + ylab(ylab_risk)
  surv_outgrph <-
    ggplot(data = x$dt_out_plot, aes_string(x = x$time_name, y = "survival", color = "legend")) +
    geom_point() + geom_line() + xlab(xlab) + ylab(ylab_surv)

  if (x$nsamples > 0 & ci_risk){
    nat_data <- x$result[x$result$Interv. == 0]
    risk_outgrph <-
      risk_outgrph + geom_errorbar(aes(ymin = rep(nat_data$`Risk lower 95% CI`, 2),
                                       ymax = rep(nat_data$`Risk upper 95% CI`, 2)),
                                   width = 0.3)
  }

  outgrphs <- list()
  if (risk){
    outgrphs[['risk']] <- risk_outgrph
  }
  if (survival){
    outgrphs[['survival']] <- surv_outgrph
  }

  return(outgrphs)
}
