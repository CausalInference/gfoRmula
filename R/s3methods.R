# Copyright (c) 2019 The President and Fellows of Harvard College

#' Print method for objects of class "gformula_survival"
#'
#' Print method for objects of class "gformula_survival".
#'
#' @param x Object of class "gformula_survival".
#' @param coefficients Logical scalar indicating whether to print the model coefficients. The default is \code{FALSE}.
#' @param stderrs Logical scalar indicating whether to print the standard error of the model coefficients. The default is \code{FALSE}.
#' @param rmses Logical scalar indicating whether to print the model root mean square errors (RMSEs). The default is \code{FALSE}.
#' @param hazardratio Logical scalar indicating whether to print the hazard ratio between two interventions (if computed). If bootstrapping was used in \code{\link{gformula_survival}}, 95\% confidence intervals will be given. The default is \code{FALSE}.
#' @param ... Other arguments.
#' @return No value is returned.
#' @seealso \code{\link{gformula_survival}}
#'
#' @export

print.gformula_survival <- function(x, coefficients = FALSE, stderrs = FALSE,
                                    rmses = FALSE, hazardratio = FALSE, ...) {
  if (!inherits(x, "gformula_survival")){
    stop("Argument 'x' must be an object of class \"gformula_survival\".")
  }

  cat(x$header, '\n\n')
  print(x$result[k == max(k)], row.names = FALSE, col.names='top')

  if (rmses){
    cat('\n\n RMSE Values\n')
    print(x$rmses)
  }
  if (coefficients){
    cat('\n\n Model Coefficients\n')
    print(x$coeffs)
  }
  if (stderrs){
    cat('\n\n Standard Errors\n')
    print(x$stderrs)
  }
  if (hazardratio & !is.na(x$hazardratio_val[1])){
    cat('\n\n Hazard ratio\n')
    print(x$hazardratio_val)
  }
}



#' Print method for objects of class "gformula_continuous_eof"
#'
#' Print method for objects of class "gformula_continuous_eof".
#'
#' @param x Object of class "gformula_continuous_eof".
#' @param coefficients Logical scalar indicating whether to  print the model coefficients. The default is \code{FALSE}.
#' @param stderrs Logical scalar indicating whether to print the standard error of the model coefficients. The default is \code{FALSE}.
#' @param rmses Logical scalar indicating whether to print the model root mean square errors (RMSEs). The default is \code{FALSE}.
#' @param ... Other arguments.
#' @return No value is returned.
#' @seealso \code{\link{gformula_continuous_eof}}
#'
#' @export

print.gformula_continuous_eof <- function(x, coefficients = FALSE,
                                          stderrs = FALSE, rmses = FALSE, ...) {
  if (!inherits(x, "gformula_continuous_eof")){
    stop("Argument 'x' must be an object of class \"gformula_continuous_eof\".")
  }

  cat(x$header, '\n\n')
  print(x$result[k == max(k)], row.names = FALSE, col.names='top')

  if (rmses){
    cat('\n\n RMSE Values\n')
    print(x$rmses)
  }
  if (coefficients){
    cat('\n\n Model Coefficients\n')
    print(x$coeffs)
  }
  if (stderrs){
    cat('\n\n Standard Errors\n')
    print(x$stderrs)
  }
}



#' Print method for objects of class "gformula_binary_eof"
#'
#' Print method for objects of class "gformula_binary_eof".
#'
#' @param x Object of class "gformula_binary_eof".
#' @param coefficients Logical scalar indicating whether to  print the model coefficients. The default is \code{FALSE}.
#' @param stderrs Logical scalar indicating whether to print the standard error of the model coefficients. The default is \code{FALSE}.
#' @param rmses Logical scalar indicating whether to print the model root mean square errors (RMSEs). The default is \code{FALSE}.
#' @param ... Other arguments.
#' @return No value is returned.
#' @seealso \code{\link{gformula_binary_eof}}
#'
#' @export

print.gformula_binary_eof <- function(x, coefficients = FALSE, stderrs = FALSE,
                                      rmses = FALSE, ...) {
  if (!inherits(x, "gformula_binary_eof")){
    stop("Argument 'x' must be an object of class \"gformula_binary_eof\".")
  }

  cat(x$header, '\n\n')
  print(x$result[k == max(k)], row.names = FALSE, col.names='top')

  if (rmses){
    cat('\n\n RMSE Values\n')
    print(x$rmses)
  }
  if (coefficients){
    cat('\n\n Model Coefficients\n')
    print(x$coeffs)
  }
  if (stderrs){
    cat('\n\n Standard Errors\n')
    print(x$stderrs)
  }
}



#' Plot method for objects of class "gformula_survival"
#'
#' This function generates graphs of the mean simulated vs. observed values at each time point of the
#' time-varying covariates, risk, and survival under the natural course. For categorical covariates,
#' the observed and simulated counts of the levels of the factors are plotted at each time point.
#'
#' @param x Object of class "gformula_survival".
#' @param covnames Vector of character strings specifying the names of the time-varying covariates to be plotted. The ordering of covariates given here is used in the plot grid. Time-varying covariates of type \code{"categorical time"} cannot be included. To plot none of the time-varying covariates, set this argument to \code{NA}. By default, this argument is set equal to the \code{covnames} argument used in \code{\link{gformula_survival}}, where covariates of type 'categorical time' are removed.
#' @param risk Logical scalar indicating whether to include a plot for the risk. The default is \code{TRUE}.
#' @param survival Logical scalar indicating whether to include a plot for the survival. The default is \code{FALSE}.
#' @param ncol Number of columns in the plot grid. By default, two columns are used when there is at least two plots.
#' @param nrow Number of rows in the plot grid. By default, a maximum of six rows is used and additional plots are included in subsequent pages.
#' @param common.legend Logical scalar indicating whether to include a legend. The default is \code{TRUE}.
#' @param legend Character string specifying the legend position. Valid values are \code{"top"}, \code{"bottom"}, \code{"left"}, \code{"right"}, and \code{"none"}. The default is \code{"bottom"}.
#' @param xlab Character string for the x axes of all plots. By default, this argument is set to the \code{time_name} argument specified in \code{\link{gformula_survival}}.
#' @param ylab_cov Vector of character strings for the y axes of the plots for the covariates. This argument must be the same length as \code{covnames}. The i-th element of this argument corresponds to the plot for the i-th element of \code{covnames}.
#' @param ylab_risk Character string for the y axis of the plot for the risk (if applicable). The default is \code{"risk"}.
#' @param ylab_surv Character string for the y axis of the plot for the survival (if applicable). The default is \code{"survival"}.
#' @param pos_risk Integer specifying the position at which to order the risk plot (if applicable). By default, this argument is set to the number of plots in the grid minus one (i.e., orders the risk plot second last).
#' @param pos_surv Integer specifying the position at which to order the survival plot (if applicable). By default, this argument is set to the number of plots in the grid (i.e., orders the survival plot last).
#' @param ci_risk Logical scalar specifying whether to include error bars for the 95\% confidence intervals of the estimated risk under the natural course. This argument is only effective if the argument \code{nsamples} was set to a positive value in \code{\link{gformula_survival}}. The default is \code{TRUE}.
#' @param ... Other arguments, which are passed to \code{\link[ggpubr]{ggarrange}}.
#' @seealso \code{\link{gformula_survival}}
#'
#' @export

plot.gformula_survival <- function(x, covnames = NULL, risk = TRUE,
                                   survival= FALSE, ncol = NULL, nrow = NULL,
                                   common.legend = TRUE, legend = 'bottom',
                                   xlab = NULL, ylab_cov = NULL,
                                   ylab_risk = 'risk', ylab_surv = 'survival',
                                   pos_risk = NULL, pos_surv = NULL,
                                   ci_risk = FALSE, ...){
  if (!inherits(x, "gformula_survival")){
    stop("Argument 'x' must be an object of class \"gformula_survival\".")
  }
  if (x$time_points == 1){
    stop("Plotting is not available when time_points = 1.")
  }

  if (is.null(covnames)){
    covnames <- x$covnames[x$covtypes != 'categorical time']
    covtypes <- x$covtypes[x$covtypes != 'categorical time']
  } else if (length(covnames)>1 || !is.na(covnames)){
    if (any(!(covnames %in% x$covnames))) {
      stop(paste0("Argument 'covnames' includes the name(s) of covariate(s) not found in x$covnames. These include: ", covnames[!(covnames %in% x$covnames)],".", collapse = ", "))
    }
    covtypes <- x$covtypes[match(covnames, x$covnames)]
    if (any(covtypes == 'categorical time')){
      stop("Argument 'covnames' includes the name of a time-varying covariate that is of type 'categorical time', which is not allowed.")
    }
  }

  if (!is.null(ylab_cov) && length(ylab_cov) != length(covnames)) {
    stop("Argument 'ylab_cov' must have the same length as 'covnames'. If the default value of 'covnames' was used, note that its length equals the number of the time-varying covariates that are not of type 'categorical time'.")
  }

  if (is.null(xlab)){
    xlab <- x$time_name
  }
  if (is.null(ylab_cov)){
    ylab_cov <- ifelse(covtypes == 'categorical', paste(covnames, "counts"), paste(covnames))
  }

  if (length(covnames)>1 || !is.na(covnames)){
    cvgrphs <- get_cvgrphs(x, covnames, covtypes, xlab, ylab_cov)
  } else {
    cvgrphs <- list()
  }
  outgrphs <- get_outgrphs(x, risk, survival, xlab, ylab_risk, ylab_surv, ci_risk)
  plotlist <- c(cvgrphs, outgrphs)

  plot_order <- rep(NA, length(plotlist))
  if (risk & !survival){
    if (is.null(pos_risk)){
      pos_risk <- length(plotlist)
    }
    plot_order[c(pos_risk)] <- c("risk")
  } else if (!risk & survival){
    if (is.null(pos_surv)){
      pos_surv <- length(plotlist)
    }
    plot_order[c(pos_surv)] <- c("survival")
  } else if (risk & survival) {
    if (is.null(pos_risk)){
      pos_risk <- length(plotlist)-1
    }
    if (is.null(pos_surv)){
      pos_surv <- length(plotlist)
    }
    plot_order[c(pos_risk, pos_surv)] <- c("risk", "survival")
  }
  plot_order[is.na(plot_order)] <- covnames
  plotlist <- plotlist[plot_order]

  if (is.null(ncol)){
    if (length(plotlist) > 1){
      ncol <- 2
    } else {
      ncol <- 1
    }
  }
  if (is.null(nrow)){
    nrow <- ceiling(length(plotlist) / ncol)
  }


  while (grDevices::dev.cur() > 1){
    grDevices::dev.off()
  }
  ggpubr::ggarrange(plotlist = plotlist, ncol = ncol, nrow = nrow,
                    common.legend = common.legend, legend = legend, ...)
}



#' Plot method for objects of class "gformula_continuous_eof"
#'
#' This function generates graphs of the mean simulated vs. observed values at each time point of the
#' time-varying covariates under the natural course. For categorical covariates,
#' the observed and simulated counts of the levels of the factors are plotted at each time point.
#'
#'
#' @param x Object of class "gformula_continuous_eof".
#' @param covnames Vector of character strings specifying the names of the time-varying covariates to be plotted. The ordering of covariates given here is used in the plot grid. Time-varying covariates of type \code{"categorical time"} cannot be included. By default, this argument is set equal to the \code{covnames} argument used in \code{\link{gformula_continuous_eof}}, where covariates of type \code{"categorical time"} are removed.
#' @param ncol Number of columns in the plot grid. By default, two columns are used when there is at least two plots.
#' @param nrow Number of rows in the plot grid. By default, a maximum of six rows is used and additional plots are included in subsequent pages.
#' @param common.legend Logical scalar indicating whether to include a legend. The default is \code{TRUE}.
#' @param legend Character string specifying the legend position. Valid values are \code{"top"}, \code{"bottom"}, \code{"left"}, \code{"right"}, and \code{"none"}. The default is \code{"bottom"}.
#' @param xlab Character string for the x axes of all plots. By default, this argument is set to the \code{time_name} argument specified in \code{\link{gformula_continuous_eof}}.
#' @param ylab_cov Vector of character strings for the y axes of the plots for the covariates. This argument must be the same length as \code{covnames}. The i-th element of this argument corresponds to the plot for the i-th element of \code{covnames}.
#' @param ... Other arguments, which are passed to \code{\link[ggpubr]{ggarrange}}.
#' @seealso \code{\link{gformula_continuous_eof}}
#'
#' @export

plot.gformula_continuous_eof <- function(x, covnames = NULL, ncol = NULL, nrow = NULL,
                                         common.legend = TRUE, legend = 'bottom',
                                         xlab = NULL, ylab_cov = NULL, ...){
  if (!inherits(x, "gformula_continuous_eof")){
    stop("Argument 'x' must be an object of class \"gformula_continuous_eof\".")
  }
  if (x$time_points == 1){
    stop("Plotting is not available when time_points = 1.")
  }

  if (is.null(covnames)){
    covnames <- x$covnames[x$covtypes != 'categorical time']
    covtypes <- x$covtypes[x$covtypes != 'categorical time']
  } else if (any(covnames %in% x$covnames)){
    if (any(!(covnames %in% x$covnames))) {
      stop(paste0("Argument 'covnames' includes the name(s) of covariate(s) not found in x$covnames. These include: ", covnames[!(covnames %in% x$covnames)],".", collapse = ", "))
    }
    covtypes <- x$covtypes[match(covnames, x$covnames)]
    if (any(covtypes == 'categorical time')){
      stop("Argument 'covnames' includes the name of a time-varying covariate that is of type 'categorical time', which is not allowed.")
    }
  } else {
    stop("Argument 'covnames' must include the name of at least one time-varying covariate.")
  }

  if (!is.null(ylab_cov) && length(ylab_cov) != length(covnames)) {
    stop("Argument 'ylab_cov' must have the same length as 'covnames'. If the default value of 'covnames' was used, note that its length equals the number of the time-varying covariates that are not of type 'categorical time'.")
  }

  if (is.null(xlab)){
    xlab <- x$time_name
  }
  if (is.null(ylab_cov)){
    ylab_cov <- ifelse(covtypes == 'categorical', paste(covnames, "counts"), paste(covnames))
  }

  plotlist <- get_cvgrphs(x, covnames, covtypes, xlab, ylab_cov)

  if (is.null(ncol)){
    if (length(plotlist) > 1){
      ncol <- 2
    } else {
      ncol <- 1
    }
  }
  if (is.null(nrow)){
    nrow <- ceiling(length(plotlist) / ncol)
  }

  while (grDevices::dev.cur() > 1){
    grDevices::dev.off()
  }
  ggpubr::ggarrange(plotlist = plotlist, ncol = ncol, nrow = nrow,
                    common.legend = common.legend, legend = legend, ...)
}



#' Plot method for objects of class "gformula_binary_eof"
#'
#' This function generates graphs of the mean simulated vs. observed values at each time point of the
#' time-varying covariates under the natural course. For categorical covariates,
#' the observed and simulated counts of the levels of the factors are plotted at each time point.
#'
#' @param x Object of class "gformula_binary_eof".
#' @param covnames Vector of character strings specifying the names of the time-varying covariates to be plotted. The ordering of covariates given here is used in the plot grid. Time-varying covariates of type \code{"categorical time"} cannot be included. By default, this argument is set equal to the \code{covnames} argument used in \code{\link{gformula_binary_eof}}, where covariates of type \code{"categorical time"} are removed.
#' @param ncol Number of columns in the plot grid. By default, two columns are used when there is at least two plots.
#' @param nrow Number of rows in the plot grid. By default, a maximum of six rows is used and additional plots are included in subsequent pages.
#' @param common.legend Logical scalar indicating whether to include a legend. The default is \code{TRUE}.
#' @param legend Character string specifying the legend position. Valid values are \code{"top"}, \code{"bottom"}, \code{"left"}, \code{"right"}, and \code{"none"}. The default is \code{"bottom"}.
#' @param xlab Character string for the x axes of all plots. By default, this argument is set to the \code{time_name} argument specified in \code{\link{gformula_binary_eof}}.
#' @param ylab_cov Vector of character strings for the y axes of the plots for the covariates. This argument must be the same length as \code{covnames}. The i-th element of this argument corresponds to the plot for the i-th element of \code{covnames}.
#' @param ... Other arguments, which are passed to \code{\link[ggpubr]{ggarrange}}.
#' @seealso \code{\link{gformula_binary_eof}}
#'
#' @export

plot.gformula_binary_eof <- function(x, covnames = NULL, ncol = NULL, nrow = NULL,
                                     common.legend = TRUE, legend = 'bottom',
                                     xlab = NULL, ylab_cov = NULL, ...){
  if (!inherits(x, "gformula_binary_eof")){
    stop("Argument 'x' must be an object of class \"gformula_binary_eof\".")
  }
  if (x$time_points == 1){
    stop("Plotting is not available when time_points = 1.")
  }

  if (is.null(covnames)){
    covnames <- x$covnames[x$covtypes != 'categorical time']
    covtypes <- x$covtypes[x$covtypes != 'categorical time']
  } else if (any(covnames %in% x$covnames)){
    if (any(!(covnames %in% x$covnames))) {
      stop(paste0("Argument 'covnames' includes the name(s) of covariate(s) not found in x$covnames. These include: ", covnames[!(covnames %in% x$covnames)],".", collapse = ", "))
    }
    covtypes <- x$covtypes[match(covnames, x$covnames)]
    if (any(covtypes == 'categorical time')){
      stop("Argument 'covnames' includes the name of a time-varying covariate that is of type 'categorical time', which is not allowed.")
    }
  } else {
    stop("Argument 'covnames' must include the name of at least one time-varying covariate.")
  }

  if (!is.null(ylab_cov) && length(ylab_cov) != length(covnames)) {
    stop("Argument 'ylab_cov' must have the same length as 'covnames'. If the default value of 'covnames' was used, note that its length equals the number of the time-varying covariates that are not of type 'categorical time'.")
  }

  if (is.null(xlab)){
    xlab <- x$time_name
  }
  if (is.null(ylab_cov)){
    ylab_cov <- ifelse(covtypes == 'categorical', paste(covnames, "counts"), paste(covnames))
  }

  plotlist <- get_cvgrphs(x, covnames, covtypes, xlab, ylab_cov)

  if (is.null(ncol)){
    if (length(plotlist) > 1){
      ncol <- 2
    } else {
      ncol <- 1
    }
  }
  if (is.null(nrow)){
    nrow <- ceiling(length(plotlist) / ncol)
  }

  while (grDevices::dev.cur() > 1){
    grDevices::dev.off()
  }
  ggpubr::ggarrange(plotlist = plotlist, ncol = ncol, nrow = nrow,
                    common.legend = common.legend, legend = legend, ...)
}
