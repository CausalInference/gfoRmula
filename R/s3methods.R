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
#' print(gform_basic)
#' }
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
#' @examples
#' ## Estimating the effect of treatment strategies on the mean of a continuous
#' ## end of follow-up outcome
#' \donttest{
#' library('Hmisc')
#' id <- 'id'
#' time_name <- 't0'
#' covnames <- c('L1', 'L2', 'A')
#' outcome_name <- 'Y'
#' covtypes <- c('categorical', 'normal', 'binary')
#' histories <- c(lagged)
#' histvars <- list(c('A', 'L1', 'L2'))
#' covparams <- list(covmodels = c(L1 ~ lag1_A + lag1_L1 + L3 + t0 +
#'                                   rcspline.eval(lag1_L2, knots = c(-1, 0, 1)),
#'                                 L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + L3 + t0,
#'                                 A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0))
#' ymodel <- Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3
#' intvars <- list('A', 'A')
#' interventions <- list(list(c(static, rep(0, 7))),
#'                       list(c(static, rep(1, 7))))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_cont_eof <- gformula_continuous_eof(obs_data = continuous_eofdata,
#'                                           id = id,
#'                                           time_name = time_name,
#'                                           covnames = covnames,
#'                                           outcome_name = outcome_name,
#'                                           covtypes = covtypes,
#'                                           covparams = covparams, ymodel = ymodel,
#'                                           intvars = intvars,
#'                                           interventions = interventions,
#'                                           int_descript = int_descript,
#'                                           histories = histories, histvars = histvars,
#'                                           basecovs = c("L3"),
#'                                           nsimul = nsimul, seed = 1234)
#' print(gform_cont_eof)
#' }
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
#' @examples
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
#' print(gform_bin_eof)
#' }
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
#' @return An object of class "ggarrange". See documentation of \code{\link[ggpubr]{ggarrange}}.
#' @seealso \code{\link{gformula_survival}}
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
#' plot(gform_basic)
#' }
#'
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
#' @return An object of class "ggarrange". See documentation of \code{\link[ggpubr]{ggarrange}}.
#' @seealso \code{\link{gformula_continuous_eof}}
#'
#' @examples
#' ## Estimating the effect of treatment strategies on the mean of a continuous
#' ## end of follow-up outcome
#' \donttest{
#' library('Hmisc')
#' id <- 'id'
#' time_name <- 't0'
#' covnames <- c('L1', 'L2', 'A')
#' outcome_name <- 'Y'
#' covtypes <- c('categorical', 'normal', 'binary')
#' histories <- c(lagged)
#' histvars <- list(c('A', 'L1', 'L2'))
#' covparams <- list(covmodels = c(L1 ~ lag1_A + lag1_L1 + L3 + t0 +
#'                                   rcspline.eval(lag1_L2, knots = c(-1, 0, 1)),
#'                                 L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + L3 + t0,
#'                                 A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0))
#' ymodel <- Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3
#' intvars <- list('A', 'A')
#' interventions <- list(list(c(static, rep(0, 7))),
#'                       list(c(static, rep(1, 7))))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_cont_eof <- gformula_continuous_eof(obs_data = continuous_eofdata,
#'                                           id = id,
#'                                           time_name = time_name,
#'                                           covnames = covnames,
#'                                           outcome_name = outcome_name,
#'                                           covtypes = covtypes,
#'                                           covparams = covparams, ymodel = ymodel,
#'                                           intvars = intvars,
#'                                           interventions = interventions,
#'                                           int_descript = int_descript,
#'                                           histories = histories, histvars = histvars,
#'                                           basecovs = c("L3"),
#'                                           nsimul = nsimul, seed = 1234)
#' plot(gform_cont_eof)
#' }
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
#' @return An object of class "ggarrange". See documentation of \code{\link[ggpubr]{ggarrange}}.
#' @seealso \code{\link{gformula_binary_eof}}
#'
#' @examples
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
#' plot(gform_bin_eof)
#' }
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



#' Coefficient method for objects of class "gformula"
#'
#' This function extracts the coefficients of the fitted models for the
#' time-varying covariates, outcome, and compevent event (if applicable).
#'
#' @param object Object of class "gformula".
#' @param ... Other arguments.
#' @return If \code{bootdiag} was set to \code{FALSE} in \code{\link{gformula}},
#' this function returns a list of the coefficients of the fitted models to the
#' observed data set. If bootstrapping was used and \code{bootdiag} was set to \code{TRUE} in
#' \code{\link{gformula}}, this function returns a list described as follows.
#' The first element (named 'Original sample') is a list of the coefficients of
#' the fitted models to the observed data set. The kth element (named 'Bootstrap
#' sample k-1') is a list of the coefficients of the fitted models corresponding
#' to the k-1th bootstrap sample.
#' @seealso \code{\link{gformula}}
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
#' intvars <- list('A', 'A')
#' interventions <- list(list(c(static, rep(0, time_points))),
#'                       list(c(static, rep(1, time_points))))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_basic <- gformula(obs_data = basicdata_nocomp, id = id,
#'                         time_points = time_points,
#'                         time_name = time_name, covnames = covnames,
#'                         outcome_name = outcome_name,
#'                         outcome_type = outcome_type, covtypes = covtypes,
#'                         covparams = covparams, ymodel = ymodel,
#'                         intvars = intvars,
#'                         interventions = interventions,
#'                         int_descript = int_descript,
#'                         histories = histories, histvars = histvars,
#'                         basecovs = c('L3'), nsimul = nsimul,
#'                         seed = 1234)
#' coefs(gform_basic)
#' }
#'
#' @export
coef.gformula <- function(object, ...){
  if (!inherits(object, "gformula")){
    stop("Argument 'object' must be an object of class \"gformula\".")
  }
  res <- object$coeffs
  if (!is.null(object$bootcoeffs)){
    res <- c(list(res), object$bootcoeffs)
    names(res)[1] <- c('Original sample')
    names(res)[2:length(res)] <- c(paste('Bootstrap sample', 1:length(object$bootcoeffs)))
  }
  return(res)
}



#' Variance-covariance method for objects of class "gformula"
#'
#' This function extracts the variance-covariance matrices of the parameters of
#' the fitted models for the time-varying covariates, outcome, and competing
#' event (if applicable).
#'
#' @param object Object of class "gformula".
#' @param ... Other arguments.
#' @return If \code{bootdiag} was set to \code{FALSE} in \code{\link{gformula}},
#' this function returns a list of the variance-covariance matrices of the
#' parameters of the fitted models to the observed data set. If bootstrapping
#' was used and \code{bootdiag} was set to \code{TRUE} in
#' \code{\link{gformula}}, this function returns a list described as follows.
#' The first element (named 'Original sample') is a list of the
#' variance-covariance matrices of the parameters of the fitted models to the
#' observed data set. The kth element (named 'Bootstrap sample k-1') is a list
#' of the variance-covariance matrices of the parameters of the fitted models
#' corresponding to the k-1th bootstrap sample.
#' @seealso \code{\link{gformula}}
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
#' intvars <- list('A', 'A')
#' interventions <- list(list(c(static, rep(0, time_points))),
#'                       list(c(static, rep(1, time_points))))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_basic <- gformula(obs_data = basicdata_nocomp, id = id,
#'                         time_points = time_points,
#'                         time_name = time_name, covnames = covnames,
#'                         outcome_name = outcome_name,
#'                         outcome_type = outcome_type, covtypes = covtypes,
#'                         covparams = covparams, ymodel = ymodel,
#'                         intvars = intvars,
#'                         interventions = interventions,
#'                         int_descript = int_descript,
#'                         histories = histories, histvars = histvars,
#'                         basecovs = c('L3'), nsimul = nsimul,
#'                         seed = 1234)
#' vcov(gform_basic)
#' }
#'
#' @export
vcov.gformula <- function(object, ...){
  if (!inherits(object, "gformula")){
    stop("Argument 'object' must be an object of class \"gformula\".")
  }
  res <- object$vcovs
  if (!is.null(object$bootvcovs)){
    res <- c(list(res), object$bootvcovs)
    names(res)[1] <- c('Original sample')
    names(res)[2:length(res)] <- c(paste('Bootstrap sample', 1:length(object$bootvcovs)))
  }
  return(res)
}
