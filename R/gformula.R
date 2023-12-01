#' Estimation of Survival Outcome, Continuous End-of-Follow-Up Outcome, or Binary End-of-Follow-Up Outcome Under the Parametric G-Formula
#'
#' Based on an observed data set, this function estimates the risk over time (for survival outcomes),
#' outcome mean at end-of-follow-up (for continuous end-of-follow-up outcomes), or outcome probability at
#' end-of-follow-up (for binary end-of-follow-up outcomes) under multiple user-specified interventions using
#' the parametric g-formula. See McGrath et al. (2020) for further details concerning the application and
#' implementation of the parametric g-formula.
#'
#' To assess model misspecification in the parametric g-formula, users can obtain inverse probability (IP) weighted estimates of the natural course risk and/or means of the time-varying covariates from the observed data.
#' See Chiu et al. (In press) for details.
#' In addition to the general requirements described in McGrath et al. (2020), the requirements for the input data set and the call to the gformula function for such analyses are described below.
#'
#' Users need to include a column in \code{obs_data} with a time-varying censoring variable.
#' Users need to indicate the name of the censoring variable and a model statement for the censoring variable with parameters \code{censor_name} and \code{censor_model}, respectively.
#' When competing events are present, users need to include a column in \code{obs_data} with a time-varying indicator of the competing event variable and need to indicate the name of the competing event variable and the corresponding model statement with parameters \code{compevent_name} and \code{compevent_model}, respectively.
#' Users need to indicate whether to treat competing events as censoring events with the \code{compevent_cens} parameter.
#' Finally, users can specify how to truncate IP weights with the \code{ipw_cutoff_quantile} or \code{ipw_cutoff_value} parameters.
#'
#' In addition to the package output described in McGrath et al. (2020), the output will display estimates of the "cumulative percent intervened on" and the "average percent intervened on". When using a custom intervention function, users need to specify whether each individual at that time point is eligible to contribute person-time to the percent intervened on calculations. Specifically, this must be specified in the \code{eligible_pt} column of \code{newdf}. By default, \code{eligible_pt} is set to \code{TRUE} for each individual at each time point in custom interventions.
#'
#' @param id                      Character string specifying the name of the ID variable in \code{obs_data}.
#' @param time_points             Number of time points to simulate. By default, this argument is set equal to the maximum
#'                                number of records that \code{obs_data} contains for any individual plus 1.
#' @param obs_data                Data table containing the observed data.
#' @param seed                    Starting seed for simulations and bootstrapping.
#' @param nsimul                  Number of subjects for whom to simulate data. By default, this argument is set
#'                                equal to the number of subjects in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param compevent_name          Character string specifying the name of the competing event variable in \code{obs_data}. Only applicable for survival outcomes.
#' @param censor_name             Character string specifying the name of the censoring variable in \code{obs_data}. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param censor_model            Model statement for the censoring variable. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param outcome_type            Character string specifying the "type" of outcome. The possible "types" are: \code{"survival"}, \code{"continuous_eof"}, and \code{"binary_eof"}.
#' @param intvars                 List, whose elements are vectors of character strings. The kth vector in \code{intvars} specifies the name(s) of the variable(s) to be intervened
#'                                on in each round of the simulation under the kth intervention in \code{interventions}.
#' @param interventions           List, whose elements are lists of vectors. Each list in \code{interventions} specifies a unique intervention on the relevant variable(s) in \code{intvars}. Each vector contains a function
#'                                implementing a particular intervention on a single variable, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#' @param int_times               List, whose elements are lists of vectors. The kth list in \code{int_times} corresponds to the kth intervention in \code{interventions}. Each vector specifies the time points in which the relevant intervention is applied on the corresponding variable in \code{intvars}.
#'                                When an intervention is not applied, the simulated natural course value is used. By default, this argument is set so that all interventions are applied in all time points.
#' @param int_descript            Vector of character strings, each describing an intervention. It must
#'                                be in same order as the entries in \code{interventions}.
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the risk ratio and risk difference. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}. The default is 0.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}. The default is \code{NA}.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}. The default is \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}. These covariates are not simulated using a model but rather carry their value over all time points from the first time point of \code{obs_data}. These covariates should not be included in \code{covnames}. The default is \code{NA}.
#' @param histvars                List of vectors. The kth vector specifies the names of the variables for which the kth history function
#'                                in \code{histories} is to be applied.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}. The default is \code{NA}.
#' @param ymodel                  Model statement for the outcome variable.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#'                                The default is \code{NA}.
#' @param compevent_restrictions  List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the competing event variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the competing event variable takes on the value in the
#'                                second entry. The default is \code{NA}. Only applicable for survival outcomes.
#' @param restrictions            List of vectors. Each vector contains as its first entry a covariate for which
#'                                \emph{a priori} knowledge of its distribution is available; its second entry a condition
#'                                under which no knowledge of its distribution is available and that must be \code{TRUE}
#'                                for the distribution of that covariate given that condition to be estimated via a parametric
#'                                model or other fitting procedure; its third entry a function for estimating the distribution
#'                                of that covariate given the condition in the second entry is false such that \emph{a priori} knowledge
#'                                of the covariate distribution is available; and its fourth entry a value used by the function in the
#'                                third entry. The default is \code{NA}.
#' @param visitprocess            List of vectors. Each vector contains as its first entry
#'                                the covariate name of a visit process; its second entry
#'                                the name of a covariate whose modeling depends on the
#'                                visit process; and its third entry the maximum number
#'                                of consecutive visits that can be missed before an
#'                                individual is censored. The default is \code{NA}.
#' @param compevent_model         Model statement for the competing event variable. The default is \code{NA}. Only applicable for survival outcomes.
#' @param compevent_cens          Logical scalar indicating whether to treat competing events as censoring events.
#'                                This argument is only applicable for survival outcomes and when a competing even model is supplied (i.e., \code{compevent_name} and \code{compevent_model} are specified).
#'                                If this argument is set to \code{TRUE}, the competing event model will only be used to construct inverse probability weights to estimate the natural course means / risk from the observed data.
#'                                If this argument is set to \code{FALSE}, the competing event model will be used in the parametric g-formula estimates of the risk and will not be used to construct inverse probability weights.
#'                                See "Details". The default is \code{FALSE}.
#' @param intcomp                 List of two numbers indicating a pair of interventions to be compared by a hazard ratio.
#'                                The default is \code{NA}, resulting in no hazard ratio calculation.
#' @param baselags                Logical scalar for specifying the convention used for lagi and lag_cumavgi terms in the model statements when pre-baseline times are not
#'                                included in \code{obs_data} and when the current time index, \eqn{t}, is such that \eqn{t < i}. If this argument is set to \code{FALSE}, the value
#'                                of all lagi and lag_cumavgi terms in this context are set to 0 (for non-categorical covariates) or the reference
#'                                level (for categorical covariates). If this argument is set to \code{TRUE}, the value of lagi and lag_cumavgi terms
#'                                are set to their values at time 0. The default is \code{FALSE}.
#' @param nsamples                Integer specifying the number of bootstrap samples to generate.
#'                                The default is 0.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of CPU cores to use in parallel
#'                                simulation. This argument is required when parallel is set to \code{TRUE}.
#'                                In many applications, users may wish to set this argument equal to \code{parallel::detectCores() - 1}.
#' @param sim_data_b              Logical scalar indicating whether to return the simulated data set. If bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0), this argument must be set to \code{FALSE}. The default is \code{FALSE}.
#' @param ci_method               Character string specifying the method for calculating the bootstrap 95\% confidence intervals, if applicable. The options are \code{"percentile"} and \code{"normal"}.
#' @param threads                 Integer specifying the number of threads to be used in \code{data.table}. See \code{\link[data.table]{setDTthreads}} for further details.
#' @param model_fits              Logical scalar indicating whether to return the fitted models. Note that if this argument is set to \code{TRUE}, the output of this function may use a lot of memory. The default is \code{FALSE}.
#' @param boot_diag               Logical scalar indicating whether to return the parametric g-formula estimates as well as the coefficients, standard errors, and variance-covariance matrices of the parameters of the fitted models in the bootstrap samples. The default is \code{FALSE}.
#' @param show_progress           Logical scalar indicating whether to print a progress bar for the number of bootstrap samples completed in the R console. This argument is only applicable when \code{parallel} is set to \code{FALSE} and bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0). The default is \code{TRUE}.
#' @param ipw_cutoff_quantile     Percentile by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param ipw_cutoff_value        Cutoff value by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param int_visit_type          Vector of logicals. The kth element is a logical specifying whether to carry forward the intervened value (rather than the natural value) of the treatment variables(s) when performing a carry forward restriction type for the kth intervention in \code{interventions}.
#'                                When the kth element is set to \code{FALSE}, the natural value of the treatment variable(s) in the kth intervention in \code{interventions} will be carried forward.
#'                                By default, this argument is set so that the intervened value of the treatment variable(s) is carried forward for all interventions.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#' @return                        An object of class \code{gformula_survival}. The object is a list with the following components:
#' \item{result}{Results table. For survival outcomes, this contains the estimated risk, risk difference, and risk ratio for all interventions (inculding the natural course) at each time point. For continuous end-of-follow-up outcomes, this contains estimated mean outcome, mean difference, and mean ratio for all interventions (inculding natural course) at the last time point. For binary end-of-follow-up outcomes, this contains the estimated outcome probability, probability difference, and probability ratio for all interventions (inculding natural course) at the last time point. For all outcome types, this also contains the "cumulative percent intervened on" and the "average percent intervened on". If bootstrapping was used, the results table includes the bootstrap risk / mean / probability difference, ratio, standard error, and 95\% confidence interval.}
#' \item{coeffs}{A list of the coefficients of the fitted models.}
#' \item{stderrs}{A list of the standard errors of the coefficients of the fitted models.}
#' \item{vcovs}{A list of the variance-covariance matrices of the parameters of the fitted models.}
#' \item{rmses}{A list of root mean square error (RMSE) values of the fitted models.}
#' \item{hazardratio_val}{Hazard ratio between two interventions (if applicable).}
#' \item{fits}{A list of the fitted models for the time-varying covariates, outcome, and competing event (if applicable). If \code{model_fits} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{sim_data}{A list of data tables of the simulated data. Each element in the list corresponds to one of the interventions. If the argument \code{sim_data_b} is set to \code{FALSE}, a value of \code{NA} is given.}
#' \item{IP_weights}{A numeric vector specifying the inverse probability weights. See "Details".}
#' \item{bootests}{A data.table containing the bootstrap replicates of the parametric g-formula estimates. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootcoeffs}{A list, where the kth element is a list containing the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootstderrs}{A list, where the kth element is a list containing the standard errors of the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootvcovs}{A list, where the kth element is a list containing the variance-covariance matrices of the parameters of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{...}{Some additional elements.}
#'
#' The results for the g-formula simulation are printed with the \code{\link{print.gformula_survival}}, \code{\link{print.gformula_continuous_eof}}, and \code{\link{print.gformula_binary_eof}} functions. To generate graphs comparing the mean estimated covariate values and risks over time and mean observed covariate values and risks over time, use the \code{\link{plot.gformula_survival}}, \code{\link{plot.gformula_continuous_eof}}, and \code{\link{plot.gformula_binary_eof}} functions.
#'
#' @references Chiu YH, Wen L, McGrath S, Logan R, Dahabreh IJ, Hernán MA. Evaluating model specification when using the parametric g-formula in the presence of censoring. American Journal of Epidemiology. In press.
#' @references McGrath S, Lin V, Zhang Z, Petito LC, Logan RW, Hernán MA, and JG Young. gfoRmula: An R package for estimating the effects of sustained treatment strategies via the parametric g-formula. Patterns. 2020;1:100008.
#' @references Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period: application to the healthy worker survivor effect. Mathematical Modelling. 1986;7:1393–1512. [Errata (1987) in Computers and Mathematics with Applications 14, 917.-921. Addendum (1987) in Computers and Mathematics with Applications 14, 923-.945. Errata (1987) to addendum in Computers and Mathematics with Applications 18, 477.].
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
#' gform_basic
#' }
#'
#'
#' ## Estimating the effect of treatment strategies on risk of a failure event
#' ## when competing events exist
#' \donttest{
#' id <- 'id'
#' time_points <- 7
#' time_name <- 't0'
#' covnames <- c('L1', 'L2', 'A')
#' outcome_name <- 'Y'
#' compevent_name <- 'D'
#' outcome_type <- 'survival'
#' covtypes <- c('binary', 'bounded normal', 'binary')
#' histories <- c(lagged, lagavg)
#' histvars <- list(c('A', 'L1', 'L2'), c('L1', 'L2'))
#' covparams <- list(covlink = c('logit', 'identity', 'logit'),
#'                   covmodels = c(L1 ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
#'                                   L3 + as.factor(t0),
#'                                 L2 ~ lag1_A + L1 + lag_cumavg1_L1 +
#'                                   lag_cumavg1_L2 + L3 + as.factor(t0),
#'                                 A ~ lag1_A + L1 + L2 + lag_cumavg1_L1 +
#'                                   lag_cumavg1_L2 + L3 + as.factor(t0)))
#' ymodel <- Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + as.factor(t0)
#' compevent_model <- D ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + as.factor(t0)
#' intvars <- list('A', 'A')
#' interventions <- list(list(c(static, rep(0, time_points))),
#'                       list(c(static, rep(1, time_points))))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_basic <- gformula(obs_data = basicdata, id = id,
#'                         time_points = time_points,
#'                         time_name = time_name, covnames = covnames,
#'                         outcome_name = outcome_name,
#'                         outcome_type = outcome_type,
#'                         compevent_name = compevent_name,
#'                         covtypes = covtypes,
#'                         covparams = covparams, ymodel = ymodel,
#'                         compevent_model = compevent_model,
#'                         intvars = intvars, interventions = interventions,
#'                         int_descript = int_descript,
#'                         histories = histories, histvars = histvars,
#'                         basecovs = c('L3'), nsimul = nsimul,
#'                         seed = 1234)
#' gform_basic
#' }
#'
#'
#' ## Estimating the effect of treatment strategies on the mean of a continuous
#' ## end of follow-up outcome
#' \donttest{
#' library('Hmisc')
#' id <- 'id'
#' time_name <- 't0'
#' covnames <- c('L1', 'L2', 'A')
#' outcome_name <- 'Y'
#' outcome_type <- 'continuous_eof'
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
#' gform_cont_eof <- gformula(obs_data = continuous_eofdata,
#'                            id = id, time_name = time_name,
#'                            covnames = covnames, outcome_name = outcome_name,
#'                            outcome_type = outcome_type, covtypes = covtypes,
#'                            covparams = covparams, ymodel = ymodel,
#'                            intvars = intvars, interventions = interventions,
#'                            int_descript = int_descript,
#'                            histories = histories, histvars = histvars,
#'                            basecovs = c("L3"), nsimul = nsimul, seed = 1234)
#' gform_cont_eof
#' }
#'
#'
#' ## Estimating the effect of threshold interventions on the mean of a binary
#' ## end of follow-up outcome
#' \donttest{
#' outcome_type <- 'binary_eof'
#' id <- 'id_num'
#' time_name <- 'time'
#' covnames <- c('cov1', 'cov2', 'treat')
#' outcome_name <- 'outcome'
#' histories <- c(lagged, cumavg)
#' histvars <- list(c('treat', 'cov1', 'cov2'), c('cov1', 'cov2'))
#' covtypes <- c('binary', 'zero-inflated normal', 'normal')
#' covparams <- list(covmodels = c(cov1 ~ lag1_treat + lag1_cov1 + lag1_cov2 +
#'                                   cov3 + time,
#'                                 cov2 ~ lag1_treat + cov1 + lag1_cov1 +
#'                                   lag1_cov2 + cov3 + time,
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
#' gform_bin_eof <- gformula(obs_data = binary_eofdata,
#'                           outcome_type = outcome_type, id = id,
#'                           time_name = time_name, covnames = covnames,
#'                           outcome_name = outcome_name, covtypes = covtypes,
#'                           covparams = covparams, ymodel = ymodel,
#'                           intvars = intvars, interventions = interventions,
#'                           int_descript = int_descript, histories = histories,
#'                           histvars = histvars, basecovs = c("cov3"),
#'                           seed = 1234, parallel = TRUE, nsamples = 5,
#'                           nsimul = nsimul, ncores = ncores)
#' gform_bin_eof
#' }
#'
#' ## Using IP weighting to estimate natural course risk
#' ## Only the natural course intervention is included for simplicity
#' \donttest{
#' covnames <- c('L', 'A')
#' histories <- c(lagged)
#' histvars <- list(c('A', 'L'))
#' ymodel <- Y ~ L + A
#' covtypes <- c('binary', 'normal')
#' covparams <- list(covmodels = c(L ~ lag1_L + lag1_A,
#'                                 A ~ lag1_L + L + lag1_A))
#' censor_name <- 'C'
#' censor_model <- C ~ L
#' res_censor <- gformula(obs_data = censor_data, id = 'id',
#'                        time_name = 't0', covnames = covnames,
#'                        outcome_name = 'Y', outcome_type = 'survival',
#'                        censor_name = censor_name, censor_model = censor_model,
#'                        covtypes = covtypes,
#'                        covparams = covparams, ymodel = ymodel,
#'                        intvars = NULL, interventions = NULL, int_descript = NULL,
#'                        histories = histories, histvars = histvars,
#'                        seed = 1234)
#' plot(res_censor)
#' }
#'
#' @import data.table
#' @export

gformula <- function(obs_data, id, time_points = NULL,
                     time_name, covnames, covtypes, covparams,
                     covfits_custom = NA, covpredict_custom = NA,
                     histvars = NULL, histories = NA, basecovs = NA,
                     outcome_name, outcome_type, ymodel,
                     compevent_name = NULL, compevent_model = NA,
                     compevent_cens = FALSE,
                     censor_name = NULL, censor_model = NA,
                     intvars = NULL, interventions = NULL,
                     int_times = NULL, int_descript = NULL, ref_int = 0, intcomp = NA,
                     visitprocess = NA, restrictions = NA,
                     yrestrictions = NA, compevent_restrictions = NA,
                     baselags = FALSE,
                     nsimul = NA, sim_data_b = FALSE, seed,
                     nsamples = 0, parallel = FALSE, ncores = NA,
                     ci_method = 'percentile', threads, model_fits = FALSE,
                     boot_diag = FALSE, show_progress = TRUE, ipw_cutoff_quantile = NULL,
                     ipw_cutoff_value = NULL, int_visit_type = NULL, ...){
  if (! outcome_type %in% c('survival', 'continuous_eof', 'binary_eof')){
    stop("outcome_type must be 'survival', 'continuous_eof', or 'binary_eof', but outcome_type was set to", outcome_type)
  }
  if (outcome_type == 'survival'){
    gformula_survival(obs_data = obs_data, id = id, time_points = time_points,
                      time_name = time_name, covnames = covnames,
                      covtypes = covtypes, covparams = covparams,
                      covfits_custom = covfits_custom,
                      covpredict_custom = covpredict_custom,
                      histvars = histvars, histories = histories,
                      basecovs = basecovs, outcome_name = outcome_name,
                      ymodel = ymodel,
                      compevent_name = compevent_name,
                      compevent_model = compevent_model,
                      compevent_cens = compevent_cens,
                      censor_name = censor_name,
                      censor_model = censor_model,
                      intvars = intvars, interventions = interventions,
                      int_times = int_times, int_descript = int_descript,
                      ref_int = ref_int, intcomp = intcomp,
                      visitprocess = visitprocess, restrictions = restrictions,
                      yrestrictions = yrestrictions,
                      compevent_restrictions = compevent_restrictions,
                      baselags = baselags,
                      nsimul = nsimul, sim_data_b = sim_data_b, seed = seed,
                      nsamples = nsamples, parallel = parallel, ncores = ncores,
                      ci_method = ci_method, threads = threads,
                      model_fits = model_fits, boot_diag = boot_diag,
                      show_progress = show_progress, ipw_cutoff_quantile = ipw_cutoff_quantile,
                      ipw_cutoff_value = ipw_cutoff_value, int_visit_type = int_visit_type, ...)
  } else if (outcome_type == 'continuous_eof'){
    gformula_continuous_eof(obs_data = obs_data, id = id,
                            time_name = time_name, covnames = covnames,
                            covtypes = covtypes,
                            covparams = covparams,
                            covfits_custom = covfits_custom,
                            covpredict_custom = covpredict_custom,
                            histvars = histvars,
                            histories = histories, basecovs = basecovs,
                            outcome_name = outcome_name, ymodel = ymodel,
                            censor_name = censor_name,
                            censor_model = censor_model,
                            intvars = intvars, interventions = interventions,
                            int_times = int_times, int_descript = int_descript,
                            ref_int = ref_int,
                            visitprocess = visitprocess,
                            restrictions = restrictions,
                            yrestrictions = yrestrictions, baselags = baselags,
                            nsimul = nsimul, sim_data_b = sim_data_b,
                            seed = seed, nsamples = nsamples,
                            parallel = parallel, ncores = ncores,
                            ci_method = ci_method, threads = threads,
                            model_fits = model_fits, boot_diag = boot_diag,
                            show_progress = show_progress, ipw_cutoff_quantile = ipw_cutoff_quantile,
                            ipw_cutoff_value = ipw_cutoff_value, int_visit_type = int_visit_type, ...)
  } else if (outcome_type == 'binary_eof'){
    gformula_binary_eof(obs_data = obs_data, id = id,
                        time_name = time_name, covnames = covnames,
                        covtypes = covtypes, covparams = covparams,
                        covfits_custom = covfits_custom,
                        covpredict_custom = covpredict_custom,
                        histvars = histvars, histories = histories,
                        basecovs = basecovs, outcome_name = outcome_name,
                        ymodel = ymodel, intvars = intvars,
                        censor_name = censor_name,
                        censor_model = censor_model,
                        interventions = interventions, int_times = int_times,
                        int_descript = int_descript,
                        ref_int = ref_int, visitprocess = visitprocess,
                        restrictions = restrictions,
                        yrestrictions = yrestrictions, baselags = baselags,
                        nsimul = nsimul, sim_data_b = sim_data_b, seed = seed,
                        nsamples = nsamples, parallel = parallel,
                        ncores = ncores,
                        ci_method = ci_method, threads = threads,
                        model_fits = model_fits, boot_diag = boot_diag,
                        show_progress = show_progress, ipw_cutoff_quantile = ipw_cutoff_quantile,
                        ipw_cutoff_value = ipw_cutoff_value, int_visit_type = int_visit_type, ...)
  }
}



#' Estimation of Survival Outcome Under the Parametric G-Formula
#'
#' Based on an observed data set, this internal function estimates the risk over time under multiple
#' user-specified interventions using the parametric g-formula. See McGrath et al. (2020) for
#' further details concerning the application and implementation of the parametric g-formula.
#'
#' To assess model misspecification in the parametric g-formula, users can obtain inverse probability (IP) weighted estimates of the natural course risk and/or means of the time-varying covariates from the observed data.
#' See Chiu et al. (In press) for details.
#' In addition to the general requirements described in McGrath et al. (2020), the requirements for the input data set and the call to the gformula function for such analyses are described below.
#'
#' Users need to include a column in \code{obs_data} with a time-varying censoring variable.
#' Users need to indicate the name of the censoring variable and a model statement for the censoring variable with parameters \code{censor_name} and \code{censor_model}, respectively.
#' When competing events are present, users need to include a column in \code{obs_data} with a time-varying indicator of the competing event variable and need to indicate the name of the competing event variable and the corresponding model statement with parameters \code{compevent_name} and \code{compevent_model}, respectively.
#' Users need to indicate whether to treat competing events as censoring events with the \code{compevent_cens} parameter.
#' Finally, users can specify how to truncate IP weights with the \code{ipw_cutoff_quantile} or \code{ipw_cutoff_value} parameters.
#'
#' In addition to the package output described in McGrath et al. (2020), the output will display estimates of the "cumulative percent intervened on" and the "average percent intervened on". When using a custom intervention function, users need to specify whether each individual at that time point is eligible to contribute person-time to the percent intervened on calculations. Specifically, this must be specified in the \code{eligible_pt} column of \code{newdf}. By default, \code{eligible_pt} is set to \code{TRUE} for each individual at each time point in custom interventions.
#'
#' @param id                      Character string specifying the name of the ID variable in \code{obs_data}.
#' @param time_points             Number of time points to simulate. By default, this argument is set equal to the maximum
#'                                number of records that \code{obs_data} contains for any individual.
#' @param obs_data                Data table containing the observed data.
#' @param seed                    Starting seed for simulations and bootstrapping.
#' @param nsimul                  Number of subjects for whom to simulate data. By default, this argument is set
#'                                equal to the number of subjects in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param compevent_name          Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param censor_name             Character string specifying the name of the censoring variable in \code{obs_data}. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param censor_model            Model statement for the censoring variable. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param intvars                 List, whose elements are vectors of character strings. The kth vector in \code{intvars} specifies the name(s) of the variable(s) to be intervened
#'                                on in each round of the simulation under the kth intervention in \code{interventions}.
#' @param interventions           List, whose elements are lists of vectors. Each list in \code{interventions} specifies a unique intervention on the relevant variable(s) in \code{intvars}. Each vector contains a function
#'                                implementing a particular intervention on a single variable, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#' @param int_times               List, whose elements are lists of vectors. The kth list in \code{int_times} corresponds to the kth intervention in \code{interventions}. Each vector specifies the time points in which the relevant intervention is applied on the corresponding variable in \code{intvars}.
#'                                When an intervention is not applied, the simulated natural course value is used. By default, this argument is set so that all interventions are applied in all time points.
#' @param int_descript            Vector of character strings, each describing an intervention. It must
#'                                be in same order as the entries in \code{interventions}.
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the risk ratio and risk difference. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}. The default is 0.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}. The default is \code{NA}.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}. The default is \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}. These covariates are not simulated using a model but rather carry their value over all time points from the first time point of \code{obs_data}. These covariates should not be included in \code{covnames}. The default is \code{NA}.
#' @param histvars                List of vectors. The kth vector specifies the names of the variables for which the kth history function
#'                                in \code{histories} is to be applied.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}. The default is \code{NA}.
#' @param ymodel                  Model statement for the outcome variable.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#'                                The default is \code{NA}.
#' @param compevent_restrictions  List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the competing event variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the competing event variable takes on the value in the
#'                                second entry. The default is \code{NA}.
#' @param restrictions            List of vectors. Each vector contains as its first entry a covariate for which
#'                                \emph{a priori} knowledge of its distribution is available; its second entry a condition
#'                                under which no knowledge of its distribution is available and that must be \code{TRUE}
#'                                for the distribution of that covariate given that condition to be estimated via a parametric
#'                                model or other fitting procedure; its third entry a function for estimating the distribution
#'                                of that covariate given the condition in the second entry is false such that \emph{a priori} knowledge
#'                                of the covariate distribution is available; and its fourth entry a value used by the function in the
#'                                third entry. The default is \code{NA}.
#' @param visitprocess            List of vectors. Each vector contains as its first entry
#'                                the covariate name of a visit process; its second entry
#'                                the name of a covariate whose modeling depends on the
#'                                visit process; and its third entry the maximum number
#'                                of consecutive visits that can be missed before an
#'                                individual is censored. The default is \code{NA}.
#' @param compevent_model         Model statement for the competing event variable. The default is \code{NA}.
#' @param compevent_cens          Logical scalar indicating whether to treat competing events as censoring events.
#'                                This argument is only applicable for survival outcomes and when a competing even model is supplied (i.e., \code{compevent_name} and \code{compevent_model} are specified).
#'                                If this argument is set to \code{TRUE}, the competing event model will only be used to construct inverse probability weights to estimate the natural course means / risk from the observed data.
#'                                If this argument is set to \code{FALSE}, the competing event model will be used in the parametric g-formula estimates of the risk and will not be used to construct inverse probability weights.
#'                                See "Details". The default is \code{FALSE}.
#' @param intcomp                 List of two numbers indicating a pair of interventions to be compared by a hazard ratio.
#'                                The default is \code{NA}, resulting in no hazard ratio calculation.
#' @param baselags                Logical scalar for specifying the convention used for lagi and lag_cumavgi terms in the model statements when pre-baseline times are not
#'                                included in \code{obs_data} and when the current time index, \eqn{t}, is such that \eqn{t < i}. If this argument is set to \code{FALSE}, the value
#'                                of all lagi and lag_cumavgi terms in this context are set to 0 (for non-categorical covariates) or the reference
#'                                level (for categorical covariates). If this argument is set to \code{TRUE}, the value of lagi and lag_cumavgi terms
#'                                are set to their values at time 0. The default is \code{FALSE}.
#' @param nsamples                Integer specifying the number of bootstrap samples to generate.
#'                                The default is 0.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of CPU cores to use in parallel
#'                                simulation. This argument is required when parallel is set to \code{TRUE}.
#'                                In many applications, users may wish to set this argument equal to \code{parallel::detectCores() - 1}.
#' @param sim_data_b              Logical scalar indicating whether to return the simulated data set. If bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0), this argument must be set to \code{FALSE}. The default is \code{FALSE}.
#' @param ci_method               Character string specifying the method for calculating the bootstrap 95\% confidence intervals, if applicable. The options are \code{"percentile"} and \code{"normal"}.
#' @param threads                 Integer specifying the number of threads to be used in \code{data.table}. See \code{\link[data.table]{setDTthreads}} for further details.
#' @param model_fits              Logical scalar indicating whether to return the fitted models. Note that if this argument is set to \code{TRUE}, the output of this function may use a lot of memory. The default is \code{FALSE}.
#' @param boot_diag               Logical scalar indicating whether to return the parametric g-formula estimates as well as the coefficients, standard errors, and variance-covariance matrices of the parameters of the fitted models in the bootstrap samples. The default is \code{FALSE}.
#' @param show_progress           Logical scalar indicating whether to print a progress bar for the number of bootstrap samples completed in the R console. This argument is only applicable when \code{parallel} is set to \code{FALSE} and bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0). The default is \code{TRUE}.
#' @param ipw_cutoff_quantile     Percentile by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param ipw_cutoff_value        Cutoff value by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param int_visit_type          Vector of logicals. The kth element is a logical specifying whether to carry forward the intervened value (rather than the natural value) of the treatment variables(s) when performing a carry forward restriction type for the kth intervention in \code{interventions}.
#'                                When the kth element is set to \code{FALSE}, the natural value of the treatment variable(s) in the kth intervention in \code{interventions} will be carried forward.
#'                                By default, this argument is set so that the intervened value of the treatment variable(s) is carried forward for all interventions.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#' @return                        An object of class \code{gformula_survival}. The object is a list with the following components:
#' \item{result}{Results table containing the estimated risk and risk ratio for all interventions (inculding the natural course) at each time point as well as the "cumulative percent intervened on" and the "average percent intervened on". If bootstrapping was used, the results table includes the bootstrap mean risk ratio, standard error, and 95\% confidence interval.}
#' \item{coeffs}{A list of the coefficients of the fitted models.}
#' \item{stderrs}{A list of the standard errors of the coefficients of the fitted models.}
#' \item{vcovs}{A list of the variance-covariance matrices of the parameters of the fitted models.}
#' \item{rmses}{A list of root mean square error (RMSE) values of the fitted models.}
#' \item{hazardratio_val}{Hazard ratio between two interventions (if applicable).}
#' \item{fits}{A list of the fitted models for the time-varying covariates, outcome, and competing event (if applicable). If \code{model_fits} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{sim_data}{A list of data tables of the simulated data. Each element in the list corresponds to one of the interventions. If the argument \code{sim_data_b} is set to \code{FALSE}, a value of \code{NA} is given.}
#' \item{IP_weights}{A numeric vector specifying the inverse probability weights. See "Details".}
#' \item{bootests}{A data.table containing the bootstrap replicates of the parametric g-formula estimates. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootcoeffs}{A list, where the kth element is a list containing the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootstderrs}{A list, where the kth element is a list containing the standard errors of the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootvcovs}{A list, where the kth element is a list containing the variance-covariance matrices of the parameters of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{...}{Some additional elements.}
#'
#' The results for the g-formula simulation under various interventions only for the first and last time points are printed with the \code{\link{print.gformula_survival}} function. To generate graphs comparing the mean estimated covariate values and risks over time and mean observed covariate values and risks over time, use the \code{\link{plot.gformula_survival}} function.
#'
#' @seealso \code{\link{gformula}}
#' @references Chiu YH, Wen L, McGrath S, Logan R, Dahabreh IJ, Hernán MA. Evaluating model specification when using the parametric g-formula in the presence of censoring. American Journal of Epidemiology. In press.
#' @references McGrath S, Lin V, Zhang Z, Petito LC, Logan RW, Hernán MA, and JG Young. gfoRmula: An R package for estimating the effects of sustained treatment strategies via the parametric g-formula. Patterns. 2020;1:100008.
#' @references Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period: application to the healthy worker survivor effect. Mathematical Modelling. 1986;7:1393–1512. [Errata (1987) in Computers and Mathematics with Applications 14, 917.-921. Addendum (1987) in Computers and Mathematics with Applications 14, 923-.945. Errata (1987) to addendum in Computers and Mathematics with Applications 18, 477.].
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
#'
#' ## Estimating the effect of treatment strategies on risk of a failure event
#' ## when competing events exist
#' \donttest{
#' id <- 'id'
#' time_points <- 7
#' time_name <- 't0'
#' covnames <- c('L1', 'L2', 'A')
#' outcome_name <- 'Y'
#' compevent_name <- 'D'
#' covtypes <- c('binary', 'bounded normal', 'binary')
#' histories <- c(lagged, lagavg)
#' histvars <- list(c('A', 'L1', 'L2'), c('L1', 'L2'))
#' covparams <- list(covlink = c('logit', 'identity', 'logit'),
#'                   covmodels = c(L1 ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
#'                                   L3 + as.factor(t0),
#'                                 L2 ~ lag1_A + L1 + lag_cumavg1_L1 +
#'                                   lag_cumavg1_L2 + L3 + as.factor(t0),
#'                                 A ~ lag1_A + L1 + L2 + lag_cumavg1_L1 +
#'                                   lag_cumavg1_L2 + L3 + as.factor(t0)))
#' ymodel <- Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + as.factor(t0)
#' compevent_model <- D ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + as.factor(t0)
#' intvars <- list('A', 'A')
#' interventions <- list(list(c(static, rep(0, time_points))),
#'                       list(c(static, rep(1, time_points))))
#' int_descript <- c('Never treat', 'Always treat')
#' nsimul <- 10000
#'
#' gform_basic <- gformula_survival(obs_data = basicdata, id = id,
#'                                  time_points = time_points,
#'                                  time_name = time_name, covnames = covnames,
#'                                  outcome_name = outcome_name,
#'                                  compevent_name = compevent_name,
#'                                  covtypes = covtypes,
#'                                  covparams = covparams, ymodel = ymodel,
#'                                  compevent_model = compevent_model,
#'                                  intvars = intvars, interventions = interventions,
#'                                  int_descript = int_descript,
#'                                  histories = histories, histvars = histvars,
#'                                  basecovs = c('L3'), nsimul = nsimul,
#'                                  seed = 1234)
#' gform_basic
#' }
#'
#' ## Using IP weighting to estimate natural course risk
#' ## Only the natural course intervention is included for simplicity
#' \donttest{
#' covnames <- c('L', 'A')
#' histories <- c(lagged)
#' histvars <- list(c('A', 'L'))
#' ymodel <- Y ~ L + A
#' covtypes <- c('binary', 'normal')
#' covparams <- list(covmodels = c(L ~ lag1_L + lag1_A,
#'                                 A ~ lag1_L + L + lag1_A))
#' censor_name <- 'C'
#' censor_model <- C ~ L
#' res_censor <- gformula(obs_data = censor_data, id = 'id',
#'                        time_name = 't0', covnames = covnames,
#'                        outcome_name = 'Y', outcome_type = 'survival',
#'                        censor_name = censor_name, censor_model = censor_model,
#'                        covtypes = covtypes,
#'                        covparams = covparams, ymodel = ymodel,
#'                        intvars = NULL, interventions = NULL, int_descript = NULL,
#'                        histories = histories, histvars = histvars,
#'                        seed = 1234)
#' plot(res_censor)
#' }
#'
#' @import data.table
#' @export
gformula_survival <- function(obs_data, id, time_points = NULL,
                              time_name, covnames, covtypes, covparams,
                              covfits_custom = NA, covpredict_custom = NA,
                              histvars = NULL, histories = NA, basecovs = NA,
                              outcome_name, ymodel,
                              compevent_name = NULL, compevent_model = NA,
                              compevent_cens = FALSE,
                              censor_name = NULL, censor_model = NA,
                              intvars = NULL, interventions = NULL,
                              int_times = NULL, int_descript = NULL, ref_int = 0, intcomp = NA,
                              visitprocess = NA, restrictions = NA,
                              yrestrictions = NA, compevent_restrictions = NA,
                              baselags = FALSE,
                              nsimul = NA, sim_data_b = FALSE, seed,
                              nsamples = 0, parallel = FALSE, ncores = NA,
                              ci_method = 'percentile', threads,
                              model_fits = FALSE, boot_diag = FALSE,
                              show_progress = TRUE, ipw_cutoff_quantile = NULL,
                              ipw_cutoff_value = NULL, int_visit_type = NULL, ...){

  lag_indicator <- lagavg_indicator <- cumavg_indicator <- c()
  lag_indicator <- update_lag_indicator(covparams$covmodels, lag_indicator)
  lagavg_indicator <- update_lagavg_indicator(covparams$covmodels, lagavg_indicator)
  cumavg_indicator <- update_cumavg_indicator(covparams$covmodels, cumavg_indicator)

  comprisk <- !(length(compevent_model) == 1 && is.na(compevent_model))
  censor <- !(length(censor_model) == 1 && is.na(censor_model))

  if (!missing(ymodel)){
    lag_indicator <- update_lag_indicator(ymodel, lag_indicator)
    lagavg_indicator <- update_lagavg_indicator(ymodel, lagavg_indicator)
    cumavg_indicator <- update_cumavg_indicator(ymodel, cumavg_indicator)
  }
  if (comprisk){
    lag_indicator <- update_lag_indicator(compevent_model, lag_indicator)
    lagavg_indicator <- update_lagavg_indicator(compevent_model, lagavg_indicator)
    cumavg_indicator <- update_cumavg_indicator(compevent_model, cumavg_indicator)
  }
  if (censor){
    lag_indicator <- update_lag_indicator(censor_model, lag_indicator)
    lagavg_indicator <- update_lagavg_indicator(censor_model, lagavg_indicator)
    cumavg_indicator <- update_cumavg_indicator(censor_model, cumavg_indicator)
  }
  histvals <- list(lag_indicator = lag_indicator, lagavg_indicator = lagavg_indicator,
                   cumavg_indicator = cumavg_indicator)


  if (!missing(threads)){
    setDTthreads(threads = threads)
  }
  else {
    threads <- getDTthreads()
  }
  outcome_type <- 'survival'
  hazardratio <- !(length(intcomp) == 1 && is.na(intcomp))

  error_catch(id = id, nsimul = nsimul, intvars = intvars, interventions = interventions,
              int_times = int_times, int_descript = int_descript,
              covnames = covnames, covtypes = covtypes, basecovs = basecovs,
              histvars = histvars, histories = histories, compevent_model = compevent_model,
              hazardratio = hazardratio, intcomp = intcomp, time_points = time_points,
              outcome_type = outcome_type, time_name = time_name,
              obs_data = obs_data, parallel = parallel, ncores = ncores,
              nsamples = nsamples, sim_data_b = sim_data_b,
              outcome_name = outcome_name, compevent_name = compevent_name,
              comprisk = comprisk, censor = censor, censor_name = censor_name,
              covmodels = covparams$covmodels,
              histvals = histvals, ipw_cutoff_quantile = ipw_cutoff_quantile,
              ipw_cutoff_value = ipw_cutoff_value)


  if (comprisk & compevent_cens){
    comprisk2 <- TRUE; compevent2_name <- compevent_name; compevent2_model <- compevent_model
    comprisk <- FALSE; compevent_name <- NULL; compevent_model <- NA
  } else {
    comprisk2 <- FALSE; compevent2_name <- NULL; compevent2_model <- NA
  }

  min_time <- min(obs_data[[time_name]])
  below_zero_indicator <- min_time < 0


  obs_data <- copy(obs_data)



  max_visits <- NA
  if (!is.na(visitprocess[[1]][[1]])){
    for (vp in visitprocess){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(vp[1], paste("lag1_ts_", vp[1], "!=", vp[3], sep = ""),
                               ### rwl paste("visit_sum_", vp[3], "_", vp[1], "!=0", sep = ""),
                               simple_restriction, 1),
                             c(vp[2], paste(vp[1], "==1", sep = ""), carry_forward)))
      if (is.na(max_visits[1])){
        max_visits <- as.numeric(vp[3])
      } else {
        max_visits <- c(max_visits, as.numeric(vp[3]))
      }
      if (is.na(histories[1])){
        histories <- c(visit_sum)
      } else {
        histories <- c(visit_sum, histories)
        histvars <- append(list(c(vp[1])),histvars)

      }
    }
  }


  for (t in 0:max(obs_data[[time_name]])) {
    make_histories(pool = obs_data, histvars = histvars, histvals = histvals,
                   histories = histories, time_name = time_name, t = t, id = id ,
                   max_visits = max_visits, baselags = baselags,
                   below_zero_indicator = below_zero_indicator)
  }

  sample_size <- length(unique(obs_data[[id]]))
  if (is.null(time_points)){
    time_points <- max(obs_data[[time_name]])+1
  }


  for (i in seq_along(covnames)){
    if (covtypes[i] == 'absorbing'){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(covnames[i], paste("lag1_", covnames[i], "==0", sep = ""),
                               carry_forward, 1)))
      covtypes[i] <- 'binary'
    }
  }


  # Create 1-indexed numerical IDs for observed datasets
  ids <- as.data.table(sort(unique(obs_data[[id]])))
  ids[, 'newid' := seq_len(.N)]
  setkeyv(obs_data, id)
  obs_data <- obs_data[J(ids), allow.cartesian = TRUE]
  obs_data_geq_0 <- obs_data[obs_data[[time_name]] >= 0]

  # Set default number of simulated individuals to equal number of individuals in
  # observed dataset
  if (is.na(nsimul)){
    nsimul <- length(unique(obs_data$newid))
  }


  # Generate seeds for simulations and bootstrapping
  set.seed(seed)
  newseeds <- sample.int(2^30, size = nsamples + 1)
  subseed <- newseeds[1]
  bootseeds <- newseeds[2:(nsamples + 1)]

  # Determine ranges of observed covariates and outcome
  ranges <- lapply(seq_along(covnames), FUN = function(i){
    if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
        covtypes[i] == 'truncated normal') {
      range(obs_data_geq_0[[covnames[i]]])
    } else if (covtypes[i] == 'zero-inflated normal'){
      range(obs_data_geq_0[obs_data_geq_0[[covnames[i]]] > 0][[covnames[i]]])
    } else {
      NA
    }
  })

  # Fit models to covariates and outcome variable
  if (time_points > 1){
    fitcov <- pred_fun_cov(covparams = covparams, covnames = covnames, covtypes = covtypes,
                           covfits_custom = covfits_custom, restrictions = restrictions,
                           time_name = time_name, obs_data = obs_data_geq_0,
                           model_fits = model_fits)
    names(fitcov) <- covnames
  } else {
    fitcov <- NULL
  }
  fitY <- pred_fun_Y(ymodel, yrestrictions, outcome_type, outcome_name, time_name, obs_data_geq_0,
                     model_fits = model_fits)

  # If competing event exists, fit model for competing event variable
  if (comprisk){
    fitD <- pred_fun_D(compevent_model, compevent_restrictions, obs_data_geq_0,
                       model_fits = model_fits)
  } else {
    fitD <- NA
  }
  if (comprisk2){
    fitD2 <- pred_fun_D(compevent2_model, NA, obs_data_geq_0,
                        model_fits = model_fits)
  } else {
    fitD2 <- NA
  }
  if (censor){
    fitC <- pred_fun_D(censor_model, NA, obs_data_geq_0, model_fits = model_fits)
  } else {
    fitC <- NA
  }


  obs_data_noresample <- copy(obs_data)
  len <- length(unique(obs_data$newid))
  # If the number of user desired simulations differs from the number of individuals in
  # the observed dataset, sample the desired number of observed IDs with replacement
  if (nsimul < len){
    ids <- as.data.table(sort(sample(unique(obs_data$newid), nsimul, replace = TRUE)))
    colnames(ids) <- "newid"
    ids[, 'sid' := seq_len(.N)]
    obs_data <- merge(ids, obs_data, all.x = TRUE, by = "newid")
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  } else if (nsimul > len){
    ids <- as.data.table(sample(unique(obs_data$newid), nsimul, replace = TRUE))
    ids[, 'newid' := 1:nsimul]
    colnames(ids) <- c("newid", "sid")
    setkeyv(obs_data, "newid")
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE]
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  }

  # Add natural course to list of interventions
  if (!is.null(interventions)){
    comb_interventions <- c(list(list(c(natural))), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    comb_interventions <- list(list(c(natural)))
    comb_intvars <- list('none')
  }

  if (is.null(int_times)){
    comb_int_times <- list()
    for (i in seq_along(comb_interventions)){
      comb_int_times[[i]] <- lapply(seq_along(comb_interventions[[i]]),
                                    FUN = function(i) {0:(time_points - 1)})
    }
  } else {
    comb_int_times <- c(list(list(0:(time_points - 1))), int_times)
  }

  if (is.null(int_visit_type)){
    int_visit_type <- rep(TRUE, length(comb_interventions))
  } else {
    int_visit_type <- c(TRUE, int_visit_type)
  }

  # Simulate pooled-over-time datasets containing covariates, outcome, and risk for each
  # subject
  if (parallel){
    cl <- prep_cluster(ncores = ncores, threads = threads , covtypes = covtypes)
    pools <- parallel::parLapply(cl, seq_along(comb_interventions), simulate,
                                 fitcov = fitcov, fitY = fitY, fitD = fitD,
                                 yrestrictions = yrestrictions,
                                 compevent_restrictions = compevent_restrictions,
                                 restrictions = restrictions,
                                 outcome_name = outcome_name, compevent_name = compevent_name,
                                 time_name = time_name,
                                 intvars = comb_intvars, interventions = comb_interventions,
                                 int_times = comb_int_times, histvars = histvars,
                                 histvals = histvals, histories = histories,
                                 covparams = covparams, covnames = covnames, covtypes = covtypes,
                                 covpredict_custom = covpredict_custom, basecovs = basecovs,
                                 comprisk = comprisk, ranges = ranges,
                                 outcome_type = outcome_type,
                                 subseed = subseed, time_points = time_points,
                                 obs_data = obs_data, parallel = parallel, max_visits = max_visits,
                                 baselags = baselags, below_zero_indicator = below_zero_indicator,
                                 min_time = min_time, show_progress = FALSE, int_visit_type = int_visit_type, ...)
    parallel::stopCluster(cl)
  } else {
    pools <- lapply(seq_along(comb_interventions), FUN = function(i){
      simulate(fitcov = fitcov, fitY = fitY, fitD = fitD,
               yrestrictions = yrestrictions,
               compevent_restrictions = compevent_restrictions,
               restrictions = restrictions,
               outcome_name = outcome_name, compevent_name = compevent_name,
               time_name = time_name,
               intvars = comb_intvars[[i]], interventions = comb_interventions[[i]],
               int_times = comb_int_times[[i]], histvars = histvars, histvals = histvals,
               histories = histories, covparams = covparams,
               covnames = covnames, covtypes = covtypes,
               covpredict_custom = covpredict_custom, basecovs = basecovs, comprisk = comprisk,
               ranges = ranges,
               outcome_type = outcome_type,
               subseed = subseed, time_points = time_points,
               obs_data = obs_data, parallel = parallel, max_visits = max_visits,
               baselags = baselags, below_zero_indicator = below_zero_indicator,
               min_time = min_time, show_progress = FALSE, int_visit_type = int_visit_type[i], ...)
    })
  }

  nat_pool <- pools[[1]] # Natural course data
  pools <- pools[-1] # List of intervention datasets

  # Initialize results matrices
  result_ratio <- result_diff <- int_result <-
    matrix(NA, nrow = length(pools) + 1, ncol = time_points)

  # Calculate mean risk over all subjects at each time for natural course
  nat_result <- tapply(nat_pool$poprisk, nat_pool[[time_name]], FUN = mean)

  if (ref_int == 0){
    # Set reference intervention to the natural course
    ref_result <- nat_result
  } else {
    # Set reference intervention as specified
    # Calculate mean risk over all subjects at each time for this intervention
    ref_result <- tapply(pools[[ref_int]]$poprisk, pools[[ref_int]][[time_name]], FUN = mean)
  }

  # Compile results
  int_result[1, ] <- nat_result
  result_ratio[1, ] <- int_result[1, ]/ref_result
  result_diff[1, ] <- int_result[1, ] - ref_result
  # Calculate mean risk over all subjects at each time for all interventions other than
  # the natural course
  if (length(comb_interventions) > 1){
    for (i in 2:(length(pools) + 1)){
      int_result[i, ] <- tapply(pools[[i - 1]]$poprisk, pools[[i - 1]][[time_name]],
                                FUN = mean)
      result_ratio[i, ] <- int_result[i, ]/ref_result
      result_diff[i, ] <- int_result[i, ] - ref_result
    }
  }

  # Get hazard ratio results
  if (hazardratio){
    # Generate dataset containing failure/censor time information for each subject
    # under each intervention
    pools_hr <- lapply(seq_along(intcomp), FUN = hr_helper, intcomp = intcomp,
                       time_name = time_name, pools = pools)
    data_hr <- rbindlist(pools_hr)
    names(data_hr)[names(data_hr) == time_name] <- "t0"
    names(data_hr)[names(data_hr) == outcome_name] <- "Y"
    # Factor event variable
    data_hr$event <- factor(data_hr$Ycomp, 0:2, labels=c("censor", "Y", "D"))

    if (comprisk){
      hr_data <- survival::finegray(survival::Surv(t0, event) ~ ., data = data_hr, etype = "Y")
      hr_res <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ regime, data = hr_data)
      hr_res <- exp(hr_res$coefficients)
    }
    else {
      hr_res <- survival::coxph(formula = survival::Surv(t0, Y == "1") ~ regime, data = data_hr)
      hr_res <- exp(hr_res$coefficients)
    }
    names(hr_res) <- "Est. HR"
  } else {
    hr_res <- NA
  }

  # Calculate percent intervened
  percent_intervened_res <- get_percent_intervened(pools = pools)

  # Calculate user specified number of bootstrap risk ratios
  if (nsamples > 0){
    if (parallel){
      cl <- prep_cluster(ncores = ncores, threads = threads , covtypes = covtypes,
                         bootstrap_option = TRUE)
      final_bs <- parallel::parLapply(cl, 1:nsamples, bootstrap_helper_with_trycatch, time_points = time_points,
                                      obs_data = obs_data_noresample, bootseeds = bootseeds,
                                      intvars = comb_intvars, interventions = comb_interventions, int_times = comb_int_times, ref_int = ref_int,
                                      covparams = covparams, covnames = covnames, covtypes = covtypes,
                                      covfits_custom = covfits_custom, covpredict_custom = covpredict_custom,
                                      basecovs = basecovs, ymodel = ymodel,
                                      histvars = histvars, histvals = histvals, histories = histories,
                                      comprisk = comprisk, compevent_model = compevent_model,
                                      yrestrictions = yrestrictions,
                                      compevent_restrictions = compevent_restrictions,
                                      restrictions = restrictions, outcome_type = outcome_type,
                                      ranges = ranges,
                                      time_name = time_name, outcome_name = outcome_name,
                                      compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                                      max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp,
                                      boot_diag = boot_diag, nsimul = nsimul, baselags = baselags,
                                      below_zero_indicator = below_zero_indicator, min_time = min_time,
                                      show_progress = FALSE, int_visit_type = int_visit_type, ...)
      parallel::stopCluster(cl)
    } else {
      if (show_progress){
        pb <- progress::progress_bar$new(total = nsamples * length(comb_interventions),
                                         clear = FALSE,
                                         format = 'Bootstrap progress [:bar] :percent, Elapsed time :elapsed, Est. time remaining :eta')
      }
      final_bs <- lapply(1:nsamples, FUN = bootstrap_helper_with_trycatch, time_points = time_points,
                         obs_data = obs_data_noresample, bootseeds = bootseeds,
                         intvars = comb_intvars, interventions = comb_interventions, int_times = comb_int_times, ref_int = ref_int,
                         covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, covpredict_custom = covpredict_custom,
                         basecovs = basecovs, ymodel = ymodel,
                         histvars = histvars, histvals = histvals, histories = histories,
                         comprisk = comprisk, compevent_model = compevent_model,
                         yrestrictions = yrestrictions,
                         compevent_restrictions = compevent_restrictions,
                         restrictions = restrictions,
                         outcome_type = outcome_type,
                         ranges = ranges,
                         time_name = time_name, outcome_name = outcome_name,
                         compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                         max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp,
                         boot_diag = boot_diag, nsimul = nsimul, baselags = baselags,
                         below_zero_indicator = below_zero_indicator, min_time = min_time,
                         show_progress = show_progress, pb = pb, int_visit_type = int_visit_type, ...)
    }

    comb_result <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$Result))
    }))
    comb_RR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultRatio))
    }))
    comb_RD <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultDiff))
    }))
    comb_HR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(m$ResultHR)
    }))

    comb_result$t0 <- comb_RR$t0 <- comb_RD$t0 <-
      rep(0:(time_points - 1), nsamples)

    se_result <- comb_result[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]
    se_RR <- comb_RR[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]
    se_RD <- comb_RD[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]
    if (hazardratio){
      hr_res[2] <- stats::sd(comb_HR$V1, na.rm = TRUE)
    }

    if (ci_method == 'normal'){
      ci_lb_result <- t(int_result) - stats::qnorm(0.975)*se_result[,-c('t0')]
      ci_lb_RR <- t(result_ratio) - stats::qnorm(0.975)*se_RR[,-c('t0')]
      ci_lb_RD <- t(result_diff) - stats::qnorm(0.975)*se_RD[,-c('t0')]
      ci_ub_result <- t(int_result) + stats::qnorm(0.975)*se_result[,-c('t0')]
      ci_ub_RR <- t(result_ratio) + stats::qnorm(0.975)*se_RR[,-c('t0')]
      ci_ub_RD <- t(result_diff) + stats::qnorm(0.975)*se_RD[,-c('t0')]
      if (hazardratio){
        hr_res[3:4] <- c(hr_res[1] - stats::qnorm(0.975)*hr_res[2],
                         hr_res[1] + stats::qnorm(0.975)*hr_res[2])
      }
    }
    if (ci_method == 'percentile') {
      ci_lb_result <- comb_result[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_lb_RR <- comb_RR[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_lb_RD <- comb_RD[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_ub_result <- comb_result[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
      ci_ub_RR <- comb_RR[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
      ci_ub_RD <- comb_RD[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
      if (hazardratio){
        hr_res[3:4] <- stats::quantile(comb_HR$V1, probs = c(0.025, 0.975), na.rm = TRUE)
      }
    }
    if (hazardratio){
      names(hr_res)[2:4] <- c('HR SE', 'HR lower 95% CI', 'HR upper 95% CI')
    }
  }
  if (nsamples > 0 & boot_diag){
    bootests <- comb_result
    if (!is.null(int_descript)){
      colnames(bootests)[1:(1 + length(interventions))] <- c('Natural course', int_descript)
    } else {
      colnames(bootests)[1:length(comb_interventions)] <- c('Natural course', paste('Intervention', 1:length(interventions)))
    }
    bootests[, 'Bootstrap replicate'] <- rep(1:nsamples, each = time_points)
    bootcoeffs <- lapply(final_bs, "[[", 'bootcoeffs')
    bootstderrs <- lapply(final_bs, "[[", 'bootstderrs')
    bootvcovs <- lapply(final_bs, "[[", 'bootvcovs')
  } else {
    bootests <- NULL
    bootcoeffs <- NULL
    bootstderrs <- NULL
    bootvcovs <- NULL
  }

  plot_info <- get_plot_info(outcome_name = outcome_name,
                             compevent_name = compevent_name,
                             compevent2_name = compevent2_name,
                             censor_name = censor_name,
                             time_name = time_name,
                             id = id,
                             time_points = time_points,
                             covnames = covnames,
                             covtypes = covtypes,
                             nat_pool = nat_pool,
                             nat_result = nat_result,
                             comprisk = comprisk,
                             comprisk2 = comprisk2,
                             censor = censor,
                             fitD2 = fitD2,
                             fitC = fitC,
                             outcome_type = outcome_type,
                             obs_data = obs_data_noresample,
                             ipw_cutoff_quantile = ipw_cutoff_quantile,
                             ipw_cutoff_value = ipw_cutoff_value)
  obs_results <- plot_info$obs_results

  # Generate results table
  if (!is.null(interventions)){
    resultdf <- lapply(1:time_points, function(i){
      if (nsamples > 0){
        rowdfs <- lapply(1:dim(int_result)[1], function(k){
          data.table(t = i - 1, Intervention = k - 1, Risk = int_result[k, ][i],
                     Risk_SE = se_result[[paste0('V',k)]][i],
                     Risk_CI_LL95 = ci_lb_result[[paste0('V',k)]][i],
                     Risk_CI_UL95 = ci_ub_result[[paste0('V',k)]][i],
                     RiskRatio = result_ratio[k, ][i],
                     RR_SE = se_RR[[paste0('V',k)]][i],
                     RR_CI_LL95 = ci_lb_RR[[paste0('V',k)]][i],
                     RR_CI_UL95 = ci_ub_RR[[paste0('V',k)]][i],
                     RiskDiff = result_diff[k, ][i],
                     RD_SE = se_RD[[paste0('V',k)]][i],
                     RD_CI_LL95 = ci_lb_RD[[paste0('V',k)]][i],
                     RD_CI_UL95 = ci_ub_RD[[paste0('V',k)]][i])

        })
      } else {
        rowdfs <- lapply(1:dim(int_result)[1], function(k){
          data.table(t = i - 1, Intervention = k - 1, Risk = int_result[k, ][i],
                     RiskRatio = result_ratio[k, ][i], RiskDiff = result_diff[k, ][i])
        })
      }
      rowdfs <- rbindlist(rowdfs)
      rowdfs[, 'obs_risk' := c(obs_results[[2]][i], rep(NA, dim(rowdfs)[1] - 1))]
      return (rowdfs)
    })
    resultdf <- rbindlist(resultdf)
  } else {
    if (nsamples > 0){
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             Risk = t(int_result),
                             Risk_SE = se_RR[,-c('t0')], Risk_CI_LL95 = ci_lb_result[,-c('t0')],
                             Risk_CI_UL95 = ci_ub_result[,-c('t0')], RiskRatio = t(result_ratio),
                             RR_SE = se_RR[,-c('t0')],
                             RR_CI_LL95 = ci_lb_RR[,-c('t0')], RR_CI_UL95 = ci_ub_RR[,-c('t0')],
                             RiskDiff = t(result_diff),
                             RD_SE = se_RD[,-c('t0')],
                             RD_CI_LL95 = ci_lb_RD[,-c('t0')], RD_CI_UL95 = ci_ub_RD[,-c('t0')])
    } else {
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             Risk = t(int_result), RiskRatio = t(result_ratio),
                             RiskDiff = t(result_diff))
    }
    resultdf[, 'obs_risk' := obs_results[[2]]]
  }
  obs_risk_name <- ifelse(censor, 'IP weighted risk', 'NP Risk')
  if (nsamples > 0){
    colnames(resultdf) <- c("k", "Interv.", "g-form risk",
                            "Risk SE", "Risk lower 95% CI",
                            "Risk upper 95% CI", "Risk ratio",
                            "RR SE",
                            "RR lower 95% CI", "RR upper 95% CI",
                            "Risk difference", "RD SE",
                            "RD lower 95% CI", "RD upper 95% CI",
                            obs_risk_name)
    setcolorder(resultdf, c("k", "Interv.", obs_risk_name, "g-form risk",
                            "Risk SE",
                            "Risk lower 95% CI", "Risk upper 95% CI",
                            "Risk ratio", "RR SE",
                            "RR lower 95% CI", "RR upper 95% CI",
                            "Risk difference", "RD SE",
                            "RD lower 95% CI", "RD upper 95% CI"))
  } else {
    colnames(resultdf) <- c("k", "Interv.", "g-form risk", "Risk ratio",
                            "Risk difference", obs_risk_name)
    setcolorder(resultdf, c("k", "Interv.", obs_risk_name, "g-form risk",
                            "Risk ratio", "Risk difference"))
  }
  resultdf[k == max(k), '% Intervened On'] <- percent_intervened_res$percent_intervened
  resultdf[k == max(k), 'Aver % Intervened On'] <- percent_intervened_res$average_percent_intervened

  if (time_points > 1){
    fits <- fitcov
    fits[[length(fits) + 1]] <- fitY
    names(fits)[length(fits)] <- outcome_name
  } else {
    fits <- list(fitY)
  }
  if (!is.na(fitD)[[1]]){
    fits[[length(fits) + 1]] <- fitD
    names(fits)[length(fits)] <- compevent_name
  }
  if (!is.na(fitC)[[1]]){
    fits[[length(fits) + 1]] <- fitC
    names(fits)[length(fits)] <- censor_name
  }
  if (!is.na(fitD2)[[1]]){
    fits[[length(fits) + 1]] <- fitD2
    names(fits)[length(fits)] <- compevent2_name
  }

  # Add list of coefficients for covariates, outcome variable, and competing event
  # variable (if any) to results output
  coeffs <- get_coeffs(fits = fits, fitD = fitD, time_points = time_points,
                       outcome_name = outcome_name, compevent_name = compevent_name,
                       covnames = covnames)
  stderrs <- get_stderrs(fits = fits, fitD = fitD, time_points = time_points,
                         outcome_name = outcome_name, compevent_name = compevent_name,
                         covnames = covnames)
  vcovs <- get_vcovs(fits = fits, fitD = fitD, time_points = time_points,
                     outcome_name = outcome_name, compevent_name = compevent_name,
                     covnames = covnames)

  rmses <- lapply(seq_along(fits), FUN = rmse_calculate, fits = fits, covnames = covnames,
                  covtypes = covtypes)
  if (!is.na(fitD)[[1]]){
    if (time_points == 1){
      rmses <- stats::setNames(rmses, c(outcome_name, compevent_name))
    } else {
      rmses <- stats::setNames(rmses, c(covnames, outcome_name, compevent_name))
    }
  }
  else {
    if (time_points == 1){
      rmses <- stats::setNames(rmses, outcome_name)
    } else {
      rmses <- stats::setNames(rmses, c(covnames, outcome_name))
    }
  }

  # Create header
  header <- get_header(int_descript, sample_size, nsimul, nsamples, ref_int)

  if (sim_data_b){
    sim_data <- c(list('Natural course' = nat_pool), pools)
    if (!is.null(int_descript)){
      names(sim_data)[2:length(sim_data)] <- int_descript
    }
  } else {
    sim_data <- NA
  }
  if (!model_fits){
    fits <- NULL
  }

  res <- list(
    result = resultdf,
    coeffs = coeffs,
    stderrs = stderrs,
    vcovs = vcovs,
    rmses = rmses,
    hazardratio_val = hr_res,
    fits = fits,
    sim_data = sim_data,
    IP_weights = obs_results$w,
    bootests = bootests,
    bootcoeffs = bootcoeffs,
    bootstderrs = bootstderrs,
    bootvcovs = bootvcovs,
    time_name = time_name,
    time_points = time_points,
    covnames = covnames,
    covtypes = covtypes,
    dt_cov_plot = plot_info$dt_cov_plot,
    dt_out_plot = plot_info$dt_out_plot,
    nsamples = nsamples,
    interventions = interventions,
    comprisk = comprisk,
    header = header
  )
  class(res) <- c("gformula_survival", "gformula")
  return (res)
}

#' Estimation of Continuous End-of-Follow-Up Outcome Under the Parametric G-Formula
#'
#' Based on an observed data set, this internal function estimates the outcome mean at end-of-follow-up under
#' multiple user-specified interventions using the parametric g-formula. See McGrath et al. (2020) for
#' further details concerning the application and implementation of the parametric g-formula.
#'
#' To assess model misspecification in the parametric g-formula, users can obtain inverse probability (IP) weighted estimates of the natural course means of the time-varying covariates from the observed data.
#' See Chiu et al. (In press) for details.
#' In addition to the general requirements described in McGrath et al. (2020), the requirements for the input data set and the call to the gformula function for such analyses are described below.
#'
#' Users need to include a column in \code{obs_data} with a time-varying censoring variable.
#' Users need to indicate the name of the censoring variable and a model statement for the censoring variable with parameters \code{censor_name} and \code{censor_model}, respectively.
#' Finally, users can specify how to truncate IP weights with the \code{ipw_cutoff_quantile} or \code{ipw_cutoff_value} parameters.
#'
#' In addition to the package output described in McGrath et al. (2020), the output will display estimates of the "cumulative percent intervened on" and the "average percent intervened on". When using a custom intervention function, users need to specify whether each individual at that time point is eligible to contribute person-time to the percent intervened on calculations. Specifically, this must be specified in the \code{eligible_pt} column of \code{newdf}. By default, \code{eligible_pt} is set to \code{TRUE} for each individual at each time point in custom interventions.
#'
#' @param id                      Character string specifying the name of the ID variable in \code{obs_data}.
#' @param obs_data                Data table containing the observed data.
#' @param seed                    Starting seed for simulations and bootstrapping.
#' @param nsimul                  Number of subjects for whom to simulate data. By default, this argument is set
#'                                equal to the number of subjects in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param censor_name             Character string specifying the name of the censoring variable in \code{obs_data}. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param censor_model            Model statement for the censoring variable. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param intvars                 List, whose elements are vectors of character strings. The kth vector in \code{intvars} specifies the name(s) of the variable(s) to be intervened
#'                                on in each round of the simulation under the kth intervention in \code{interventions}.
#' @param interventions           List, whose elements are lists of vectors. Each list in \code{interventions} specifies a unique intervention on the relevant variable(s) in \code{intvars}. Each vector contains a function
#'                                implementing a particular intervention on a single variable, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#' @param int_times               List, whose elements are lists of vectors. The kth list in \code{int_times} corresponds to the kth intervention in \code{interventions}. Each vector specifies the time points in which the relevant intervention is applied on the corresponding variable in \code{intvars}.
#'                                When an intervention is not applied, the simulated natural course value is used. By default, this argument is set so that all interventions are applied in all time points.
#' @param int_descript            Vector of character strings, each describing an intervention. It must
#'                                be in same order as the entries in \code{interventions}.
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the end-of-follow-up mean ratio and mean difference. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}. The default is 0.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}. The default is \code{NA}.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}. The default is \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}. These covariates are not simulated using a model but rather carry their value over all time points from the first time point of \code{obs_data}. These covariates should not be included in \code{covnames}. The default is \code{NA}.
#' @param histvars                List of vectors. The kth vector specifies the names of the variables for which the kth history function
#'                                in \code{histories} is to be applied.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}. The default is \code{NA}.
#' @param ymodel                  Model statement for the outcome variable.
#' @param visitprocess            List of vectors. Each vector contains as its first entry
#'                                the covariate name of a visit process; its second entry
#'                                the name of a covariate whose modeling depends on the
#'                                visit process; and its third entry the maximum number
#'                                of consecutive visits that can be missed before an
#'                                individual is censored. The default is \code{NA}.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#'                                The default is \code{NA}.
#' @param restrictions            List of vectors. Each vector contains as its first entry a covariate for which
#'                                \emph{a priori} knowledge of its distribution is available; its second entry a condition
#'                                under which no knowledge of its distribution is available and that must be \code{TRUE}
#'                                for the distribution of that covariate given that condition to be estimated via a parametric
#'                                model or other fitting procedure; its third entry a function for estimating the distribution
#'                                of that covariate given the condition in the second entry is false such that \emph{a priori} knowledge
#'                                of the covariate distribution is available; and its fourth entry a value used by the function in the
#'                                third entry. The default is \code{NA}.
#' @param baselags                Logical scalar for specifying the convention used for lagi and lag_cumavgi terms in the model statements when pre-baseline times are not
#'                                included in \code{obs_data} and when the current time index, \eqn{t}, is such that \eqn{t < i}. If this argument is set to \code{FALSE}, the value
#'                                of all lagi and lag_cumavgi terms in this context are set to 0 (for non-categorical covariates) or the reference
#'                                level (for categorical covariates). If this argument is set to \code{TRUE}, the value of lagi and lag_cumavgi terms
#'                                are set to their values at time 0. The default is \code{FALSE}.
#' @param nsamples                Integer specifying the number of bootstrap samples to generate.
#'                                The default is 0.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of CPU cores to use in parallel
#'                                simulation. This argument is required when parallel is set to \code{TRUE}.
#'                                In many applications, users may wish to set this argument equal to \code{parallel::detectCores() - 1}.
#' @param sim_data_b              Logical scalar indicating whether to return the simulated data set. If bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0), this argument must be set to \code{FALSE}. The default is \code{FALSE}.
#' @param ci_method               Character string specifying the method for calculating the bootstrap 95\% confidence intervals, if applicable. The options are \code{"percentile"} and \code{"normal"}.
#' @param threads                 Integer specifying the number of threads to be used in \code{data.table}. See \code{\link[data.table]{setDTthreads}} for further details.
#' @param model_fits              Logical scalar indicating whether to return the fitted models. Note that if this argument is set to \code{TRUE}, the output of this function may use a lot of memory. The default is \code{FALSE}.
#' @param boot_diag               Logical scalar indicating whether to return the parametric g-formula estimates as well as the coefficients, standard errors, and variance-covariance matrices of the parameters of the fitted models in the bootstrap samples. The default is \code{FALSE}.
#' @param show_progress           Logical scalar indicating whether to print a progress bar for the number of bootstrap samples completed in the R console. This argument is only applicable when \code{parallel} is set to \code{FALSE} and bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0). The default is \code{TRUE}.
#' @param ipw_cutoff_quantile     Percentile by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param ipw_cutoff_value        Cutoff value by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param int_visit_type          Vector of logicals. The kth element is a logical specifying whether to carry forward the intervened value (rather than the natural value) of the treatment variables(s) when performing a carry forward restriction type for the kth intervention in \code{interventions}.
#'                                When the kth element is set to \code{FALSE}, the natural value of the treatment variable(s) in the kth intervention in \code{interventions} will be carried forward.
#'                                By default, this argument is set so that the intervened value of the treatment variable(s) is carried forward for all interventions.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#'
#' @return                        An object of class \code{gformula_continuous_eof}. The object is a list with the following components:
#' \item{result}{Results table containing the estimated mean outcome for all interventions (inculding natural course) at the last time point as well as the "cumulative percent intervened on" and the "average percent intervened on". If bootstrapping was used, the results table includes the bootstrap end-of-follow-up mean ratio, standard error, and 95\% confidence interval.}
#' \item{coeffs}{A list of the coefficients of the fitted models.}
#' \item{stderrs}{A list of the standard errors of the coefficients of the fitted models.}
#' \item{vcovs}{A list of the variance-covariance matrices of the parameters of the fitted models.}
#' \item{rmses}{A list of root mean square error (RMSE) values of the fitted models.}
#' \item{fits}{A list of the fitted models for the time-varying covariates and outcome. If \code{model_fits} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{sim_data}{A list of data tables of the simulated data. Each element in the list corresponds to one of the interventions. If the argument \code{sim_data_b} is set to \code{FALSE}, a value of \code{NA} is given.}
#' \item{IP_weights}{A numeric vector specifying the inverse probability weights. See "Details".}
#' \item{bootests}{A data.table containing the bootstrap replicates of the parametric g-formula estimates. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootcoeffs}{A list, where the kth element is a list containing the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootstderrs}{A list, where the kth element is a list containing the standard errors of the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootvcovs}{A list, where the kth element is a list containing the variance-covariance matrices of the parameters of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{...}{Some additional elements.}
#'
#' The results for the g-formula simulation under various interventions for the last time point are printed with the \code{\link{print.gformula_continuous_eof}} function. To generate graphs comparing the mean estimated and observed covariate values over time, use the \code{\link{print.gformula_continuous_eof}} function.
#'
#' @seealso \code{\link{gformula}}
#' @references Chiu YH, Wen L, McGrath S, Logan R, Dahabreh IJ, Hernán MA. Evaluating model specification when using the parametric g-formula in the presence of censoring. American Journal of Epidemiology. In press.
#' @references McGrath S, Lin V, Zhang Z, Petito LC, Logan RW, Hernán MA, and JG Young. gfoRmula: An R package for estimating the effects of sustained treatment strategies via the parametric g-formula. Patterns. 2020;1:100008.
#' @references Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period: application to the healthy worker survivor effect. Mathematical Modelling. 1986;7:1393–1512. [Errata (1987) in Computers and Mathematics with Applications 14, 917.-921. Addendum (1987) in Computers and Mathematics with Applications 14, 923-.945. Errata (1987) to addendum in Computers and Mathematics with Applications 18, 477.].
#' @examples
#'
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
#' gform_cont_eof
#' }
#'
#' @import data.table
#' @export
gformula_continuous_eof <- function(obs_data, id,
                                    time_name, covnames, covtypes,
                                    covparams, covfits_custom = NA,
                                    covpredict_custom = NA, histvars = NULL,
                                    histories = NA, basecovs = NA,
                                    outcome_name, ymodel,
                                    censor_name = NULL, censor_model = NA,
                                    intvars = NULL, interventions = NULL,
                                    int_times = NULL, int_descript = NULL, ref_int = 0,
                                    visitprocess = NA, restrictions = NA,
                                    yrestrictions = NA, baselags = FALSE, nsimul = NA,
                                    sim_data_b = FALSE,  seed, nsamples = 0,
                                    parallel = FALSE, ncores = NA,
                                    ci_method = 'percentile', threads,
                                    model_fits = FALSE, boot_diag = FALSE,
                                    show_progress = TRUE, ipw_cutoff_quantile = NULL,
                                    ipw_cutoff_value = NULL, int_visit_type = NULL, ...){

  lag_indicator <- lagavg_indicator <- cumavg_indicator <- c()
  lag_indicator <- update_lag_indicator(covparams$covmodels, lag_indicator)
  lagavg_indicator <- update_lagavg_indicator(covparams$covmodels, lagavg_indicator)
  cumavg_indicator <- update_cumavg_indicator(covparams$covmodels, cumavg_indicator)

  censor <- !(length(censor_model) == 1 && is.na(censor_model))

  if (!missing(ymodel)){
    lag_indicator <- update_lag_indicator(ymodel, lag_indicator)
    lagavg_indicator <- update_lagavg_indicator(ymodel, lagavg_indicator)
    cumavg_indicator <- update_cumavg_indicator(ymodel, cumavg_indicator)
  }
  if (censor){
    lag_indicator <- update_lag_indicator(censor_model, lag_indicator)
    lagavg_indicator <- update_lagavg_indicator(censor_model, lagavg_indicator)
    cumavg_indicator <- update_cumavg_indicator(censor_model, cumavg_indicator)
  }
  histvals <- list(lag_indicator = lag_indicator, lagavg_indicator = lagavg_indicator,
                   cumavg_indicator = cumavg_indicator)

  comprisk <- FALSE; comprisk2 <- FALSE

  if (!missing(threads)){
    setDTthreads(threads = threads)
  }
  else {
    threads <- getDTthreads()
  }
  outcome_type <- 'continuous_eof'
  compevent_model <- NA; compevent2_model <- NA
  compevent_name <- NULL; compevent2_name <- NULL
  compevent_restrictions <- NA
  hazardratio <- FALSE
  intcomp <- NA

  extra_args <- list(...)
  if ('time_points' %in% names(extra_args)){
    stop('Argument time_points cannot be supplied in this function. For end of follow up outcomes, the mean is calculated at the last time point in obs_data')
  }


  error_catch(id = id, nsimul = nsimul, intvars = intvars, interventions = interventions,
              int_times = int_times, int_descript = int_descript,
              covnames = covnames, covtypes = covtypes, basecovs = basecovs,
              histvars = histvars, histories = histories, compevent_model = compevent_model,
              hazardratio = hazardratio, intcomp = intcomp, time_points = NULL,
              outcome_type = outcome_type, time_name = time_name,
              obs_data = obs_data, parallel = parallel, ncores = ncores,
              nsamples = nsamples, sim_data_b = sim_data_b,
              outcome_name = outcome_name, compevent_name = compevent_name,
              comprisk = comprisk, censor = censor, censor_name = censor_name,
              covmodels = covparams$covmodels,
              histvals = histvals, ipw_cutoff_quantile = ipw_cutoff_quantile,
              ipw_cutoff_value = ipw_cutoff_value)

  min_time <- min(obs_data[[time_name]])
  below_zero_indicator <- min_time < 0

  obs_data <- copy(obs_data)




  max_visits <- NA
  if (!is.na(visitprocess[[1]][[1]])){
    for (vp in visitprocess){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(vp[1], paste("lag1_ts_", vp[1], "!=", vp[3], sep = ""),
                               ### rwl paste("visit_sum_", vp[3], "_", vp[1], "!=0", sep = ""),
                               simple_restriction, 1),
                             c(vp[2], paste(vp[1], "==1", sep = ""), carry_forward)))
      if (is.na(max_visits[1])){
        max_visits <- as.numeric(vp[3])
      } else {
        max_visits <- c(max_visits, as.numeric(vp[3]))
      }
      if (is.na(histories[1])){
        histories <- c(visit_sum)
      } else {
        histories <- c(visit_sum, histories)
        histvars <- append(list(c(vp[1])),histvars)

      }
    }
  }




  for (t in 0:max(obs_data[[time_name]])) {
    make_histories(pool = obs_data, histvars = histvars, histvals = histvals,
                   histories = histories, time_name = time_name, t = t, id = id ,
                   max_visits = max_visits, baselags = baselags,
                   below_zero_indicator = below_zero_indicator)
  }

  sample_size <- length(unique(obs_data[[id]]))
  time_points <- max(obs_data[[time_name]])+1

  for (i in seq_along(covnames)){
    if (covtypes[i] == 'absorbing'){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(covnames[i], paste("lag1_", covnames[i], "==0", sep = ""),
                               carry_forward, 1)))
      covtypes[i] <- 'binary'
    }
  }

  # Create 1-indexed numerical IDs for observed datasets
  ids <- as.data.table(sort(unique(obs_data[[id]])))
  ids[, 'newid' := seq_len(.N)]
  setkeyv(obs_data, id)
  obs_data <- obs_data[J(ids), allow.cartesian = TRUE]
  obs_data_geq_0 <- obs_data[obs_data[[time_name]] >= 0]

  # Set default number of simulated individuals to equal number of individuals in
  # observed dataset
  if (is.na(nsimul)){
    nsimul <- length(unique(obs_data$newid))
  }


  # Generate seeds for simulations and bootstrapping
  set.seed(seed)
  newseeds <- sample.int(2^30, size = nsamples + 1)
  subseed <- newseeds[1]
  bootseeds <- newseeds[2:(nsamples + 1)]

  # Determine ranges of observed covariates and outcome
  ranges <- lapply(seq_along(covnames), FUN = function(i){
    if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
        covtypes[i] == 'truncated normal') {
      range(obs_data_geq_0[[covnames[i]]])
    } else if (covtypes[i] == 'zero-inflated normal'){
      range(obs_data_geq_0[obs_data_geq_0[[covnames[i]]] > 0][[covnames[i]]])
    } else {
      NA
    }
  })

  # Fit models to covariates and outcome variable
  if (time_points > 1){
    fitcov <- pred_fun_cov(covparams = covparams, covnames = covnames, covtypes = covtypes,
                           covfits_custom = covfits_custom, restrictions = restrictions,
                           time_name = time_name, obs_data = obs_data_geq_0,
                           model_fits = model_fits)
    names(fitcov) <- covnames
  } else {
    fitcov <- NULL
  }
  fitY <- pred_fun_Y(ymodel, yrestrictions, outcome_type, outcome_name, time_name, obs_data_geq_0,
                     model_fits = model_fits)
  if (censor){
    fitC <- pred_fun_D(censor_model, NA, obs_data_geq_0, model_fits = model_fits)
  } else {
    fitC <- NA
  }

  obs_data_noresample <- copy(obs_data)
  len <- length(unique(obs_data$newid))
  # If the number of user desired simulations differs from the number of individuals in
  # the observed dataset, sample the desired number of observed IDs with replacement
  if (nsimul < len){
    ids <- as.data.table(sort(sample(unique(obs_data$newid), nsimul, replace = TRUE)))
    colnames(ids) <- "newid"
    ids[, 'sid' := seq_len(.N)]
    obs_data <- merge(ids, obs_data, all.x = TRUE, by = "newid")
    obs_data$newid <- obs_data$sid
    obs_data$sid <- NULL
  } else if (nsimul > len){
    ids <- as.data.table(sample(unique(obs_data$newid), nsimul, replace = TRUE))
    ids$newid <- 1:nsimul
    colnames(ids) <- c("newid", "sid")
    setkeyv(obs_data, "newid")
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE]  # create the new data set names "sample"
    obs_data$newid <- obs_data$sid
    obs_data$sid <- NULL
  }

  # Add natural course to list of interventions
  if (!is.null(interventions)){
    comb_interventions <- c(list(list(c(natural))), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    comb_interventions <- list(list(c(natural)))
    comb_intvars <- list('none')
  }

  if (is.null(int_times)){
    comb_int_times <- list()
    for (i in seq_along(comb_interventions)){
      comb_int_times[[i]] <- lapply(seq_along(comb_interventions[[i]]),
                                    FUN = function(i) {0:(time_points - 1)})
    }
  } else {
    comb_int_times <- c(list(list(0:(time_points - 1))), int_times)
  }

  if (is.null(int_visit_type)){
    int_visit_type <- rep(TRUE, length(comb_interventions))
  } else {
    int_visit_type <- c(TRUE, int_visit_type)
  }

  if (parallel){
    cl <- prep_cluster(ncores = ncores, threads = threads,  covtypes = covtypes)
    pools <- parallel::parLapply(cl, seq_along(comb_interventions), simulate,
                                 fitcov = fitcov, fitY = fitY, fitD = NA,
                                 yrestrictions = yrestrictions,
                                 compevent_restrictions = compevent_restrictions,
                                 restrictions = restrictions,
                                 outcome_name = outcome_name, compevent_name = compevent_name,
                                 time_name = time_name,
                                 intvars = comb_intvars, interventions = comb_interventions,
                                 int_times = comb_int_times, histvars = histvars, histvals = histvals, histories = histories,
                                 covparams = covparams, covnames = covnames, covtypes = covtypes,
                                 covpredict_custom = covpredict_custom, basecovs = basecovs,
                                 comprisk = comprisk, ranges = ranges,
                                 outcome_type = outcome_type,
                                 subseed = subseed, time_points = time_points,
                                 obs_data = obs_data, parallel = parallel,
                                 baselags = baselags, below_zero_indicator = below_zero_indicator,
                                 min_time = min_time, show_progress = FALSE, int_visit_type = int_visit_type, ...)
    parallel::stopCluster(cl)
  } else {
    pools <- lapply(seq_along(comb_interventions), FUN = function(i){
      simulate(fitcov = fitcov, fitY = fitY, fitD = NA,
               yrestrictions = yrestrictions,
               compevent_restrictions = compevent_restrictions,
               restrictions = restrictions,
               outcome_name = outcome_name, compevent_name = compevent_name,
               time_name = time_name,
               intvars = comb_intvars[[i]], interventions = comb_interventions[[i]],
               int_times = comb_int_times[[i]], histvars = histvars, histvals = histvals, histories = histories,
               covparams = covparams, covnames = covnames, covtypes = covtypes,
               covpredict_custom = covpredict_custom, basecovs = basecovs, comprisk = comprisk,
               ranges = ranges,
               outcome_type = outcome_type,
               subseed = subseed, time_points = time_points,
               obs_data = obs_data, parallel = parallel,
               baselags = baselags, below_zero_indicator = below_zero_indicator,
               min_time = min_time, show_progress = FALSE, int_visit_type = int_visit_type[i], ...)
    })
  }

  nat_pool <- pools[[1]] # Natural course data
  pools <- pools[-1] # List of intervention datasets

  # Initialize results matrices
  result_ratio <- result_diff <- int_result <- rep(NA, length(pools) + 1)

  # Calculate mean outcome over all subjects at each time for natural course
  nat_result <- mean(nat_pool$Ey, na.rm = TRUE)

  if (ref_int == 0){
    # Set reference intervention to the natural course
    ref_mean <- nat_result
  } else {
    # Set reference intervention as specified
    # Calculate mean outcome over all subjects at each time for this intervention
    ref_mean <- mean(pools[[ref_int]]$Ey, na.rm = TRUE)
  }

  # Compile results
  int_result[1] <- nat_result
  # Calculate mean risk over all subjects at each time for all interventions other than
  # the natural course
  int_result[-1] <- sapply(pools, FUN = function(pool){mean(pool$Ey, na.rm = TRUE)})
  result_ratio <- int_result / ref_mean
  result_diff <- int_result - ref_mean

  # Calculate percent intervened
  percent_intervened_res <- get_percent_intervened(pools = pools)

  if (nsamples > 0){
    if (parallel){
      cl <- prep_cluster(ncores = ncores, threads = threads, covtypes = covtypes,
                         bootstrap_option = TRUE)
      final_bs <- parallel::parLapply(cl, 1:nsamples, bootstrap_helper_with_trycatch, time_points = time_points,
                                      obs_data = obs_data_noresample, bootseeds = bootseeds,
                                      intvars = comb_intvars, interventions = comb_interventions, int_times = comb_int_times, ref_int = ref_int,
                                      covparams = covparams, covnames = covnames, covtypes = covtypes,
                                      covfits_custom = covfits_custom, covpredict_custom = covpredict_custom,
                                      basecovs = basecovs, ymodel = ymodel,
                                      histvars = histvars, histvals = histvals, histories = histories,
                                      comprisk = comprisk, compevent_model = compevent_model,
                                      yrestrictions = yrestrictions,
                                      compevent_restrictions = compevent_restrictions,
                                      restrictions = restrictions, outcome_type = outcome_type,
                                      ranges = ranges,
                                      time_name = time_name, outcome_name = outcome_name,
                                      compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                                      max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp,
                                      boot_diag = boot_diag, nsimul = nsimul, baselags = baselags,
                                      below_zero_indicator = below_zero_indicator, min_time = min_time,
                                      show_progress = FALSE, int_visit_type = int_visit_type, ...)
      parallel::stopCluster(cl)

    } else {
      if (show_progress){
        pb <- progress::progress_bar$new(total = nsamples * length(comb_interventions),
                                         clear = FALSE,
                                         format = 'Bootstrap progress [:bar] :percent, Elapsed time :elapsed, Est. time remaining :eta')
      }
      final_bs <- lapply(1:nsamples, FUN = bootstrap_helper_with_trycatch, time_points = time_points,
                         obs_data = obs_data_noresample, bootseeds = bootseeds,
                         intvars = comb_intvars, interventions = comb_interventions, int_times = comb_int_times, ref_int = ref_int,
                         covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, covpredict_custom = covpredict_custom,
                         basecovs = basecovs, ymodel = ymodel,
                         histvars = histvars, histvals = histvals, histories = histories,
                         comprisk = comprisk, compevent_model = compevent_model,
                         yrestrictions = yrestrictions,
                         compevent_restrictions = compevent_restrictions,
                         restrictions = restrictions,
                         outcome_type = outcome_type,
                         ranges = ranges,
                         time_name = time_name, outcome_name = outcome_name,
                         compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                         max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp,
                         boot_diag = boot_diag, nsimul = nsimul, baselags = baselags,
                         below_zero_indicator = below_zero_indicator, min_time = min_time,
                         show_progress = show_progress, pb = pb, int_visit_type = int_visit_type, ...)
    }
    comb_result <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$Result))
    }))
    comb_MR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultRatio))
    }))
    comb_MD <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultDiff))
    }))

    comb_result$t0 <- comb_MR$t0 <- comb_MD$t0 <- time_points

    se_result <- comb_result[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]
    se_MR <- comb_MR[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]
    se_MD <- comb_MD[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]

    if (ci_method == 'normal'){
      ci_lb_result <- t(int_result) - stats::qnorm(0.975)*se_result[,-c('t0')]
      ci_lb_MR <- t(int_result) - stats::qnorm(0.975)*se_MR[,-c('t0')]
      ci_lb_MD <- t(int_result) - stats::qnorm(0.975)*se_MD[,-c('t0')]
      ci_ub_result <- t(int_result) + stats::qnorm(0.975)*se_result[,-c('t0')]
      ci_ub_MR <- t(int_result) + stats::qnorm(0.975)*se_MR[,-c('t0')]
      ci_ub_MD <- t(int_result) + stats::qnorm(0.975)*se_MD[,-c('t0')]
    }
    if (ci_method == 'percentile') {
      ci_lb_result <- comb_result[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_lb_MR <- comb_MR[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_lb_MD <- comb_MD[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_ub_result <- comb_result[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
      ci_ub_MR <- comb_MR[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
      ci_ub_MD <- comb_MD[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
    }
  }
  if (nsamples > 0 & boot_diag){
    bootests <- comb_result
    if (!is.null(int_descript)){
      colnames(bootests)[1:(1 + length(interventions))] <- c('Natural course', int_descript)
    } else {
      colnames(bootests)[1:length(comb_interventions)] <- c('Natural course', paste('Intervention', 1:length(interventions)))
    }
    bootests[, 'Bootstrap replicate'] <- 1:nsamples
    bootcoeffs <- lapply(final_bs, "[[", 'bootcoeffs')
    bootstderrs <- lapply(final_bs, "[[", 'bootstderrs')
    bootvcovs <- lapply(final_bs, "[[", 'bootvcovs')
  } else {
    bootests <- NULL
    bootcoeffs <- NULL
    bootstderrs <- NULL
    bootvcovs <- NULL
  }

  plot_info <- get_plot_info(outcome_name = outcome_name,
                             compevent_name = compevent_name,
                             compevent2_name = compevent2_name,
                             censor_name = censor_name,
                             time_name = time_name,
                             id = id,
                             time_points = time_points,
                             covnames = covnames,
                             covtypes = covtypes,
                             nat_pool = nat_pool,
                             nat_result = nat_result,
                             comprisk = comprisk,
                             comprisk2 = comprisk2,
                             censor = censor,
                             fitC = fitC,
                             outcome_type = outcome_type,
                             obs_data = obs_data_noresample,
                             ipw_cutoff_quantile = ipw_cutoff_quantile,
                             ipw_cutoff_value = ipw_cutoff_value)
  obs_results <- plot_info$obs_results

  # Generate results table
  if (!is.null(interventions)){
    resultdf <- lapply(seq_along(int_result), function(k){
      if (nsamples > 0){
        data.table(t = time_points - 1, Intervention = k - 1,
                   EOFMean = int_result[k],
                   EOFMean_SE = se_result[[paste0('V',k)]],
                   EOFMean_CI_LL95 = ci_lb_result[[paste0('V',k)]],
                   EOFMean_CI_UL95 = ci_ub_result[[paste0('V',k)]],
                   EOFMeanRatio = result_ratio[k],
                   MR_SE = se_MR[[paste0('V',k)]],
                   MR_CI_LL95 = ci_lb_MR[[paste0('V',k)]],
                   MR_CI_UL95 = ci_ub_MR[[paste0('V',k)]],
                   EOFMeanDiff = result_diff[k],
                   MD_SE = se_MD[[paste0('V',k)]],
                   MD_CI_LL95 = ci_lb_MD[[paste0('V',k)]],
                   MD_CI_UL95 = ci_ub_MD[[paste0('V',k)]])
      } else {
        data.table(t = time_points - 1, Intervention = k - 1, OutcomeEOFMean = int_result[k],
                   EOFMeanRatio = result_ratio[k], EOFMeanDiff = result_diff[k])
      }
    })
    resultdf <- rbindlist(resultdf)
  } else {
    if (nsamples > 0){
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             EOFMean = int_result,
                             EOFMean_SE = se_MD[,-c('t0')], EOFMean_CI_LL95 = ci_lb_result[,-c('t0')],
                             EOFMean_CI_UL95 = ci_ub_result[,-c('t0')], EOFMeanRatio = result_ratio,
                             MR_SE = se_MR[,-c('t0')],
                             MR_CI_LL95 = ci_lb_MR[,-c('t0')], MR_CI_UL95 = ci_ub_MR[,-c('t0')],
                             EOFMeanDiff = result_diff,
                             EOFMD_se = se_MD[,-c('t0')],
                             MD_CI_LL95 = ci_lb_MD[,-c('t0')], MD_CI_UL95 = ci_ub_MD[,-c('t0')])
    } else {
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             OutcomeEOFMean = int_result, EOFMeanRatio = result_ratio,
                             EOFMeanDiff = result_diff)
    }
  }

  resultdf$obs_risk <- c(utils::tail(obs_results[[2]], 1), rep(NA, dim(resultdf)[1] - 1))
  obs_mean_name <- ifelse(censor, 'IP weighted mean', 'NP mean')
  if (nsamples > 0){
    colnames(resultdf) <- c("k", "Interv.", "g-form mean",
                            "Mean SE",
                            "Mean lower 95% CI", "Mean upper 95% CI",
                            "Mean ratio",
                            "MR SE", "MR lower 95% CI", "MR upper 95% CI",
                            "Mean difference",
                            "MD SE", "MD lower 95% CI", "MD upper 95% CI",
                            obs_mean_name)
    setcolorder(resultdf, c("k", "Interv.", obs_mean_name, "g-form mean",
                            "Mean SE",
                            "Mean lower 95% CI", "Mean upper 95% CI",
                            "Mean ratio",
                            "MR SE", "MR lower 95% CI", "MR upper 95% CI",
                            "Mean difference",
                            "MD SE", "MD lower 95% CI", "MD upper 95% CI"))
  } else {
    colnames(resultdf) <- c("k", "Interv.", "g-form mean", "Mean ratio",
                            "Mean difference", obs_mean_name)
    setcolorder(resultdf, c("k", "Interv.", obs_mean_name, "g-form mean",
                            "Mean ratio", "Mean difference"))
  }

  resultdf[k == max(k), '% Intervened On'] <- percent_intervened_res$percent_intervened
  resultdf[k == max(k), 'Aver % Intervened On'] <- percent_intervened_res$average_percent_intervened

  if (time_points > 1){
    fits <- fitcov
    fits[[length(fits) + 1]] <- fitY
    names(fits)[length(fits)] <- outcome_name
  } else {
    fits <- list(fitY)
  }
  if (!is.na(fitC)[[1]]){
    fits[[length(fits) + 1]] <- fitC
    names(fits)[length(fits)] <- censor_name
  }

  # Add list of coefficients for covariates, outcome variable, and competing event
  # variable (if any) to results output
  coeffs <- get_coeffs(fits = fits, fitD = NA, time_points = time_points,
                       outcome_name = outcome_name, compevent_name = compevent_name,
                       covnames = covnames)
  stderrs <- get_stderrs(fits = fits, fitD = NA, time_points = time_points,
                         outcome_name = outcome_name, compevent_name = compevent_name,
                         covnames = covnames)
  vcovs <- get_vcovs(fits = fits, fitD = NA, time_points = time_points,
                     outcome_name = outcome_name, compevent_name = compevent_name,
                     covnames = covnames)


  rmses <- lapply(seq_along(fits), FUN = rmse_calculate, fits = fits, covnames = covnames,
                  covtypes = covtypes)

  if (time_points == 1){
    rmses <- stats::setNames(rmses, outcome_name)
  } else {
    rmses <- stats::setNames(rmses, c(covnames, outcome_name))
  }


  # Create header
  header <- get_header(int_descript, sample_size, nsimul, nsamples, ref_int)

  if (sim_data_b){
    sim_data <- c(list('Natural course' = nat_pool), pools)
    if (!is.null(int_descript)){
      names(sim_data)[2:length(sim_data)] <- int_descript
    }
  } else {
    sim_data <- NA
  }
  if (!model_fits){
    fits <- NULL
  }

  res <- list(
    result = resultdf,
    coeffs = coeffs,
    stderrs = stderrs,
    vcovs = vcovs,
    rmses = rmses,
    fits = fits,
    sim_data = sim_data,
    IP_weights = obs_results$w,
    bootests = bootests,
    bootcoeffs = bootcoeffs,
    bootstderrs = bootstderrs,
    bootvcovs = bootvcovs,
    time_name = time_name,
    time_points = time_points,
    covnames = covnames,
    covtypes = covtypes,
    dt_cov_plot = plot_info$dt_cov_plot,
    dt_out_plot = plot_info$dt_out_plot,
    header = header
  )
  class(res) <- c("gformula_continuous_eof", "gformula")
  return (res)
}

#' Estimation of Binary End-of-Follow-Up Outcome Under the Parametric G-Formula
#'
#' Based on an observed data set, this internal function estimates the outcome probability at
#' end-of-follow-up under multiple user-specified interventions using the parametric g-formula. See McGrath et al. (2020) for
#' further details concerning the application and implementation of the parametric g-formula.
#'
#' To assess model misspecification in the parametric g-formula, users can obtain inverse probability (IP) weighted estimates of the natural course means of the time-varying covariates from the observed data.
#' See Chiu et al. (In press) for details.
#' In addition to the general requirements described in McGrath et al. (2020), the requirements for the input data set and the call to the gformula function for such analyses are described below.
#'
#' Users need to include a column in \code{obs_data} with a time-varying censoring variable.
#' Users need to indicate the name of the censoring variable and a model statement for the censoring variable with parameters \code{censor_name} and \code{censor_model}, respectively.
#' Finally, users can specify how to truncate IP weights with the \code{ipw_cutoff_quantile} or \code{ipw_cutoff_value} parameters.
#'
#' In addition to the package output described in McGrath et al. (2020), the output will display estimates of the "cumulative percent intervened on" and the "average percent intervened on". When using a custom intervention function, users need to specify whether each individual at that time point is eligible to contribute person-time to the percent intervened on calculations. Specifically, this must be specified in the \code{eligible_pt} column of \code{newdf}. By default, \code{eligible_pt} is set to \code{TRUE} for each individual at each time point in custom interventions.
#'
#' @param id                      Character string specifying the name of the ID variable in \code{obs_data}.
#' @param obs_data                Data table containing the observed data.
#' @param seed                    Starting seed for simulations and bootstrapping.
#' @param nsimul                  Number of subjects for whom to simulate data. By default, this argument is set
#'                                equal to the number of subjects in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param censor_name             Character string specifying the name of the censoring variable in \code{obs_data}. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param censor_model            Model statement for the censoring variable. Only applicable when using inverse probability weights to estimate the natural course means / risk from the observed data. See "Details".
#' @param intvars                 List, whose elements are vectors of character strings. The kth vector in \code{intvars} specifies the name(s) of the variable(s) to be intervened
#'                                on in each round of the simulation under the kth intervention in \code{interventions}.
#' @param interventions           List, whose elements are lists of vectors. Each list in \code{interventions} specifies a unique intervention on the relevant variable(s) in \code{intvars}. Each vector contains a function
#'                                implementing a particular intervention on a single variable, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#' @param int_times               List, whose elements are lists of vectors. The kth list in \code{int_times} corresponds to the kth intervention in \code{interventions}. Each vector specifies the time points in which the relevant intervention is applied on the corresponding variable in \code{intvars}.
#'                                When an intervention is not applied, the simulated natural course value is used. By default, this argument is set so that all interventions are applied in all time points.
#' @param int_descript            Vector of character strings, each describing an intervention. It must
#'                                be in same order as the entries in \code{interventions}.
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the end-of-follow-up mean ratio and mean difference. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}. The default is 0.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, \code{"absorbing"}, \code{"categorical time"}, and \code{"custom"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}. The default is \code{NA}.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}. The default is \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}. These covariates are not simulated using a model but rather carry their value over all time points from the first time point of \code{obs_data}. These covariates should not be included in \code{covnames}. The default is \code{NA}.
#' @param histvars                List of vectors. The kth vector specifies the names of the variables for which the kth history function
#'                                in \code{histories} is to be applied.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}. The default is \code{NA}.
#' @param ymodel                  Model statement for the outcome variable.
#' @param visitprocess            List of vectors. Each vector contains as its first entry
#'                                the covariate name of a visit process; its second entry
#'                                the name of a covariate whose modeling depends on the
#'                                visit process; and its third entry the maximum number
#'                                of consecutive visits that can be missed before an
#'                                individual is censored. The default is \code{NA}.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#'                                The default is \code{NA}.
#' @param restrictions            List of vectors. Each vector contains as its first entry a covariate for which
#'                                \emph{a priori} knowledge of its distribution is available; its second entry a condition
#'                                under which no knowledge of its distribution is available and that must be \code{TRUE}
#'                                for the distribution of that covariate given that condition to be estimated via a parametric
#'                                model or other fitting procedure; its third entry a function for estimating the distribution
#'                                of that covariate given the condition in the second entry is false such that \emph{a priori} knowledge
#'                                of the covariate distribution is available; and its fourth entry a value used by the function in the
#'                                third entry. The default is \code{NA}.
#' @param baselags                Logical scalar for specifying the convention used for lagi and lag_cumavgi terms in the model statements when pre-baseline times are not
#'                                included in \code{obs_data} and when the current time index, \eqn{t}, is such that \eqn{t < i}. If this argument is set to \code{FALSE}, the value
#'                                of all lagi and lag_cumavgi terms in this context are set to 0 (for non-categorical covariates) or the reference
#'                                level (for categorical covariates). If this argument is set to \code{TRUE}, the value of lagi and lag_cumavgi terms
#'                                are set to their values at time 0. The default is \code{FALSE}.
#' @param nsamples                Integer specifying the number of bootstrap samples to generate.
#'                                The default is 0.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of CPU cores to use in parallel
#'                                simulation. This argument is required when parallel is set to \code{TRUE}.
#'                                In many applications, users may wish to set this argument equal to \code{parallel::detectCores() - 1}.
#' @param sim_data_b              Logical scalar indicating whether to return the simulated data set. If bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0), this argument must be set to \code{FALSE}. The default is \code{FALSE}.
#' @param ci_method               Character string specifying the method for calculating the bootstrap 95\% confidence intervals, if applicable. The options are \code{"percentile"} and \code{"normal"}.
#' @param threads                 Integer specifying the number of threads to be used in \code{data.table}. See \code{\link[data.table]{setDTthreads}} for further details.
#' @param model_fits              Logical scalar indicating whether to return the fitted models. Note that if this argument is set to \code{TRUE}, the output of this function may use a lot of memory. The default is \code{FALSE}.
#' @param boot_diag               Logical scalar indicating whether to return the parametric g-formula estimates as well as the coefficients, standard errors, and variance-covariance matrices of the parameters of the fitted models in the bootstrap samples. The default is \code{FALSE}.
#' @param show_progress           Logical scalar indicating whether to print a progress bar for the number of bootstrap samples completed in the R console. This argument is only applicable when \code{parallel} is set to \code{FALSE} and bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0). The default is \code{TRUE}.
#' @param ipw_cutoff_quantile     Percentile by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param ipw_cutoff_value        Cutoff value by which to truncate inverse probability weights. The default is \code{NULL} (i.e., no truncation). See "Details".
#' @param int_visit_type          Vector of logicals. The kth element is a logical specifying whether to carry forward the intervened value (rather than the natural value) of the treatment variables(s) when performing a carry forward restriction type for the kth intervention in \code{interventions}.
#'                                When the kth element is set to \code{FALSE}, the natural value of the treatment variable(s) in the kth intervention in \code{interventions} will be carried forward.
#'                                By default, this argument is set so that the intervened value of the treatment variable(s) is carried forward for all interventions.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#'
#' @return An object of class \code{gformula_binary_eof}. The object is a list with the following components:
#' \item{result}{Results table containing the estimated outcome probability for all interventions (inculding natural course) at the last time point as well as the "cumulative percent intervened on" and the "average percent intervened on". If bootstrapping was used, the results table includes the bootstrap end-of-follow-up mean ratio, standard error, and 95\% confidence interval.}
#' \item{coeffs}{A list of the coefficients of the fitted models.}
#' \item{stderrs}{A list of the standard errors of the coefficients of the fitted models.}
#' \item{vcovs}{A list of the variance-covariance matrices of the parameters of the fitted models.}
#' \item{rmses}{A list of root mean square error (RMSE) values of the fitted models.}
#' \item{fits}{A list of the fitted models for the time-varying covariates and outcome. If \code{model_fits} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{sim_data}{A list of data tables of the simulated data. Each element in the list corresponds to one of the interventions. If the argument \code{sim_data_b} is set to \code{FALSE}, a value of \code{NA} is given.}
#' \item{IP_weights}{A numeric vector specifying the inverse probability weights. See "Details".}
#' \item{bootests}{A data.table containing the bootstrap replicates of the parametric g-formula estimates. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootcoeffs}{A list, where the kth element is a list containing the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootstderrs}{A list, where the kth element is a list containing the standard errors of the coefficients of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{bootvcovs}{A list, where the kth element is a list containing the variance-covariance matrices of the parameters of the fitted models corresponding to the kth bootstrap sample. If \code{boot_diag} is set to \code{FALSE}, a value of \code{NULL} is given.}
#' \item{...}{Some additional elements.}
#'
#' The results for the g-formula simulation under various interventions for the last time point are printed with the \code{\link{print.gformula_binary_eof}} function. To generate graphs comparing the mean estimated and observed covariate values over time, use the \code{\link{plot.gformula_binary_eof}} function.
#'
#' @seealso \code{\link{gformula}}
#' @references Chiu YH, Wen L, McGrath S, Logan R, Dahabreh IJ, Hernán MA. Evaluating model specification when using the parametric g-formula in the presence of censoring. American Journal of Epidemiology. In press.
#' @references McGrath S, Lin V, Zhang Z, Petito LC, Logan RW, Hernán MA, and JG Young. gfoRmula: An R package for estimating the effects of sustained treatment strategies via the parametric g-formula. Patterns. 2020;1:100008.
#' @references Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period: application to the healthy worker survivor effect. Mathematical Modelling. 1986;7:1393–1512. [Errata (1987) in Computers and Mathematics with Applications 14, 917.-921. Addendum (1987) in Computers and Mathematics with Applications 14, 923-.945. Errata (1987) to addendum in Computers and Mathematics with Applications 18, 477.].
#' @examples
#'
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
#' gform_bin_eof
#' }
#'
#' @import data.table
#' @export
gformula_binary_eof <- function(obs_data, id,
                                time_name, covnames, covtypes, covparams,
                                covfits_custom = NA, covpredict_custom = NA,
                                histvars = NULL, histories = NA, basecovs = NA,
                                censor_name = NULL, censor_model = NA,
                                outcome_name, ymodel, intvars = NULL,
                                interventions = NULL, int_times = NULL, int_descript = NULL,
                                ref_int = 0, visitprocess = NA,
                                restrictions = NA, yrestrictions = NA, baselags = FALSE,
                                nsimul = NA, sim_data_b = FALSE, seed,
                                nsamples = 0, parallel = FALSE, ncores = NA,
                                ci_method = 'percentile', threads,
                                model_fits = FALSE, boot_diag = FALSE,
                                show_progress = TRUE, ipw_cutoff_quantile = NULL,
                                ipw_cutoff_value = NULL, int_visit_type = NULL, ...){

  lag_indicator <- lagavg_indicator <- cumavg_indicator <- c()
  lag_indicator <- update_lag_indicator(covparams$covmodels, lag_indicator)
  lagavg_indicator <- update_lagavg_indicator(covparams$covmodels, lagavg_indicator)
  cumavg_indicator <- update_cumavg_indicator(covparams$covmodels, cumavg_indicator)

  censor <- !(length(censor_model) == 1 && is.na(censor_model))

  if (!missing(ymodel)){
    lag_indicator <- update_lag_indicator(ymodel, lag_indicator)
    lagavg_indicator <- update_lagavg_indicator(ymodel, lagavg_indicator)
    cumavg_indicator <- update_cumavg_indicator(ymodel, cumavg_indicator)
  }
  if (censor){
    lag_indicator <- update_lag_indicator(censor_model, lag_indicator)
    lagavg_indicator <- update_lagavg_indicator(censor_model, lagavg_indicator)
    cumavg_indicator <- update_cumavg_indicator(censor_model, cumavg_indicator)
  }
  histvals <- list(lag_indicator = lag_indicator, lagavg_indicator = lagavg_indicator,
                   cumavg_indicator = cumavg_indicator)
  comprisk <- FALSE; comprisk2 <- FALSE

  if (!missing(threads)){
    setDTthreads(threads = threads)
  }
  else {
    threads <- getDTthreads()
  }
  outcome_type <- 'binary_eof'
  compevent_model <- NA; compevent2_model <- NA
  compevent_name <- NULL; compevent2_name <- NULL
  compevent_restrictions <- NA
  hazardratio <- FALSE
  intcomp <- NA

  extra_args <- list(...)
  if ('time_points' %in% names(extra_args)){
    stop('Argument time_points cannot be supplied in this function. For end of follow up outcomes, the mean is calculated at the last time point in obs_data')
  }


  error_catch(id = id, nsimul = nsimul, intvars = intvars, interventions = interventions,
              int_times = int_times, int_descript = int_descript,
              covnames = covnames, covtypes = covtypes, basecovs = basecovs,
              histvars = histvars, histories = histories, compevent_model = compevent_model,
              hazardratio = hazardratio, intcomp = intcomp, time_points = NULL,
              outcome_type = outcome_type, time_name = time_name,
              obs_data = obs_data, parallel = parallel, ncores = ncores,
              nsamples = nsamples, sim_data_b = sim_data_b,
              outcome_name = outcome_name, compevent_name = compevent_name,
              comprisk = comprisk, censor = censor, censor_name = censor_name,
              covmodels = covparams$covmodels,
              histvals = histvals, ipw_cutoff_quantile = ipw_cutoff_quantile,
              ipw_cutoff_value = ipw_cutoff_value)

  min_time <- min(obs_data[[time_name]])
  below_zero_indicator <- min_time < 0

  obs_data <- copy(obs_data)



  max_visits <- NA
  if (!is.na(visitprocess[[1]][[1]])){
    for (vp in visitprocess){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(vp[1], paste("lag1_ts_", vp[1], "!=",vp[3], sep = ""),
                               ### rwl paste("visit_sum_", vp[3], "_", vp[1], "!=0", sep = ""),
                               simple_restriction, 1),
                             c(vp[2], paste(vp[1], "==1", sep = ""), carry_forward)))
      if (is.na(max_visits[1])){
        max_visits <- as.numeric(vp[3])
      } else {
        max_visits <- c(max_visits, as.numeric(vp[3]))
      }
      if (is.na(histories[1])){
        histories <- c(visit_sum)
      } else {
        histories <- c(visit_sum, histories)
        histvars <- append(list(c(vp[1])),histvars)

      }
    }
  }


  for (t in 0:max(obs_data[[time_name]])) {
    make_histories(pool = obs_data, histvars = histvars, histvals = histvals,
                   histories = histories, time_name = time_name, t = t, id = id,
                   baselags = baselags, below_zero_indicator = below_zero_indicator)
  }

  sample_size <- length(unique(obs_data[[id]]))
  time_points <- max(obs_data[[time_name]])+1

  for (i in seq_along(covnames)){
    if (covtypes[i] == 'absorbing'){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(covnames[i], paste("lag1_", covnames[i], "==0", sep = ""),
                               carry_forward, 1)))
      covtypes[i] <- 'binary'
    }
  }

  # Create 1-indexed numerical IDs for observed datasets
  ids <- as.data.table(sort(unique(obs_data[[id]])))
  ids[, "newid" := seq_len(.N)]
  setkeyv(obs_data, id)
  obs_data <- obs_data[J(ids), allow.cartesian = TRUE]
  obs_data_geq_0 <- obs_data[obs_data[[time_name]] >= 0]

  # Set default number of simulated individuals to equal number of individuals in
  # observed dataset
  if (is.na(nsimul)){
    nsimul <- length(unique(obs_data$newid))
  }

  # Generate seeds for simulations and bootstrapping
  set.seed(seed)
  newseeds <- sample.int(2^30, size = nsamples + 1)
  subseed <- newseeds[1]
  bootseeds <- newseeds[2:(nsamples + 1)]

  # Determine ranges of observed covariates and outcome
  ranges <- lapply(seq_along(covnames), FUN = function(i){
    if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
        covtypes[i] == 'truncated normal') {
      range(obs_data_geq_0[[covnames[i]]])
    } else if (covtypes[i] == 'zero-inflated normal'){
      range(obs_data_geq_0[obs_data_geq_0[[covnames[i]]] > 0][[covnames[i]]])
    } else {
      NA
    }
  })

  # Fit models to covariates and outcome variable
  if (time_points > 1){
    fitcov <- pred_fun_cov(covparams = covparams, covnames = covnames, covtypes = covtypes,
                           covfits_custom = covfits_custom, restrictions = restrictions,
                           time_name = time_name, obs_data = obs_data_geq_0,
                           model_fits = model_fits)
    names(fitcov) <- covnames
  } else {
    fitcov <- NULL
  }
  fitY <- pred_fun_Y(ymodel, yrestrictions, outcome_type, outcome_name, time_name, obs_data_geq_0,
                     model_fits = model_fits)
  if (censor){
    fitC <- pred_fun_D(censor_model, NA, obs_data_geq_0, model_fits = model_fits)
  } else {
    fitC <- NA
  }

  obs_data_noresample <- copy(obs_data)
  len <- length(unique(obs_data$newid))
  # If the number of user desired simulations differs from the number of individuals in
  # the observed dataset, sample the desired number of observed IDs with replacement
  if (nsimul < len){
    ids <- as.data.table(sort(sample(unique(obs_data$newid), nsimul, replace = TRUE)))
    colnames(ids) <- "newid"
    ids[, 'sid' := seq_len(.N)]
    obs_data <- merge(ids, obs_data, all.x = TRUE, by = "newid")
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  } else if (nsimul > len){
    ids <- as.data.table(sample(unique(obs_data$newid), nsimul, replace = TRUE))
    ids[, 'newid' := 1:nsimul]
    colnames(ids) <- c("newid", "sid")
    setkeyv(obs_data, "newid")
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE]  # create the new data set names "sample"
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  }

  # Add natural course to list of interventions
  if (!is.null(interventions)){
    comb_interventions <- c(list(list(c(natural))), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    comb_interventions <- list(list(c(natural)))
    comb_intvars <- list('none')
  }

  if (is.null(int_times)){
    comb_int_times <- list()
    for (i in seq_along(comb_interventions)){
      comb_int_times[[i]] <- lapply(seq_along(comb_interventions[[i]]),
                                    FUN = function(i) {0:(time_points - 1)})
    }
  } else {
    comb_int_times <- c(list(list(0:(time_points - 1))), int_times)
  }

  if (is.null(int_visit_type)){
    int_visit_type <- rep(TRUE, length(comb_interventions))
  } else {
    int_visit_type <- c(TRUE, int_visit_type)
  }

  if (parallel){
    cl <- prep_cluster(ncores = ncores, threads = threads , covtypes = covtypes)
    pools <- parallel::parLapply(cl, seq_along(comb_interventions), simulate,
                                 fitcov = fitcov, fitY = fitY, fitD = NA,
                                 yrestrictions = yrestrictions,
                                 compevent_restrictions = compevent_restrictions,
                                 restrictions = restrictions,
                                 outcome_name = outcome_name, compevent_name = compevent_name,
                                 time_name = time_name,
                                 intvars = comb_intvars, interventions = comb_interventions,
                                 int_times = comb_int_times, histvars = histvars, histvals = histvals, histories = histories,
                                 covparams = covparams, covnames = covnames, covtypes = covtypes,
                                 covpredict_custom = covpredict_custom, basecovs = basecovs,
                                 comprisk = comprisk, ranges = ranges,
                                 outcome_type = outcome_type,
                                 subseed = subseed, time_points = time_points,
                                 obs_data = obs_data, parallel = parallel,
                                 baselags = baselags, below_zero_indicator = below_zero_indicator,
                                 min_time = min_time, show_progress = FALSE, int_visit_type = int_visit_type, ...)
    parallel::stopCluster(cl)

  } else {
    pools <- lapply(seq_along(comb_interventions), FUN = function(i){
      simulate(fitcov = fitcov, fitY = fitY, fitD = NA,
               yrestrictions = yrestrictions,
               compevent_restrictions = compevent_restrictions,
               restrictions = restrictions,
               outcome_name = outcome_name, compevent_name = compevent_name,
               time_name = time_name,
               intvars = comb_intvars[[i]], interventions = comb_interventions[[i]],
               int_times = comb_int_times[[i]], histvars = histvars, histvals = histvals, histories = histories,
               covparams = covparams, covnames = covnames, covtypes = covtypes,
               covpredict_custom = covpredict_custom, basecovs = basecovs, comprisk = comprisk,
               ranges = ranges,
               outcome_type = outcome_type,
               subseed = subseed, time_points = time_points,
               obs_data = obs_data, parallel = parallel,
               baselags = baselags, below_zero_indicator = below_zero_indicator,
               min_time = min_time, show_progress = FALSE, int_visit_type = int_visit_type[i], ...)
    })
  }

  nat_pool <- pools[[1]] # Natural course data
  pools <- pools[-1] # List of intervention datasets

  # Initialize results matrices
  result_ratio <- result_diff <- int_result <- rep(NA, length(pools) + 1)

  # Calculate mean outcome over all subjects at each time for natural course
  nat_result <- mean(nat_pool$Py, na.rm = TRUE)

  if (ref_int == 0){
    # Set reference intervention to the natural course
    ref_mean <- nat_result
  } else {
    # Set reference intervention as specified
    # Calculate mean outcome over all subjects at each time for this intervention
    ref_mean <- mean(pools[[ref_int]]$Py, na.rm = TRUE)
  }

  # Compile results
  int_result[1] <- nat_result
  # Calculate mean risk over all subjects at each time for all interventions other than
  # the natural course
  int_result[-1] <- sapply(pools, FUN = function(pool){mean(pool$Py, na.rm = TRUE)})
  result_ratio <- int_result / ref_mean
  result_diff <- int_result - ref_mean

  # Calculate percent intervened
  percent_intervened_res <- get_percent_intervened(pools = pools)

  if (nsamples > 0){
    if (parallel){
      cl <- prep_cluster(ncores = ncores, threads = threads , covtypes = covtypes,
                         bootstrap_option = TRUE)
      final_bs <- parallel::parLapply(cl, 1:nsamples, bootstrap_helper_with_trycatch, time_points = time_points,
                                      obs_data = obs_data_noresample, bootseeds = bootseeds,
                                      intvars = comb_intvars, interventions = comb_interventions, int_times = comb_int_times, ref_int = ref_int,
                                      covparams = covparams, covnames = covnames, covtypes = covtypes,
                                      covfits_custom = covfits_custom, covpredict_custom = covpredict_custom,
                                      basecovs = basecovs, ymodel = ymodel,
                                      histvars = histvars, histvals = histvals, histories = histories,
                                      comprisk = comprisk, compevent_model = compevent_model,
                                      yrestrictions = yrestrictions,
                                      compevent_restrictions = compevent_restrictions,
                                      restrictions = restrictions, outcome_type = outcome_type,
                                      ranges = ranges,
                                      time_name = time_name, outcome_name = outcome_name,
                                      compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                                      max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp,
                                      boot_diag = boot_diag, nsimul = nsimul, baselags = baselags,
                                      below_zero_indicator = below_zero_indicator, min_time = min_time,
                                      show_progress = FALSE, int_visit_type = int_visit_type, ...)
      parallel::stopCluster(cl)

    } else {
      if (show_progress){
        pb <- progress::progress_bar$new(total = nsamples * length(comb_interventions),
                                         clear = FALSE,
                                         format = 'Bootstrap progress [:bar] :percent, Elapsed time :elapsed, Est. time remaining :eta')
      }
      final_bs <- lapply(1:nsamples, FUN = bootstrap_helper_with_trycatch, time_points = time_points,
                         obs_data = obs_data_noresample, bootseeds = bootseeds,
                         intvars = comb_intvars, interventions = comb_interventions, int_times = comb_int_times, ref_int = ref_int,
                         covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, covpredict_custom = covpredict_custom,
                         basecovs = basecovs, ymodel = ymodel,
                         histvars = histvars, histvals = histvals, histories = histories,
                         comprisk = comprisk, compevent_model = compevent_model,
                         yrestrictions = yrestrictions,
                         compevent_restrictions = compevent_restrictions,
                         restrictions = restrictions,
                         outcome_type = outcome_type,
                         ranges = ranges,
                         time_name = time_name, outcome_name = outcome_name,
                         compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                         max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp,
                         boot_diag = boot_diag, nsimul = nsimul, baselags = baselags,
                         below_zero_indicator = below_zero_indicator, min_time = min_time,
                         show_progress = show_progress, pb = pb, int_visit_type = int_visit_type, ...)
    }
    comb_result <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$Result))
    }))
    comb_MR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultRatio))
    }))
    comb_MD <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultDiff))
    }))

    comb_result$t0 <- comb_MR$t0 <- comb_MD$t0 <- time_points

    se_result <- comb_result[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]
    se_MR <- comb_MR[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]
    se_MD <- comb_MD[, lapply(.SD, stats::sd, na.rm = TRUE), by = t0]

    if (ci_method == 'normal'){
      ci_lb_result <- t(int_result) - stats::qnorm(0.975)*se_result[,-c('t0')]
      ci_lb_MR <- t(int_result) - stats::qnorm(0.975)*se_MR[,-c('t0')]
      ci_lb_MD <- t(int_result) - stats::qnorm(0.975)*se_MD[,-c('t0')]
      ci_ub_result <- t(int_result) + stats::qnorm(0.975)*se_result[,-c('t0')]
      ci_ub_MR <- t(int_result) + stats::qnorm(0.975)*se_MR[,-c('t0')]
      ci_ub_MD <- t(int_result) + stats::qnorm(0.975)*se_MD[,-c('t0')]
    }
    if (ci_method == 'percentile') {
      ci_lb_result <- comb_result[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_lb_MR <- comb_MR[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_lb_MD <- comb_MD[, lapply(.SD, stats::quantile, probs = 0.025, na.rm = TRUE), by = t0]
      ci_ub_result <- comb_result[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
      ci_ub_MR <- comb_MR[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
      ci_ub_MD <- comb_MD[, lapply(.SD, stats::quantile, probs = 0.975, na.rm = TRUE), by = t0]
    }
  }
  if (nsamples > 0 & boot_diag){
    bootests <- comb_result
    if (!is.null(int_descript)){
      colnames(bootests)[1:(1 + length(interventions))] <- c('Natural course', int_descript)
    } else {
      colnames(bootests)[1:length(comb_interventions)] <- c('Natural course', paste('Intervention', 1:length(interventions)))
    }
    bootests[, 'Bootstrap replicate'] <- 1:nsamples
    bootcoeffs <- lapply(final_bs, "[[", 'bootcoeffs')
    bootstderrs <- lapply(final_bs, "[[", 'bootstderrs')
    bootvcovs <- lapply(final_bs, "[[", 'bootvcovs')
  } else {
    bootests <- NULL
    bootcoeffs <- NULL
    bootstderrs <- NULL
    bootvcovs <- NULL
  }

  plot_info <- get_plot_info(outcome_name = outcome_name,
                             compevent_name = compevent_name,
                             compevent2_name = compevent2_name,
                             censor_name = censor_name,
                             time_name = time_name,
                             id = id,
                             time_points = time_points,
                             covnames = covnames,
                             covtypes = covtypes,
                             nat_pool = nat_pool,
                             nat_result = nat_result,
                             comprisk = comprisk,
                             comprisk2 = comprisk2,
                             censor = censor,
                             fitC = fitC,
                             outcome_type = outcome_type,
                             obs_data = obs_data_noresample,
                             ipw_cutoff_quantile = ipw_cutoff_quantile,
                             ipw_cutoff_value = ipw_cutoff_value)
  obs_results <- plot_info$obs_results

  # Generate results table
  if (!is.null(interventions)){
    resultdf <- lapply(seq_along(int_result), function(k){
      if (nsamples > 0){
        data.table(t = time_points - 1, Intervention = k - 1,
                   EOFMean = int_result[k],
                   EOFMean_SE = se_result[[paste0('V',k)]],
                   EOFMean_CI_LL95 = ci_lb_result[[paste0('V',k)]],
                   EOFMean_CI_UL95 = ci_ub_result[[paste0('V',k)]],
                   EOFMeanRatio = result_ratio[k],
                   MR_SE = se_MR[[paste0('V',k)]],
                   MR_CI_LL95 = ci_lb_MR[[paste0('V',k)]],
                   MR_CI_UL95 = ci_ub_MR[[paste0('V',k)]],
                   EOFMeanDiff = result_diff[k],
                   MD_SE = se_MD[[paste0('V',k)]],
                   MD_CI_LL95 = ci_lb_MD[[paste0('V',k)]],
                   MD_CI_UL95 = ci_ub_MD[[paste0('V',k)]])
      } else {
        data.table(t = time_points - 1, Intervention = k - 1, OutcomeEOFMean = int_result[k],
                   EOFMeanRatio = result_ratio[k], EOFMeanDiff = result_diff[k])
      }
    })
    resultdf <- rbindlist(resultdf)
  } else {
    if (nsamples > 0){
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             EOFMean = int_result,
                             EOFMean_SE = se_MD[,-c('t0')], EOFMean_CI_LL95 = ci_lb_result[,-c('t0')],
                             EOFMean_CI_UL95 = ci_ub_result[,-c('t0')], EOFMeanRatio = result_ratio,
                             MR_SE = se_MR[,-c('t0')],
                             MR_CI_LL95 = ci_lb_MR[,-c('t0')], MR_CI_UL95 = ci_ub_MR[,-c('t0')],
                             EOFMeanDiff = result_diff,
                             EOFMD_se = se_MD[,-c('t0')],
                             MD_CI_LL95 = ci_lb_MD[,-c('t0')], MD_CI_UL95 = ci_ub_MD[,-c('t0')])
    } else {
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             OutcomeEOFMean = int_result, EOFMeanRatio = result_ratio,
                             EOFMeanDiff = result_diff)
    }
  }

  resultdf[, 'obs_risk' := c(utils::tail(obs_results[[2]], 1), rep(NA, dim(resultdf)[1] - 1))]
  obs_mean_name <- ifelse(censor, 'IP weighted mean', 'NP mean')
  if (nsamples > 0){
    colnames(resultdf) <- c("k", "Interv.", "g-form mean",
                            "Mean SE",
                            "Mean lower 95% CI", "Mean upper 95% CI",
                            "Mean ratio",
                            "MR SE", "MR lower 95% CI", "MR upper 95% CI",
                            "Mean difference",
                            "MD SE", "MD lower 95% CI", "MD upper 95% CI",
                            obs_mean_name)
    setcolorder(resultdf, c("k", "Interv.", obs_mean_name, "g-form mean",
                            "Mean SE",
                            "Mean lower 95% CI", "Mean upper 95% CI",
                            "Mean ratio",
                            "MR SE", "MR lower 95% CI", "MR upper 95% CI",
                            "Mean difference",
                            "MD SE", "MD lower 95% CI", "MD upper 95% CI"))
  } else {
    colnames(resultdf) <- c("k", "Interv.", "g-form mean", "Mean ratio",
                            "Mean difference", obs_mean_name)
    setcolorder(resultdf, c("k", "Interv.", obs_mean_name, "g-form mean",
                            "Mean ratio", "Mean difference"))
  }

  resultdf[k == max(k), '% Intervened On'] <- percent_intervened_res$percent_intervened
  resultdf[k == max(k), 'Aver % Intervened On'] <- percent_intervened_res$average_percent_intervened

  if (time_points > 1){
    fits <- fitcov
    fits[[length(fits) + 1]] <- fitY
    names(fits)[length(fits)] <- outcome_name
  } else {
    fits <- list(fitY)
  }
  if (!is.na(fitC)[[1]]){
    fits[[length(fits) + 1]] <- fitC
    names(fits)[length(fits)] <- censor_name
  }

  # Add list of coefficients for covariates, outcome variable, and competing event
  # variable (if any) to results output
  coeffs <- get_coeffs(fits = fits, fitD = NA, time_points = time_points,
                       outcome_name = outcome_name, compevent_name = compevent_name,
                       covnames = covnames)
  stderrs <- get_stderrs(fits = fits, fitD = NA, time_points = time_points,
                         outcome_name = outcome_name, compevent_name = compevent_name,
                         covnames = covnames)
  vcovs <- get_vcovs(fits = fits, fitD = NA, time_points = time_points,
                     outcome_name = outcome_name, compevent_name = compevent_name,
                     covnames = covnames)


  rmses <- lapply(seq_along(fits), FUN = rmse_calculate, fits = fits, covnames = covnames,
                  covtypes = covtypes)

  if (time_points == 1){
    rmses <- stats::setNames(rmses, outcome_name)
  } else {
    rmses <- stats::setNames(rmses, c(covnames, outcome_name))
  }


  # Create header
  header <- get_header(int_descript, sample_size, nsimul, nsamples, ref_int)

  if (sim_data_b){
    sim_data <- c(list('Natural course' = nat_pool), pools)
    if (!is.null(int_descript)){
      names(sim_data)[2:length(sim_data)] <- int_descript
    }
  } else {
    sim_data <- NA
  }
  if (!model_fits){
    fits <- NULL
  }

  res <- list(
    result = resultdf,
    coeffs = coeffs,
    stderrs = stderrs,
    vcovs = vcovs,
    rmses = rmses,
    fits = fits,
    sim_data = sim_data,
    IP_weights = obs_results$w,
    bootests = bootests,
    bootcoeffs = bootcoeffs,
    bootstderrs = bootstderrs,
    bootvcovs = bootvcovs,
    time_name = time_name,
    time_points = time_points,
    covnames = covnames,
    covtypes = covtypes,
    dt_cov_plot = plot_info$dt_cov_plot,
    dt_out_plot = plot_info$dt_out_plot,
    header = header
  )
  class(res) <- c("gformula_binary_eof", "gformula")
  return (res)
}
