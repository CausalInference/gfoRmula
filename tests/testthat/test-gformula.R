test_that("Survival outcomes", {
  id <- 'id'
  time_points <- 7
  time_name <- 't0'
  covnames <- c('L1', 'L2', 'A')
  outcome_name <- 'Y'
  outcome_type <- 'survival'
  covtypes <- c('binary', 'bounded normal', 'binary')
  histories <- c(lagged, lagavg)
  histvars <- list(c('A', 'L1', 'L2'), c('L1', 'L2'))
  covparams <- list(covmodels = c(L1 ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
                                    L3 + t0,
                                  L2 ~ lag1_A + lag_cumavg1_L1 +
                                    lag_cumavg1_L2 + L3 + t0,
                                  A ~ lag1_A + L1 + L2 + lag_cumavg1_L1 +
                                    lag_cumavg1_L2 + L3 + t0))
  ymodel <- Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0

  expect_no_error(
    gform_basic <- gformula(obs_data = basicdata_nocomp, id = id,
                            time_points = time_points,
                            time_name = time_name, covnames = covnames,
                            outcome_name = outcome_name,
                            outcome_type = outcome_type, covtypes = covtypes,
                            covparams = covparams, ymodel = ymodel,
                            intervention1.A = list(static, rep(0, time_points)),
                            intervention2.A = list(static, rep(1, time_points)),
                            int_descript = c('Never treat', 'Always treat'),
                            histories = histories, histvars = histvars,
                            basecovs = c('L3'), intcomp = c(1,2),
                            nsimul = 2501,
                            seed = 1234, nsamples = 2,
                            restrictions = list(c('L2', 'L1 == 0', simple_restriction, 0)))
  )

  compevent_name <- 'D'
  compevent_model <- D ~ A + L1 + L2 + lag1_A + lag2_A + lag3_A
  expect_no_error(
    gform_basic <- gformula(obs_data = basicdata, id = id,
                            time_points = time_points,
                            time_name = time_name, covnames = covnames,
                            outcome_name = outcome_name,
                            outcome_type = outcome_type,
                            compevent_name = compevent_name,
                            covtypes = covtypes,
                            covparams = covparams, ymodel = ymodel,
                            compevent_model = compevent_model,
                            intervention1.A = list(static, rep(0, time_points)),
                            intervention2.A = list(static, rep(1, time_points)),
                            int_descript = c('Never treat', 'Always treat'),
                            histories = histories, histvars = histvars,
                            basecovs = c('L3'), intcomp = c(1,2),
                            nsimul = 2500,
                            seed = 1234, nsamples = 2,
                            ci_method = 'normal',
                            model_fits = TRUE,
                            sim_trunc = FALSE,
                            restrictions = list(c('L2', 't0%in%c(0, 3, 6)',
                                                  carry_forward)))
  )
})


test_that("Continuous outcomes", {

  library('Hmisc')
  id <- 'id'
  time_name <- 't0'
  covnames <- c('L1', 'L2', 'A')
  outcome_name <- 'Y'
  outcome_type <- 'continuous_eof'
  covtypes <- c('categorical', 'normal', 'binary')
  histories <- c(lagged)
  histvars <- list(c('A', 'L1', 'L2'))
  covparams <- list(covmodels = c(L1 ~ lag1_A + lag1_L1 + t0 +
                                    rcspline.eval(lag1_L2, knots = c(-1, 0, 1)),
                                  L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + t0,
                                  A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + t0))
  ymodel <- Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2

  expect_no_error(
    gform_cont_eof <- gformula(obs_data = continuous_eofdata,
                             id = id, time_name = time_name,
                             covnames = covnames, outcome_name = outcome_name,
                             outcome_type = outcome_type, covtypes = covtypes,
                             covparams = covparams, ymodel = ymodel,
                             intervention1.A = list(static, rep(0, 7), int_times = 0:5),
                             intervention2.A = list(static, rep(1, 7)),
                             int_descript = c('Never treat', 'Always treat'),
                             histories = histories, histvars = histvars,
                             nsimul = 2500, seed = 1234,
                             nsamples = 2)
  )

})


test_that("Binary outcomes", {

  binary_eofdata$time_f <- ifelse(binary_eofdata$time <= 1, 0,
                                  ifelse(binary_eofdata$time <= 3, 1,
                                         ifelse(binary_eofdata$time <= 5, 2, 3)))
  binary_eofdata$time_f <- as.factor(binary_eofdata$time_f)

  outcome_type <- 'binary_eof'
  id <- 'id_num'
  time_name <- 'time'
  covnames <- c('cov1', 'cov2', 'treat', 'time_f')
  outcome_name <- 'outcome'
  histories <- c(lagged, cumavg)
  histvars <- list(c('treat', 'cov1', 'cov2'), c('cov1', 'cov2'))
  covtypes <- c('binary', 'zero-inflated normal', 'normal', 'categorical time')
  covparams <- list(covmodels = c(cov1 ~ lag1_treat + lag1_cov1 + lag1_cov2 +
                                    cov3 + time_f,
                                  cov2 ~ lag1_treat + cov1 + lag1_cov1 +
                                    lag1_cov2 + cov3 + time_f,
                                  treat ~ lag1_treat + cumavg_cov1 +
                                    cumavg_cov2 + cov3 + time_f,
                                  NA))
  ymodel <- outcome ~  treat + cov1 + cov2 + lag1_cov1 + lag1_cov2 + cov3

  expect_no_error(
    gform_bin_eof <- gformula(obs_data = binary_eofdata,
                              outcome_type = outcome_type, id = id,
                              time_name = time_name, covnames = covnames,
                              outcome_name = outcome_name, covtypes = covtypes,
                              covparams = covparams, ymodel = ymodel,
                              intervention1.treat = list(threshold, 1, Inf),
                              int_descript = 'Threshold - lower bound 1',
                              histories = histories,
                              histvars = histvars, basecovs = c("cov3"),
                              seed = 1234, parallel = TRUE,
                              nsimul = 2501, ncores = 2)
  )

})


test_that("IPCW", {

  covnames <- c('L', 'A')
  histories <- c(lagged)
  histvars <- list(c('A', 'L'))
  ymodel <- Y ~ L + A
  covtypes <- c('binary', 'normal')
  covparams <- list(covmodels = c(L ~ lag1_L + lag1_A,
                                  A ~ lag1_L + L + lag1_A))
  censor_name <- 'C'
  censor_model <- C ~ L

  expect_no_error(
    res_censor <- gformula(obs_data = censor_data, id = 'id',
                          time_name = 't0', covnames = covnames,
                          outcome_name = 'Y', outcome_type = 'survival',
                          censor_name = censor_name, censor_model = censor_model,
                          covtypes = covtypes,
                          covparams = covparams, ymodel = ymodel,
                          histories = histories, histvars = histvars,
                          seed = 1234, sim_data_b = T)
  )
})



