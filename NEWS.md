Package Updates
---------------

### Changes in Version 1.0.0 (TBD)

-   Added option for users to specify censoring models to compute
    inverse probability weights for estimating the natural course means
    / risk from the observed data
-   Fixed an error in calculating the means of the time-varying
    covariates under the natural course for survival outcomes
-   Fixed an error in calculating the observed risk estimates when
    competing events are not treated like censoring events
-   Fixed an error in calculating the g-formula survival estimates when
    competing events are not treated like censoring events
-   For categorical time-varying covariates, the
    `plot.gformula_survival()`, `gformula_continuous_eof()`, and
    `gformula_binary_eof()` functions now display the nonparametric/IP
    weighted and parametric g-formula estimates of the probability of
    observing each level of the covariate. Previously, these functions
    displayed the counts of categorical variables.

### Changes in Version 0.3.2 (2021-07-13)

-   Updated computation of (lagged) cumulative averages to use the
    recursive formula. There should be a noticeable improvement in the
    computation time when using several (lagged) cumulative average
    terms and when the number of time points is large.
-   Fixed an error for covariates of type `truncated normal` (Thanks to
    @publichealthstudent)
-   Updates to the documentation

### Changes in Version 0.3.1 (2020-03-22)

-   Fixed error in the `coef.gformula()` example

### Changes in Version 0.3.0 (2020-01-30)

-   Added wrapper function called `gformula()` for the
    `gformula_survival()`, `gformula_continuous_eof()`, and
    `gformula_binary_eof()` functions. Users should now use the more
    general `gformula()` function to apply the g-formula.
-   Added option for users to specify the values for lags at
    pre-baseline times by including rows at time -1, -2, …, -i.
-   Added an example data set called `continuous_eofdata_pb`, which
    illustrates how to prepare a data set with pre-baseline times
-   Added option for users to pass in “control parameters” (e.g.,
    maximum number of iterations, maxit, in glm.control) when fitting
    models for time-varying covariates via the `covparams$control`
    argument. (Thanks to @jerzEG for the suggestion)
-   Added option for users to access the fitted models for the
    time-varying covariates, outcome, and competing event (if
    applicable). See `model_fits` argument of the `gformula()` function
-   Added simulated data under the natural course to the `sim_data`
    component of the output of the `gformula()` function
-   Added a progress bar for the number of bootstrap samples completed.
    See the `show_progress` argument of the `gformula()` function for
    further details
-   Added `summary()`, `coef()`, and `vcov()` S3 methods for objects of
    class ‘gformula’
-   Added argument `fits` in the `print.gformula_survival()`,
    `print.gformula_continuous_eof()`, and `print.gformula_binary_eof()`
    functions. Added argument `all_times` in the
    `print.gformula_survival()` function
-   Fixed minor bug in the `lagavg()` function
-   Fixed bug occuring when not using lags of the intervention
    variable(s)
-   Fixed bug occuring in the truncation beyond covariate ranges.
    (Thanks to Louisa Smith)
-   Updates to the documentation

### Changes in Version 0.2.1 (2019-08-24)

-   First version released on CRAN
    (<https://CRAN.R-project.org/package=gfoRmula>)
-   Updates to the documentation

### Changes in Version 0.2.0 (2019-08-22)

-   Removed `example_intervention1()`, `example_intervention2()`, and
    `visit_sum_orig()`, as these functions are not used internally and
    users should not directly apply them
-   Removed export of `visit_sum()` and `natural()`, as these functions
    are used internally and users should not directly apply them
-   Updates to the documentation

### Changes in Version 0.1.1 (2019-08-21)

-   Minor updates to the documentation

### Changes in Version 0.1.0 (2019-08-17)

-   First version released on GitHub
    (<https://github.com/CausalInference/gfoRmula>)
