
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DPH

<!-- badges: start -->

[![R build
status](https://github.com/celehs/DPH/workflows/R-CMD-check/badge.svg)](https://github.com/celehs/DPH/actions)
<!-- badges: end -->

DPH provides functions to fit standard discrete proportional hazards
(DPH) model and adjusted DPH models for time-to-event data with
mismeasured outcomes.

## Installation

You can install the development version from
[GitHub](https://github.com/celehs/DPH) with:

``` r
# install.packages("devtools")
devtools::install_github("celehs/DPH")
```

## Example

``` r
library(DPH)
```

Here is a simulated example which shows you how to use DPH:

``` r
set.seed(1)
n <- 100
time <- sample(1:5, n, replace = TRUE)
status <- rbinom(n, 1, 0.5)
covariates <- rnorm(n)
```

Fit standard DPH model:

``` r
dph_fit0(time, status, covariates)
#> $fit
#> 
#> Call:  stats::glm(formula = y ~ z + Z - 1, family = stats::binomial("cloglog"))
#> 
#> Coefficients:
#>      z1       z2       z3       z4       z5        Z  
#> -1.9772  -2.5258  -1.5271  -0.9350  -0.8946   0.1250  
#> 
#> Degrees of Freedom: 295 Total (i.e. Null);  289 Residual
#> Null Deviance:       535.9 
#> Residual Deviance: 251.5     AIC: 263.5
```

## References

  - Evaluating Biomarkers using Panel Current Status Data with
    Misclassication Errors. (Working Paper)

  - Meier, A. S., Richardson, B. A., & Hughes, J. P. (2003). Discrete
    proportional hazards models for mismeasured outcomes. Biometrics,
    59(4), 947–954. <https://doi.org/10.1111/j.0006-341x.2003.00109.x>

  - [Complementary Log-Log Model for Interval-Censored Survival
    Times](https://documentation.sas.com/?cdcId=statcdc&cdcVersion=14.2&docsetId=statug&docsetTarget=statug_logistic_examples19.htm&locale=en#statug.logistic.logx12codea)
