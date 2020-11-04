
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DPH

<!-- badges: start -->

[![R build
status](https://github.com/celehs/DPH/workflows/R-CMD-check/badge.svg)](https://github.com/celehs/DPH/actions)
<!-- badges: end -->

DPH provides functions to fit standard discrete proportional hazards
(DPH) model and adjusted DPH models for time-to-event data with
mismeasured
outcomes.

![](https://github.com/celehs/misclassification/blob/master/flowchart/flowchart-misclassification.jpg?raw=true)

## Example

You can install the development version from
[GitHub](https://github.com/celehs/DPH) with:

``` r
# install.packages("devtools")
devtools::install_github("celehs/DPH")
```

``` r
library(DPH)
```

Here is a simulated example which shows you how to use DPH:

``` r
sim_data
#> # A tibble: 1,000 x 4
#>     time status pred1 pred2
#>    <int>  <int> <dbl> <int>
#>  1     5      1  1.37     1
#>  2     5      1  2.18     0
#>  3     1      1  1.16     0
#>  4     1      1  3.60     0
#>  5     2      1  2.33     1
#>  6     3      1  1.18     1
#>  7     1      1  2.49     0
#>  8     1      1  2.74     0
#>  9     2      1  2.58     1
#> 10     3      1  1.69     1
#> # … with 990 more rows
```

Fit standard DPH model:

``` r
dph_fit0(
  time = sim_data$time, 
  status = sim_data$status, 
  predictors = sim_data[, -(1:2)])
#> $fit
#> 
#> Call:  stats::glm(formula = y ~ z + Z - 1, family = stats::binomial("cloglog"))
#> 
#> Coefficients:
#>       z1        z2        z3        z4        z5    Zpred1    Zpred2  
#> -2.52124  -1.86924  -1.59106  -1.49228  -1.45003   0.67569  -0.01484  
#> 
#> Degrees of Freedom: 2490 Total (i.e. Null);  2483 Residual
#> Null Deviance:       3990 
#> Residual Deviance: 2866  AIC: 2880
```

## References

  - Evaluating Biomarkers using Panel Current Status Data with
    Misclassication Errors. (Working Paper)

  - Meier, A. S., Richardson, B. A., & Hughes, J. P. (2003). Discrete
    proportional hazards models for mismeasured outcomes. Biometrics,
    59(4), 947–954. <https://doi.org/10.1111/j.0006-341x.2003.00109.x>

  - [Complementary Log-Log Model for Interval-Censored Survival
    Times](https://documentation.sas.com/?cdcId=statcdc&cdcVersion=14.2&docsetId=statug&docsetTarget=statug_logistic_examples19.htm&locale=en#statug.logistic.logx12codea)
