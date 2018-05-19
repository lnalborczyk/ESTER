
Efficient Sequential Testing with Evidence Ratios
=================================================

[![CRAN status](https://www.r-pkg.org/badges/version/ESTER)](https://cran.r-project.org/package=ESTER) [![Build Status](https://travis-ci.org/lnalborczyk/ESTER.svg?branch=master)](https://travis-ci.org/lnalborczyk/ESTER)

The `ESTER` package implements sequential testing based on evidence ratios computed from the weights of a set of models. These weights correspond to Akaike weights when based on either the Akaike Information Criterion (AIC) or the Bayesian Information Criterion (BIC), and to pseudo-BMA weights when computed from the Widely Applicable Information Criterion (WAIC).

Installation
------------

You can install the latest published version from CRAN using:

``` r
install.packages("ESTER")
```

Or the development version (recommended) from Github with:

``` r
if (!require("devtools") ) install.packages("devtools")
devtools::install_github("lnalborczyk/ESTER", dependencies = TRUE)
```

How to use the package ?
------------------------

### Model comparison

The `ictab` function takes as input a named list of models to be compared, and returns a dataframe with the given information criterion and the weight of each model.

``` r
library(ESTER)
data(mtcars)

mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + vs, mtcars)
mod3 <- lm(mpg ~ cyl * vs, mtcars)

mods <- list(mod1 = mod1, mod2 = mod2, mod3 = mod3)

ictab(mods, aic)
#>   modnames       ic k delta_ic  ic_wt
#> 1     mod1 170.1636 3   0.0000 0.7029
#> 2     mod2 172.5400 4   2.3765 0.2142
#> 3     mod3 174.4389 5   4.2753 0.0829
```

### Sequential testing

You can study the evolution of sequential ERs using the `seqtest` function.

``` r
data(mtcars)

mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)

seqtest(ic = aic, mod1, mod2, nmin = 10)
#>    ERi ppt         ER
#> 1   er  10 0.05008926
#> 2   er  11 0.07979700
#> 3   er  12 0.10363095
#> 4   er  13 0.11969544
#> 5   er  14 0.14678899
#> 6   er  15 1.10881485
#> 7   er  16 8.47867299
#> 8   er  17 6.75795190
#> 9   er  18 2.08641975
#> 10  er  19 2.19488818
#> 11  er  20 1.86450874
#> 12  er  21 2.38524035
#> 13  er  22 2.90472472
#> 14  er  23 3.39367311
#> 15  er  24 3.90918017
#> 16  er  25 1.67451190
#> 17  er  26 1.87769784
#> 18  er  27 2.10269935
#> 19  er  28 2.03674461
#> 20  er  29 2.20410125
#> 21  er  30 1.65463233
#> 22  er  31 1.73448182
#> 23  er  32 2.15556958
```

More detailed information can be found in the main vignette, available online [here](https://rawgit.com/lnalborczyk/ESTER/master/inst/doc/ESTER.html), or by typing `vignette("ESTER")` in the console.
