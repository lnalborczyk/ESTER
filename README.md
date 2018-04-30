
Efficient Sequential Testing with Evidence Ratios
=================================================

[![CRAN status](https://www.r-pkg.org/badges/version/ESTER)](https://cran.r-project.org/package=ESTER) [![Build Status](https://travis-ci.org/lnalborczyk/ESTER.svg?branch=master)](https://travis-ci.org/lnalborczyk/ESTER)

The `ESTER` package implements sequential testing based on evidence ratios computed from the weights of a set of models. These weights are being computed using either the Akaike Information Criterion (AIC) or the Bayesian Information Criterion (BIC).

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
```

More detailed information can be found in the main vignette, available online [here](https://rawgit.com/lnalborczyk/ESTER/master/inst/doc/ESTER.html), or by typing `vignette("ESTER")` in the console.
