
Efficient Sequential Testing with Evidence Ratios
=================================================

[![CRAN status](https://www.r-pkg.org/badges/version/ESTER)](https://cran.r-project.org/package=ESTER) [![Build Status](https://travis-ci.org/lnalborczyk/sticer.svg?branch=master)](https://travis-ci.org/lnalborczyk/ESTER)

The `ESTER` package implements sequential testing based on evidence ratios computed from the weights of a set of models. These weights are being computed using either the Akaike Information Criterion (AIC) or the Bayesian Information Criterion (BIC).

Installation
------------

You can install the latest published version from CRAN using:

``` r
install.packages("ESTER")
```

Or the development version from Github with:

``` r
if (!require("devtools") ) install.packages("devtools")
devtools::install_github("lnalborczyk/ESTER", dependencies = TRUE)
```

How to use the package ?
------------------------

### Model comparison

The `ictab` function...

``` r
library(ESTER)
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + vs, mtcars)
mod3 <- lm(mpg ~ cyl * vs, mtcars)
mods <- list(mod1 = mod1, mod2 = mod2, mod3 = mod3)
ictab(mods, aic)
ictab(mods, bic)
```

Sequential testing
------------------

On the other hand (and perhaps more interestingly), `ESTER` can be used to do sequential testing on your own data. You can study the evolution of sequential ERs using the `seqER` function.

``` r
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
seqER(ic = aic, mod1, mod2, nmin = 10)
```

In addition, `seqER` allows you to study the behaviour of sequential ERs computed on your own data, along with sequential ERs computed on permutation samples. This feature might be useful to study to what extent the evolution of evidence ratios you observed on the original sample is dependent to the order of the observations.

``` r
seqER(ic = aic, mod1, mod2, nmin = 10, nsims = 10)
```

More detailed information can be found in the main vignette, available online [here](https://rawgit.com/lnalborczyk/ESTER/master/inst/doc/ESTER.html), or by typing `vignette("ESTER")` in the console.
