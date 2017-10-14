
ESTER: Efficient Sequential Testing with Evidence Ratios
========================================================

[![CRAN status](http://www.r-pkg.org/badges/version/ESTER)](https://cran.r-project.org/package=ESTER) [![Build Status](https://travis-ci.org/lnalborczyk/ESTER.svg?branch=master)](https://travis-ci.org/lnalborczyk/ESTER)

The `ESTER` package implements sequential testing based on evidence ratios computed from the Akaike weights of a set of models. These weights are being computed using either the Akaike Information Criterion (AIC) or the Bayesian Information Criterion (BIC).

Installation
------------

You can install ESTER from github with:

``` r
if(!require(devtools)){install.packages("devtools")}
devtools::install_github("lnalborczyk/ESTER")
```

Different questions
-------------------

1.  **Simulation**. Given an expected effect size and a given sample size, what evolution of evidence ratios should I reasonnably expect ?

2.  **Observed data**. When to stop recruiting participants ?

Simulation
----------

This first function runs a simulated study in which we compare two independant groups, for various effect sizes and sample sizes. The `nmin` argument serves to specify from which participant we want to start doing sequential testing (we usually recommand to avoid `nmin` &lt; 10).

``` r
library(ESTER)
simER(cohensd = 0.6, nmin = 20, n = 100, ic = aic, plot = TRUE)
```

We also can study the distribution of evidence ratios for `nSims` simulations ran with the previous function using `distER`.

``` r
distER(cohensd = 0.6, nmin = 20, n = 100, ic = aic, nsims = 100)
```

Observed data
-------------

On the other hand (and perhaps more interestingly), `ESTER` can be used to do sequential testing on your own data. You can study the evolution of sequential ERs using the `seqER` function.

``` r
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
seqER(ic = bic, mod1, mod2, nmin = 10)
```

In addition, `seqERboot` allows you to study the behavior of sequential ERs computed on your own data, along with sequential ERs computed on permutation samples. This feature might be useful to study to what extent the evolution of evidence ratios you observed on the original sample is dependant to the order of the observations.

``` r
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
seqERboot(ic = bic, mod1, mod2, nmin = 10, order_nb = 20)
```

More detailed information can be found by typing `vignette("ESTER")` in the console.
