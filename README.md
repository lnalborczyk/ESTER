
ESTER: Efficient Sequential Testing with Evidence Ratios
========================================================

[![CRAN status](http://www.r-pkg.org/badges/version/ESTER)](https://cran.r-project.org/package=ESTER) [![Build Status](https://travis-ci.org/lnalborczyk/ESTER.svg?branch=master)](https://travis-ci.org/lnalborczyk/ESTER) <!--[![Coverage Status](https://codecov.io/github/paul-buerkner/brms/coverage.svg?branch=master)](https://codecov.io/github/paul-buerkner/brms?branch=master)-->

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

1.  **Simulation**. Given an expected effect size and a given sample size, what evolution of the evidence ratios should I reasonnably expect ?

2.  **Observed data**. When to stop recruiting participants ?

Simulation
----------

This first function runs a simulated study in which we compare two independant groups, for various effect sizes and sample sizes. The `nmin` argument serves to specify from which participant we want to start doing sequential testing (usually we recommand to avoid `nmin` &lt; 20).

``` r
library(ESTER)
ER <- simER(cohensd = 0.6, nmin = 20, n = 100, ic = aic, plot = TRUE)
```

![](README-unnamed-chunk-2-1.png)

We also can study the distribution of evidence ratios for `nSims` simulations ran with the previous function using `distER`.

``` r
ER <- distER(cohensd = 0.6, nmin = 20, n = 100, ic = aic, nsims = 100)
```

![](README-unnamed-chunk-3-1.png)

Observed data
-------------

On the other hand (and perhaps more interestingly), `ESTER` can be used to do sequential testing on your own data. You can study the evolution of sequential ERs using the `seqER` function.

``` r
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
plot(seqER(ic = bic, mod1, mod2, nmin = 10) )
```

![](README-unnamed-chunk-4-1.png)

In addition, `seqERboot` allows you to study the behavior of sequential ERs computed on your own data, along with sequential ERs computed on bootstrapped samples from your data.

``` r
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
plot(seqERboot(ic = bic, mod1, mod2, nmin = 10, order_nb = 20, replace = TRUE) )
```

![](README-unnamed-chunk-5-1.png)
