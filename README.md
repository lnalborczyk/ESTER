ESTER: Efficient Sequential Testing with Evidence Ratios
===

[![CRAN status](http://www.r-pkg.org/badges/version/ESTER)](https://cran.r-project.org/package=ESTER) [![Build Status](https://travis-ci.org/lnalborczyk/ESTER.svg?branch=master)](https://travis-ci.org/lnalborczyk/ESTER)

## Installation

Development version from Github:

``` r
if(!require(devtools)){install.packages("devtools")}
devtools::install_github("lnalborczyk/ESTER")
```
## Different questions

1. **Simulation**. Given an expected effect size and sample size, what ER evolution should I reasonnably expect ?

2. **Observed data**. When to stop recruiting participants ?

## 1. Simulation

This first function runs a simulated study in which we compare two independant groups, for various effect size and sample size. The `nmin` argument serves to specify from which participant we want to start doing sequential testing (usually we recommand to avoid `nmin` < 20).

```r
library(ESTER)
ER <- simER(cohensd = 0.6, nmin = 20, n = 100, plot = TRUE)
```

We also can study the distribution of evidence ratios for `nSims` simulations runned with the previous function using `distER`, where the plotted vertical dashed line represents the median of the ERs distribution.

```r
ER <- distER(cohensd = 0.6, nmin = 20, n = 100, nSims = 100)
```

## 2. Observed data

On the other hand (and perhaps more interestingly), `ESTER` can be used to do sequential testing on your own data. You can study the evolution of sequentials ERs using the `seqER` function.

```r
rm(list = ls() )
library(ESTER)
```

```r
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
plot(seqER(mod1, mod2, nmin = 10) )
```

In addition, `seqERboot` allows you to study the behavior of sequential ERs computed on your own data, along with sequential ERs computed on bootstrapped samples from your data.

```r
library(ESTER)
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
```

```r
plot(seqERboot(mod1, mod2, nmin = 10, order_nb = 20, replace = TRUE) )
```

Note that under the hood, `ESTER` uses the `AICcmodavg` package (available on [CRAN](https://cran.r-project.org/web/packages/AICcmodavg/index.html); Mazerolle, 2016) to compute AICc and Akaike Weights.

For more information and theoretic backgroun, read the [manual](https://rawgit.com/lnalborczyk/ESTER/master/vignettes/ESTER.html).
