## ---- eval = FALSE-------------------------------------------------------
#  if(!require(devtools)){install.packages("devtools")}
#  devtools::install_github("lnalborczyk/ESTER")

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5----
library(ESTER)
ER <- simER(cohensd = 0.6, nmin = 20, n = 100, ic = aic, plot = TRUE)

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 7.5----
ER <- distER(cohensd = 0.6, nmin = 20, n = 100, nsims = 100, ic = bic)

## ----echo = FALSE, eval = TRUE, warning = FALSE, results = "hide"--------
rm(list = ls() )
library(ESTER)

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
plot(seqER(ic = aic, mod1, mod2, nmin = 10) )

## ----echo = FALSE, eval = TRUE, warning = FALSE, results = "hide"--------
library(ESTER)
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
plot(seqERboot(ic = bic, mod1, mod2, nmin = 10, order_nb = 20, replace = TRUE) )

