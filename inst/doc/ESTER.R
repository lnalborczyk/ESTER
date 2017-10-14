## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5----
library(ESTER)
ER <- simER(cohensd = 0.6, nmin = 20, n = 100, ic = bic, plot = TRUE)

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 7.5----
ER <- distER(cohensd = 0.6, nmin = 20, n = 100, nsims = 100, ic = bic)

## ----echo = TRUE, eval = TRUE--------------------------------------------
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
mod3 <- lm(mpg ~ cyl * disp, mtcars)
ictab(aic, mod1, mod2, mod3)

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
data(mtcars)
mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
plot(seqER(ic = bic, mod1, mod2, nmin = 10) )

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
plot(seqERboot(ic = bic, mod1, mod2, nmin = 10, order_nb = 20) )

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
library(lme4)
data(sleepstudy)
mod1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
mod2 <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject), sleepstudy)
plot(seqERboot(ic = bic, mod1, mod2, nmin = 10, id = "Subject", order_nb = 20) )

