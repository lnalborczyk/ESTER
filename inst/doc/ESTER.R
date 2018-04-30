## ----echo = TRUE, eval = TRUE--------------------------------------------
library(ESTER)
data(mtcars)

mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)
mod3 <- lm(mpg ~ cyl * disp, mtcars)

mods <- list(m1 = mod1, m2 = mod2, m3 = mod3)

ictab(mods, aic)

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
data(mtcars)

mod1 <- lm(mpg ~ cyl, mtcars)
mod2 <- lm(mpg ~ cyl + disp, mtcars)

plot(seqtest(ic = aic, mod1, mod2, nmin = 10) )

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
plot(seqtest(ic = aic, mod1, mod2, nmin = 10, nsims = 10) )

## ----echo = TRUE, eval = TRUE, fig.align = "center", fig.height = 5, fig.width = 5, warning = FALSE, results = "hide"----
library(lme4)
data(sleepstudy)

mod1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
mod2 <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject), sleepstudy)

plot(seqtest(ic = aic, mod1, mod2, nmin = 10, id = "Subject", nsims = 10) )

