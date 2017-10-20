#' Computes Akaike weights or pseudo-BMA weights for a set of models
#'
#' Returns a table with weights of a set of models, based on various
#' information criteria. Currently, \code{ictab} supports the computation of
#' Akaike weights from the \code{aic} or the \code{bic} computed on \code{lm}
#' or \code{merMod} models, as well as the computation of pseudo-BMA weights,
#' computed from the WAIC or LOOIC of \code{brmsfit} models.
#'
#' @param ic Indicates which information criterion to use. Current supported
#' information criteria include \code{aic} and \code{bic} for \code{lm} and
#' \code{merMod} models, as well as \code{WAIC} and \code{LOO} for
#' \code{brmsfit} models.
#' @param ... A set of models of class \code{lm}, \code{merMod} or
#' \code{brmsfit}.
#'
#' @importFrom stats as.formula logLik
#' @importFrom rlang dots_list f_lhs
#' @importFrom dplyr n_distinct
#' @importFrom brms WAIC LOO
#'
#' @examples
#' data(mtcars)
#' mod1 <- lm(mpg ~ cyl, mtcars)
#' mod2 <- lm(mpg ~ cyl + vs, mtcars)
#' mod3 <- lm(mpg ~ cyl + vs + I(vs^2), mtcars)
#' mod4 <- lm(mpg ~ cyl * vs, mtcars)
#' ictab(aic, mod1, mod2, mod3, mod4)
#' ictab(bic, mod1, mod2, mod3, mod4)
#'
#' \dontrun{
#' mod1 <- brm(mpg ~ cyl, mtcars)
#' mod2 <- brm(mpg ~ cyl + vs, mtcars)
#' ictab(LOO, mod1, mod2)
#' }
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{aic}}, \code{\link{bic}}
#'
#' @references Burnham, K. P., \& Anderson, D. R. (2002). Model Selection
#' and Multimodel Inference: A Practical Information-Theoretical Approach.
#' 2d ed. New York: Springer-Verlag.
#'
#' @references Burnham, K. P., \& Anderson, D. R. (2004). Multimodel
#' inference: Understanding AIC and BIC in model selection. Sociological
#' Methods and Research, 33(2), 261-304.
#'
#' @references Yao, Y. P., Vehtari, A., Simpson, D., \& Gelman, A. (2017).
#' Using stacking to average Bayesian predictive distributions.
#'
#' @export

ictab <- function(ic, ... ) {

    mods <- dots_list(...)
    modnames <- unlist(lapply(eval(substitute(alist(...) ) ), deparse) )

    check.resp <-
        lapply(mods, FUN = function(x) f_lhs(as.formula(formula(x) ) ) )

    if (n_distinct(check.resp) > 1)
        stop("\n All models should be fitted on the same data \n")

    res <- data.frame(modnames = modnames)

    if (identical(ic, aic) | identical(ic, bic) ) {

        res$ic <- unlist(lapply(mods, function(x) ic(x)[1]) )
        res$k <- unlist(lapply(mods, function(x) attr(logLik(x), "df") ) )
        res$delta_ic <- res$ic - min(res$ic)
        res$mod_lik <- exp(-0.5 * res$delta_ic)
        res$ic_wt <- res$mod_lik / sum(res$mod_lik)

    } else if (identical(ic, WAIC) | identical(ic, LOO) ) {

        res$ic <- unlist(lapply(mods, function(x) ic(x)[[3]]) )
        res$p_ic <- unlist(lapply(mods, function(x) ic(x)[[2]]) )
        res$elpd <- unlist(lapply(mods, function(x) ic(x)[[1]]) )
        res$un_wt <- exp(res$elpd - max(res$elpd) )
        res$ic_wt <- res$un_wt / sum(res$un_wt)

    }

    res <- res[rev(order(res$ic_wt) ), ]

    return(res[, -5])

}
