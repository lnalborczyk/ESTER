#' Returns a table with AICs or AICc, and Akaike weights of a set of models
#'
#' Returns a table with AICs or AICc (depending on the number of parameters
#' and the number of observations, see \code{?aic}), and Akaike weights of a set of models.
#'
#' @param ... A set of models of class lm or merMod.
#'
#' @importFrom rlang dots_list f_lhs
#'
#' @examples
#' library(ESTER)
#' data(mtcars)
#' mod1 <- lm(mpg ~ cyl, mtcars)
#' mod2 <- lm(mpg ~ cyl + vs, mtcars)
#' aictab(mod1, mod2)
#'
#' @references Burnham, K. P., \& Anderson, D. R. (2002). Model Selection
#' and Multimodel Inference: A Practical Information-Theoretical Approach.
#' 2d ed. New York: Springer-Verlag.
#'
#' @references Burnham, K. P., \& Anderson, D. R. (2004). Multimodel
#' inference: Understanding AIC and BIC in model selection. Sociological
#' Methods and Research, 33(2), 261-304.
#'
#' @export

aictab <- function( ... ){

    mods <- dots_list(...)
    modnames <- unlist(lapply(eval(substitute(alist(...) ) ), deparse) )

    check.resp <- lapply(mods, FUN = function(x) f_lhs(formula(x) ) )

    if (length(unique(check.resp) ) > 1)
        stop("\nYou must use the same response variable for all models\n")

    res <- data.frame(modnames = modnames)
    res$aic <- unlist(lapply(mods, function(x) aic(x)[1]) )
    res$k <- unlist(lapply(mods, function(x) aic(x)[2]) )
    res$delta_aic <- res$aic - min(res$aic)
    res$mod_lik <- exp(-0.5 * res$delta_aic)
    res$aic_wt <- res$mod_lik / sum(res$mod_lik)

    res <- res[rev(order(res$aic_wt) ), ]

    return(res[, -5])

}
