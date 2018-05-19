#' Computes the Akaike Information Criterion
#'
#' Computes the Akaike Information Criterion (AIC) of a model, or the second-order
#' bias correction for small samples (AICc), as suggested by
#' Burnham & Anderson (2002, 2004).
#'
#' @param mod A fitted model of class \code{lm} or \code{merMod}.
#' @param correction Should we apply the second-order correction (default to TRUE) ?
#'
#' @importFrom stats logLik nobs
#'
#' @examples
#' data(mtcars)
#' mod1 <- lm(mpg ~ cyl, mtcars)
#' aic(mod1)
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{bic}}, \code{\link{ictab}}
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

aic <- function(mod, correction = TRUE) {

    n <- nobs(mod)
    ll <- logLik(mod)[1]
    k <- attr(logLik(mod), "df")

    if (correction) {

        aic <- -2 * ll + 2 * k * (n / (n - k - 1) )

    } else {

        aic <- -2 * ll + 2 * k
    }

    return(aic)

}
