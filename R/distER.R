#' Sequential ER for independant two-groups comparisons, as a function of sample size and Cohen's d.
#'
#' \code{distER} computes evidence ratios (ER) as a function of sample size and cohen's d
#'
#' @inheritParams simER
#' @param nSims Number of experiments to simulate.
#'
#' @importFrom stats lm rnorm t.test
#' @importFrom AICcmodavg aictab
#' @importFrom graphics abline points hist
#'
#' @examples
#' library(ESTER)
#' distER(cohensd = 0.2, n = 100, nmin = 20, nSims = 20)
#'
#' @export distER

distER <- function(cohensd = 0, n = 100, nmin = 20, nSims = 20) {

        options(scipen = 999) # disable scientific notation for numbers

        ER <- numeric(nSims) # set up empty variable to store all simulated ER

        for(i in 1:nSims){

                x <- rnorm(n = n, mean = 0, sd = 1) # produce N simulated participants
                y <- rnorm(n = n, mean = 0 + cohensd, sd = 1) # produce N simulated participants
                z <- t.test(x,y) # perform a t-test for independant samples

                ER[i] <- simER(cohensd, n, nmin)

        }


        hist(ER, breaks = 20, border = FALSE, main = "ER distribution",
                xlim = c(0, max(ER)), ylim = c(0, nSims), col = "steelblue")

        return(ER)

        }
