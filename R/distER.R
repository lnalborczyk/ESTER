#' Simulates experiments and sequential evidence ratios.
#'
#' \code{distER} computes evidence ratios (ER) as a function of sample size and cohen's d
#'
#' @inheritParams simER
#' @param nSims Number of experiments to simulate.
#'
#' @importFrom stats lm rnorm t.test median
#' @importFrom AICcmodavg aictab
#' @importFrom graphics abline points hist
#'
#' @examples
#' library(ESTER)
#' distER(cohensd = 0.4, n = 100, nmin = 20, nSims = 100)
#'
#' @export distER

distER <- function(cohensd = 0, n = 100, nmin = 20, nSims = 100) {

        options(scipen = 999) # disable scientific notation for numbers

        ER <- numeric(nSims) # set up empty variable to store all simulated ER

        for(i in 1:nSims){

                x <- rnorm(n = n, mean = 0, sd = 1) # produce N simulated participants
                y <- rnorm(n = n, mean = 0 + cohensd, sd = 1) # produce N simulated participants
                z <- t.test(x,y) # perform a t-test for independant samples

                ER[i] <- tail(simER(cohensd, n, nmin, plot = FALSE),1)

        }

        xlim = c(0, mean(ER) + IQR(ER)/2)

        hist(ER, breaks = 1000, border = FALSE, main = paste0("ER distribution with Cohen's d = ",
                cohensd, ", n = ", n, ", nSims = ", nSims),
                xlab = expression(Evidence~ ~Ratio~ ~(ER[10])),
                xlim = xlim, ylim = c(0, nSims), col = "steelblue")

        rug(ER, side = 3)

        text(median(ER), nSims/2, as.character(paste0("median = ", round(median(ER), 3))))

        return(ER)

        }
