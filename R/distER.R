#' Simulates experiments and sequential evidence ratios.
#'
#' \code{distER} computes evidence ratios (ER) as a function of sample size and cohen's d
#'
#' @inheritParams simER
#' @param nSims Number of experiments to simulate.
#'
#' @importFrom stats lm rnorm t.test median IQR
#' @importFrom AICcmodavg aictab
#' @importFrom graphics abline points hist axis
#'
#' @examples
#' library(ESTER)
#' ER <- distER(cohensd = 0.6, n = 100, nmin = 20, nSims = 100)
#'
#' @export distER

distER <- function(cohensd = 0.6, n = 100, nmin = 20, nSims = 100) {

        ER <- numeric()

        for(i in 1:nSims){

                x <- rnorm(n = n, mean = 0, sd = 1)
                y <- rnorm(n = n, mean = 0 + cohensd, sd = 1)

                ER[i] <- tail(simER(cohensd, n, nmin, plot = FALSE), 1) %>% invisible

        }

        xlab <- c("[0,1/10]","[1/10,1/6]","[1/6,1/3]","[1/3,1]","[1,3]","[3,6]",
                "[6,10]",expression(paste("[10,",infinity,"]")))

        ER2 <- ER %>% cut(c(0,1/10,1/6,1/3,1,3,6,10,Inf) ) %>% table

        ER2 %>% plot(main = paste0("ER distribution with Cohen's d = ",
                cohensd, ", n = ", n, ", nSims = ", nSims),
                xlab = expression(Evidence~ ~Ratio~ ~(ER[10])),
                ylab = "Frequency", las = 1, col = "steelblue",
                lwd = 10, xaxt = "n")

        axis(1, at = 1:length(ER2), labels = xlab)

        return(ER2)

        }
