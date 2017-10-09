#' Simulates many sequential evidence ratios, keeps the last of each simulation,
#' and plot their distribution.
#'
#' @inheritParams simER
#' @param nSims Number of experiments to simulate.
#'
#' @importFrom stats lm rnorm t.test median IQR
#' @importFrom AICcmodavg aictab
#' @import ggplot2
#'
#' @examples
#' library(ESTER)
#' ER <- distER(cohensd = 0.6, n = 100, nmin = 20, nSims = 100)
#'
#' @export

distER <- function(cohensd, nmin, n, nSims) {

        ER <- numeric()

        for(i in 1:nSims){

                x <- rnorm(n = n, mean = 0, sd = 1)
                y <- rnorm(n = n, mean = cohensd, sd = 1)

                ER[i] <- tail(simER(cohensd, n, nmin, plot = FALSE), 1) %>% invisible

        }

        print(
                qplot(x = ER, geom = "histogram", bins = sqrt(nSims),
                        alpha = 0.75, log = "x", show.legend = FALSE,
                        xlab = expression(log-Evidence~ ~Ratio~ ~(ER[10]) ) ) +
                        theme_bw(base_size = 12) +
                        geom_vline(aes(xintercept = median(ER) ), linetype = "dashed") )

        return(ER)

}
