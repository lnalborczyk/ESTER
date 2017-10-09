#' Simulates many sequential evidence ratios, keeps the last of each simulation,
#' and plot their distribution.
#'
#' @inheritParams simER
#' @param nsims Number of experiments to simulate.
#'
#' @importFrom AICcmodavg aictab
#' @importFrom stats lm rnorm median
#' @import ggplot2
#'
#' @examples
#' library(ESTER)
#' ER <- distER(cohensd = 0.6, n = 100, nmin = 20, nsims = 100)
#'
#' @export

distER <- function (cohensd, nmin, n, nsims) {

        ER <- numeric()

        for (i in 1:nsims) {

                ER[i] <-
                        tail(simER(cohensd, n, nmin, plot = FALSE), 1) %>%
                        invisible

        }

        print(
                qplot(x = ER, geom = "histogram", bins = sqrt(nsims),
                        alpha = 0.75, log = "x", show.legend = FALSE,
                        xlab = expression(Evidence~ ~Ratio~ ~ (ER[10]) ) ) +
                        theme_bw(base_size = 12) )

        return (ER)

}
