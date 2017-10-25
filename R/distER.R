#' Simulates many sequential evidence ratios to obtain their distribution
#'
#' Simulates many sequential evidence ratios using \code{simER}, keeps the last
#' of each simulation, and plots their distribution.
#'
#' @inheritParams simER
#' @param nsims Number of experiments to simulate.
#'
#' @examples
#' \dontrun{distER(cohensd = 0.6, nmin = 20, nmax = 100, nsims = 100, ic = bic)}
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{simER}}
#'
#' @export

distER <- function(cohensd, nmin, nmax, nsims, ic = bic) {

    ER <- vector(mode = "numeric", length = nsims)

    for (i in 1:nsims) {

        ER[i] <- tail(simER(cohensd, nmin, nmax, ic, plot = FALSE), 1)

    }

    print(
        qplot(x = ER, geom = "histogram", bins = sqrt(nsims),
        alpha = 0.75, log = "x", show.legend = FALSE,
        xlab = expression(Evidence~ ~Ratio~ ~ (ER[10]) ) ) +
        theme_bw(base_size = 12)
        )

    return (ER)

}
