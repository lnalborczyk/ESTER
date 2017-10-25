#' Simulates a sequential testing with evidence ratios for independent two-groups comparisons
#'
#' Simulates a sequential testing with evidence ratios for independent two-groups
#' comparisons, as a function of sample size and standardized mean difference
#' (Cohen's d).
#'
#' @param ic Indicates whether to use the aic or the bic.
#' @param cohensd Expected effect size
#' @param nmin Minimum sample size from which start computing ERs
#' @param nmax Total sample size
#' @param plot If TRUE, produces a plot of the evolution of the ERs
#'
#' @importFrom magrittr %>% set_names
#' @importFrom stats lm rnorm update
#' @importFrom brms brm WAIC LOO
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' simER(cohensd = 0.6, nmin = 20, nmax = 200, ic = aic, plot = TRUE)
#' simER(cohensd = 0, nmin = 20, nmax = 200, ic = bic, plot = TRUE)
#'
#' \dontrun{simER(cohensd = 0, nmin = 20, nmax = 60, ic = LOO, plot = TRUE)}
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{distER}}
#'
#' @export

simER <- function(cohensd, nmin, nmax, ic = bic, plot = TRUE) {

    if (nmin == 0) {

        stop("nmin should be a positive integer")

    }

    if (nmin > nmax) {

        stop("n should be superior to nmin")

    }

    if (nmin < 10) {

        warning("nmin should usually be set above 10...")

    }

    x <- cbind(rnorm(nmax, 0, 1), rep("x", nmax) )
    y <- cbind(rnorm(nmax, cohensd, 1), rep("y", nmax) )

    df_pop <-
        rbind(y, x) %>%
        as.data.frame %>%
        set_names(c("value", "group") ) %>%
        mutate_(value = ~value %>% as.character %>% as.numeric) %>%
        sample_n(nrow(.) )

    ER_comp <- numeric()

    if (identical(ic, aic) | identical(ic, bic) ) {

        for (i in nmin:nmax) {

            mod1 <- lm(value ~ 1, data = df_pop[1:i, ] )
            mod2 <- lm(value ~ group, data = df_pop[1:i, ] )

            model_comp <- ictab(list(mod1 = mod1, mod2 = mod2), ic)

            ER_comp[i] <-
                model_comp$ic_wt[model_comp$modnames == "mod2"] /
                model_comp$ic_wt[model_comp$modnames == "mod1"]

        }

    } else if (identical(ic, WAIC) | identical(ic, LOO) ) {

        mod1 <- brm(value ~ 1, data = df_pop[1:nmin, ] )
        mod2 <- brm(value ~ group, data = df_pop[1:nmin, ] )

        model_comp <- ictab(list(mod1 = mod1, mod2 = mod2), ic)

        ER_comp[nmin] <-
            model_comp$ic_wt[model_comp$modnames == "mod2"] /
            model_comp$ic_wt[model_comp$modnames == "mod1"]

        for (i in (nmin + 1):nmax) {

            mod1 <- update(mod1, newdata = df_pop[1:i, ] )
            mod2 <- update(mod2, newdata = df_pop[1:i, ] )

            model_comp <- ictab(list(mod1 = mod1, mod2 = mod2), ic)

            ER_comp[i] <-
                model_comp$ic_wt[model_comp$modnames == "mod2"] /
                model_comp$ic_wt[model_comp$modnames == "mod1"]

        }

    }

    if (plot == TRUE) {

        print(
            qplot(
                nmin - 1 + seq_along(ER_comp[nmin:nmax]), ER_comp[nmin:nmax],
                log = "y", geom = "line",
                xlab = "Sample size",
                ylab = expression(Evidence~ ~Ratio~ ~ (ER[10]) ) ) +
                theme_bw(base_size = 12)
            )

    }

    return (ER_comp[nmin:nmax])

}
