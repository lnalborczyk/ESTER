#' Simulates a sequential ER for independant two-groups comparisons
#'
#' Simulates a sequential ER for independant two-groups comparisons,
#' as a function of sample size and effect size (Cohen's d).
#'
#' @param cohensd Expected effect size
#' @param nmin Minimum sample size from which start computing ERs
#' @param n Total sample size
#' @param plot If TRUE, produces a plot of the evolution of the ERs
#'
#' @importFrom magrittr %>% set_names
#' @importFrom AICcmodavg aictab
#' @importFrom stats lm rnorm
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' library(ESTER)
#' ER <- simER(cohensd = 0.6, nmin = 20, n = 100, plot = TRUE)
#'
#' @export

simER <- function(cohensd, nmin, n,  plot = TRUE) {

  x <- cbind(rnorm(n = n, mean = 0, sd = 1), rep("x", n) )
  y <- cbind(rnorm(n = n, mean = cohensd, sd = 1), rep("y", n) )

  df_pop <-
    rbind(y, x) %>%
    as.data.frame %>%
    set_names(c("value", "group") ) %>%
    mutate_(value = ~value %>% as.character %>% as.numeric) %>%
    sample_n(nrow(.) )

  ER_comp <- numeric()

  for (i in nmin:n) {

    model_1 <- lm(value ~ 1, data = df_pop[1:i, ] )
    model_2 <- lm(value ~ group, data = df_pop[1:i, ] )

    model_comp <- as.data.frame(aictab(list(model_1, model_2),
            modnames = c("model_1", "model_2"), sort = FALSE) )

    rownames(model_comp) <- c("model_1", "model_2")

    ER_comp[i] <-
      model_comp["model_2", "AICcWt"] / model_comp["model_1", "AICcWt"]

  }

  if (plot == TRUE) {

    print(
      qplot(
        nmin - 1 + seq_along(ER_comp[nmin:n]), ER_comp[nmin:n],
        log = "y", geom = "line",
        xlab = "Sample size",
        ylab = expression(Evidence~ ~Ratio~ ~ (ER[10]) ) ) +
        theme_bw(base_size = 12)
      )

  }

  return (ER_comp[nmin:n])

}
