#' Simulates a sequential ER for independant two-groups comparisons, as a function of sample size and Cohen's d.
#'
#' \code{simER} computes evidence ratios (ER) as a function of sample size and cohen's d
#'
#' @param cohensd Expected effect size.
#' @param n Sample size.
#' @param nmin Minimum sample size from which start to compute ER.
#' @param plot If TRUE, produces a nice plot of the evolution of the ER.
#'
#' @importFrom stats lm rnorm
#' @importFrom AICcmodavg aictab
#' @importFrom graphics abline points
#' @importFrom magrittr %>%
#' @importFrom dplyr sample_n
#'
#' @examples
#' library(ESTER)
#' ER <- simER(cohensd = 0.6, n = 100, nmin = 20, plot = TRUE)
#'
#' @export simER

simER <- function(cohensd = 0.6, n = 100, nmin = 20, plot = TRUE) {

        x <- cbind( rnorm(n = n, mean = 0, sd = 1), rep("x", n) )
        y <- cbind( rnorm(n = n, mean = mean(as.numeric(x[,1])) + cohensd, sd = 1), rep("y", n) )

        df_pop <- rbind(y, x) %>% as.data.frame
        colnames(df_pop) <- c("value", "group")
        df_pop$value <- df_pop$value %>% as.character %>% as.numeric

        df_pop <- df_pop %>% sample_n(nrow(df_pop) )

        ER_comp <- nmin %>% as.numeric

        for(i in nmin:n){

                model_1 <- lm(value ~ 1, data = df_pop[1:i, ] )
                model_2 <- lm(value ~ group, data = df_pop[1:i, ] )

                model_comp <- as.data.frame(aictab(list(model_1, model_2),
                        modnames = c("model_1", "model_2"), sort = FALSE) )

                rownames(model_comp) <- c("model_1", "model_2")

                ER_comp[i] <- model_comp["model_2", "AICcWt"] / model_comp["model_1", "AICcWt"]

        }

        if(plot == TRUE){

                suppressWarnings(plot(ER_comp, type = "l", col = "steelblue",
                        lwd = 2, xlab = "Sample size",
                        ylab = expression(Evidence~ ~Ratio~ ~(ER[10])),
                        xlim = c(nmin, n), las = 1,
                        main = paste0("Cohen's d = ", cohensd, ", ", "n = ", n),
                        log = "y", bty = "l") )

                abline(h = 1, lty = 2)
        }

        return(ER_comp[nmin:n])

        }
