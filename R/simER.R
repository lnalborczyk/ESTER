#' Sequential ER for independant two-groups comparisons, as a function of sample size and Cohen's d.
#'
#' \code{simER} computes evidence ratios (ER) as a function of sample size and cohen's d
#'
#' @param cohensd Expected effect size.
#' @param n Sample size.
#' @param nmin Minimum sample size from which start to compute ER.
#'
#' @importFrom stats lm rnorm
#' @importFrom AICcmodavg aictab
#' @importFrom graphics abline points
#'
#' @examples
#' library(ESTER)
#' simER(cohensd = 0.2, n = 100, nmin = 20)
#'
#' @export simER

simER <- function(cohensd = 0, n = 100, nmin = 20) {

        options(scipen = 999) # disable scientific notation for numbers

        x <- cbind( rnorm(n = n, mean = 0, sd = 1), rep("x", n) )
        y <- cbind( rnorm(n = n, mean = mean(as.numeric(x[,1])) + cohensd, sd = 1), rep("y", n) )

        df_pop <- as.data.frame( rbind( y, x ) )
        colnames(df_pop) <- c("value", "group")
        df_pop$value <- as.numeric(as.character(df_pop$value) )

        df_pop <- df_pop[sample(nrow(df_pop), replace = FALSE), ]

        ########################################################
        # initialise plot and stuffs
        ############################################

        plot(1, type = "n", xlab = "sample size", ylab = "ER", xlim = c(0, n), log = "y",
                main = paste0("Cohen's d = ", cohensd, ", ", "n = ", n))

        abline(h = 1, lty = 3)

        ER_comp <- 0

        pb <- txtProgressBar(min = 0, max = n, initial = 0, style = 3) # initialise progression bar

        for(i in nmin:n){

                model_1 <- lm(value ~ 1, data = df_pop[1:i,])
                model_2 <- lm(value ~ group, data = df_pop[1:i,])

                model_comp <- as.data.frame(aictab(list(model_1, model_2),
                        modnames = c("model_1", "model_2"), sort = FALSE) )

                rownames(model_comp) <- c("model_1", "model_2")

                ER_comp <- model_comp["model_2", "AICcWt"] / model_comp["model_1", "AICcWt"]

                points(i, ER_comp, pch = 20)

                setTxtProgressBar(pb,i)

        }

        cat(paste("Final ER = ", round(ER_comp,4))) # print final ER

        }
