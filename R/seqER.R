#' Computes evidence ratios (ER) as a function of sample size,
#' for a given data set.
#'
#' @param mod1 A mathematical model, of class "lm" or "lmerMod".
#' @param mod2 A mathematical model, of class "lm" or "lmerMod" (of the same class of mod1).
#' @param nmin Minimum sample size from which start to compute iterative evidence ratios (ER).
#' @param samplecol If applicable (e.g., repeated measures), name of the subject column of your
#' dataframe, as a character vector.
#'
#' @importFrom stats aggregate family formula lm
#' @importFrom AICcmodavg aictab
#' @importFrom lme4 lmer glmer
#' @importFrom magrittr %>%
#' @importFrom rlang f_lhs
#' @import utils
#' @import ggplot2
#'
#' @examples
#' data(mtcars)
#' mod1 <- lm(mpg ~ cyl, mtcars)
#' mod2 <- lm(mpg ~ cyl + disp, mtcars)
#' seq_mtcars <- seqER(mod1, mod2, nmin = 10)
#'
#' @export

seqER <- function(mod1, mod2, nmin, samplecol = NULL) {

        if (!class(mod1) == class(mod2) ) {

                stop("Error: mod1 and mod2 have to be of the same class")

        }

        if (class(mod1) == "lm") {

                data <- data.frame(eval(mod1$call[["data"]] ) )
        }

        if (class(mod1) == "glmerMod" | class(mod1) == "lmerMod") {

                data <- data.frame(eval(mod1@call$data) )
        }

        if (is.null(samplecol) == TRUE) {

                samplecol <- deparse(f_lhs(formula(mod1) ) )
                nobs <- 1
                data$ppt <- rep(seq(1, length(data[, samplecol]), 1),
                        each = nobs)

        }else{

                # count frequencies
                count <- plyr::count(data[, samplecol], 1)

                # count number of observations by subject
                nobs <- max(count$freq)
                # identify subjects with less than n observations
                a <- as.vector(count$x[count$freq < nobs])

                data$ppt <- rep(seq(1, length (unique(data[, samplecol]) ), 1),
                        each = nobs)

                if (length(a) > 0) {

                        # if needed, remove subjects with less than nobs
                        for (i in 1:length(a) ) {
                                data <- data[!data[, samplecol] ==
                                                as.numeric(a[i]), ]
                        }

                        warning("Different numbers of observation by subject.
                                Subjects with less than max(nobs)
                                have been removed.")
                }

        }

        # check the row number of the nmin
        startrow <- which(as.numeric(as.character(data$ppt) ) == nmin) %>% min

        # check the row number of the last subject
        endrow <- data[, samplecol] %>% length

        pb <- txtProgressBar(min = 0, max = endrow, initial = 0, style = 3)

        for (i in seq(startrow, endrow, nobs) ) {

                maxrow <- i - 1 + nobs

                if ( (class(mod1) == "glmerMod") ) {

                        mod1 <- glmer(formula(mod1),
                                family = family(mod1)$family,
                                data[1:maxrow, ])

                        mod2 <- glmer(formula(mod2),
                                family = family(mod2)$family,
                                data[1:maxrow, ])
                }

                if ( (class(mod1) == "lmerMod") ) {

                        mod1 <- lmer(formula(mod1),
                                REML = FALSE, data[1:maxrow, ])

                        mod2 <- lmer(formula(mod2),
                                REML = FALSE, data[1:maxrow, ])
                }

                if ( (class(mod1) == "lm") ) {

                        mod1 <- lm(formula(mod1), data[1:maxrow, ])

                        mod2 <- lm(formula(mod2), data[1:maxrow, ])
                }

                tabtab <- aictab(list(mod1, mod2),
                        modnames = c("mod1", "mod2"), sort = FALSE)

                tempER <- data.frame(cbind(
                        tabtab$AICcWt[tabtab$Modnames == "mod2"] /
                                tabtab$AICcWt[tabtab$Modnames == "mod1"],
                        data$ppt[i] ) )

                # if the merged dataset doesn't exist, create it, else, rbind it
                if (!exists("ER") ) ER <- tempER else ER <- rbind(ER, tempER)

                rm(tempER)

                setTxtProgressBar(pb, i)

        }

        ER <- data.frame(ER[c(2, 1)] )
        colnames(ER) <- c("ppt", "ER")

        class(ER) <- c("seqER", "data.frame")

        return(ER)

}

#' @export

plot.seqER <- function(x, ... ) {

        qplot(x$ppt, x$ER,
                log = "y", geom = "line",
                xlab = "Sample size",
                ylab = expression(Evidence~ ~Ratio~ ~ (ER[10]) ) ) +
                theme_bw(base_size = 12)

}
