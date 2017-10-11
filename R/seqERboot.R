#' Computes sequential evidence ratios for a given data set and bootstrapped samples
#'
#' Computes sequential evidence ratios for a given data set as well as
#' for \code{order_nb} random rearrangments of this dataset.
#'
#' @inheritParams seqER
#' @param order_nb Number of random rearrangments to evaluate.
#' @param replace If TRUE, corresponds to bootstrap with replacement.
#'
#' @importFrom stats family formula lm
#' @importFrom lme4 lmer glmer
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @examples
#' data(mtcars)
#' mod1 <- lm(mpg ~ cyl, mtcars)
#' mod2 <- lm(mpg ~ cyl + disp, mtcars)
#' seq_boot_mtcars <- seqERboot(mod1, mod2, nmin = 10, order_nb = 20, replace = FALSE)
#'
#' # Example with repeated measures
#' library(lme4)
#' data(sleepstudy)
#' mod1 <- lm(Reaction ~ Days, sleepstudy)
#' mod2 <- lm(Reaction ~ Days + I(Days^2), sleepstudy)
#' seqERboot(mod1, mod2, nmin = 10, id = "Subject", order_nb = 20, replace = TRUE)
#'
#' @export

seqERboot <- function(
    mod1, mod2, nmin, id = NULL, order_nb, replace = FALSE) {

    if (!class(mod1) == class(mod2) ) {

        stop("Error: mod1 and mod2 have to be of the same class")

    }

    if (class(mod1) == "lm") {

        data1 <- data.frame(eval(mod1$call[["data"]]) )

    }

    if (class(mod1) == "lmerMod" | class(mod1) == "glmerMod") {

        data1 <- data.frame(eval(mod1@call$data) )

    }

    order_nb <- order_nb + 1

    if (is.null(id) == TRUE) {

        id <- as.character(formula(mod1)[[2]] )
        nobs <- 1
        data1$ppt <- rep(seq(1, length(data1[, id]), 1), each = nobs)

    } else {

        # count frequencies
        count <- data.frame(table(data1[, id]) )

        # count number of observations by subject
        nobs <- max(count$Freq)

        # identify subjects with less than nobs
        a <- as.vector(count$Var1[count$Freq < nobs])

        data1$ppt <-
            rep(seq(1, length(unique(data1[, id]) ), 1), each = nobs)

        if (length(a) > 0) {

            # if needed, remove subjects with less than nobs
            for (i in 1:length(a) ) {

                data1 <- data1[!data1[, id] == as.numeric(a[i]), ]

            }

        }

    }

    for (i in (order_nb - order_nb + 2):order_nb) {

        assign(paste0("data", i),
            data1[sample(nrow(data1), replace = replace), ])

    }

    list <- ls(pattern = "data*")

    if (nobs > 1) {

        pair <- function(data) {

            data <-
                data[order(factor(data$ppt, levels = unique(data$ppt) ) ), ]

            return (data)

        }

        for (i in 1:length(list) ) {

            assign(list[i], pair(get(list[i]) ) )

            }

    }

    startrow <- nmin * nobs

    rand_er <- function(data) {

        endrow <- as.numeric(nrow(data) )

        for (i in seq(startrow, endrow, nobs) ) {

            if ( (class(mod1) == "glmerMod") ) {

                mod1 <- lme4::lmer(formula(mod1),
                    family = family(mod1)$family, data[1:i, ])

                mod2 <- lme4::lmer(formula(mod2),
                    family = family(mod2)$family, data[1:i, ])

            }

            if ( (class(mod1) == "lmerMod") ) {

                mod1 <- lme4::lmer(formula(mod1), REML = FALSE, data[1:i, ])

                mod2 <- lme4::lmer(formula(mod2), REML = FALSE, data[1:i, ])

            }

            if ( (class(mod1) == "lm") ) {

                mod1 <- lm(formula(mod1), data[1:i, ])

                mod2 <- lm(formula(mod2), data[1:i, ])

            }

            tabtab <- aictab(mod1, mod2)

            temp_er <-
                tabtab$aic_wt[tabtab$modnames == "mod2"] /
                tabtab$aic_wt[tabtab$modnames == "mod1"]

            if (!exists("er") ) er <- temp_er else er <- rbind(er, temp_er)

            rm(temp_er)

        }

        er <-
            data.frame(
                cbind(seq(nmin, max(data1$ppt), 1), er),
                row.names = NULL)

        colnames(er) <- c("ppt", "ER")

        return(er)

    }

    pb <- txtProgressBar(min = 0, max = order_nb, initial = 0, style = 3)

    for (i in 1:order_nb) {

        assign(paste0("ER", i), rand_er(get(paste0("data", i) ) ) )
        setTxtProgressBar(pb, i)

    }

    for (i in 1:order_nb) {

        ERi <- rep(paste0(paste0("ER", i) ), nrow(get(paste0("ER", i) ) ) )

        temp_er <- cbind(get(paste0("ER", i) ), ERi)
        #colnames(temp_er) <- c("ppt", "ER", "ERi")
        temp_er <- temp_er[, c(3, 1, 2)]

        if (!exists("er") ) er <- temp_er else er <- rbind(er, temp_er)

    }

    class(er) <- c("ERboot", "data.frame")

    return(er)

}

#' @export

plot.ERboot <- function(x, ... ) {

    raw <- x[, 2:3][x[, 1] == "ER1", ]

    ggplot(x, aes_string(x = "ppt", y = "ER", group = "ERi") ) +
        scale_y_log10() +
        geom_line(alpha = 0.2) +
        geom_line(aes_string(x = "ppt", y = "ER", group = NULL),
            data = raw, size = 0.75) +
        theme_bw(base_size = 12) +
        xlab("Sample size") +
        ylab(expression(Evidence~ ~Ratio~ ~ (ER[10]) ) )

}
