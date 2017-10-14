#' Computes sequential evidence ratios
#'
#' Computes sequential evidence ratios, either based on the AIC or the BIC.
#' Supported models currently include \code{lm} or \code{merMod} models.
#' When data involve repeated measures (and so multiple lines per subject),
#' a column indicating the subject "id" should be provided to the \code{id} argument.
#' If nothing is passed to the \code{id} argument, \code{seqER} will suppose
#' that there is only one observation (i.e., one line) per subject.
#'
#' @param ic Indicates whether to use the aic or the bic.
#' @param mod1 A model of class \code{lm} or \code{lmerMod}.
#' @param mod2 A model of class \code{lm} or \code{lmerMod} (of the same class of mod1).
#' @param nmin Minimum sample size from which start to compute sequential evidence ratios.
#' @param id If applicable (i.e., repeated measures), name of the "id" column of your
#' dataframe, in character string.
#'
#' @importFrom stats family formula lm
#' @importFrom lme4 lmer glmer
#' @importFrom magrittr %>%
#' @importFrom rlang f_lhs
#' @import ggplot2
#' @import utils
#'
#' @examples
#' data(mtcars)
#' mod1 <- lm(mpg ~ cyl, mtcars)
#' mod2 <- lm(mpg ~ cyl + disp, mtcars)
#' seq_mtcars <- seqER(ic = aic, mod1, mod2, nmin = 10)
#'
#' # Example with repeated measures
#' library(lme4)
#' data(sleepstudy)
#' mod1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
#' mod2 <- lmer(Reaction ~ Days + I(Days^2) + (1|Subject), sleepstudy)
#' seqER(ic = aic, mod1, mod2, nmin = 10, id = "Subject")
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{seqERboot}}
#'
#' @export

seqER <- function(ic, mod1, mod2, nmin, id = NULL) {

    if (!class(mod1) == class(mod2) ) {

        stop("Error: mod1 and mod2 have to be of the same class")

    }

    if (nmin < 10) {

        warning("nmin should usually be set above 10...")

    }

    if (class(mod1) == "lm") {

        data <- data.frame(eval(mod1$call[["data"]], envir = parent.frame() ) )

    }

    if (class(mod1) == "glmerMod" | class(mod1) == "lmerMod") {

        data <- data.frame(eval(mod1@call$data, envir = parent.frame() ) )

    }

    if (is.null(id) == TRUE) {

        id <- deparse(f_lhs(formula(mod1) ) )
        nobs <- 1
        data$ppt <- rep(seq(1, length(data[, id]), 1), each = nobs)

    } else {

        # count frequencies
        count <- data.frame(table(data[, id]) )

        # count number of observations by subject
        nobs <- max(count$Freq)

        # identify subjects with less than nobs
        a <- as.vector(count$Var1[count$Freq < nobs])

        data$ppt <-
            rep(seq(1, length (unique(data[, id]) ), 1), each = nobs)

        if (length(a) > 0) {

            # if needed, remove subjects with less than nobs
            for (i in 1:length(a) ) {

                data <- data[!data[, id] == as.numeric(a[i]), ]

            }

            warning("Different numbers of observation by subject.
                          Subjects with less than max(nobs)
                          have been removed.")
            }

    }

    startrow <- min(which(as.numeric(as.character(data$ppt) ) == nmin) )
    endrow <- nrow(data)

    for (i in seq(startrow, endrow, nobs) ) {

        maxrow <- i - 1 + nobs

        if ( (class(mod1) == "glmerMod") ) {

            mod1 <- glmer(formula(mod1),
                family = family(mod1)$family, data[1:maxrow, ])

            mod2 <- glmer(formula(mod2),
                family = family(mod2)$family, data[1:maxrow, ])

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

        tabtab <- ictab(ic, mod1, mod2)

        temp_er <- data.frame(cbind(data$ppt[i],
            tabtab$ic_wt[tabtab$modnames == "mod2"] /
                tabtab$ic_wt[tabtab$modnames == "mod1"]) )

        if (!exists("er") ) er <- temp_er else er <- rbind(er, temp_er)

        rm(temp_er)

    }

    colnames(er) <- c("ppt", "ER")
    class(er) <- c("seqER", "data.frame")

    return(er)

}

#' @export

plot.seqER <- function(x, ... ) {

    qplot(x$ppt, x$ER,
        log = "y", geom = "line",
        xlab = "Sample size",
        ylab = expression(Evidence~ ~Ratio~ ~ (ER[10]) ) ) +
        theme_bw(base_size = 12)

}
