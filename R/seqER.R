#' Sequential evidence ratios.
#'
#' \code{seqER} computes evidence ratios (ER) as a function of sample size, for a given data set.
#'
#' @param mod1 A mathematical model, of class "lm" or "lmerMod".
#' @param mod2 A mathematical model, of class "lm" or "lmerMod" (of the same class of mod1).
#' @param nmin Minimum sample size from which start to compute iterative evidence ratios (ER).
#' @param samplecol Name of the subject/observation column of your dataframe, as a character vector.
#' "samplecol" has to be a column of the passed as the "data" argument of the mod1 and mod2 calls.
#'
#' @importFrom stats aggregate family formula lm
#' @importFrom AICcmodavg aictab
#' @importFrom lme4 lmer glmer
#' @importFrom utils tail
#' @importFrom graphics grid plot text
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' library(lme4)
#' data <- sleepstudy
#' mod1 <- lm(Reaction ~ 1, data)
#' mod2 <- lm(Reaction ~ Days, data)
#' seqER(mod1, mod2, samplecol = "Subject", nmin = 10)
#'
#' @export seqER

seqER <- function(mod1, mod2, samplecol, nmin) {

        if(!class(mod1)==class(mod2)){

                stop("Error: mod1 and mod2 have to be of the same class")

        }

        if(class(mod1) == "lm"){

                data <- data.frame(eval(mod1$call[["data"]]))
        }

        if(class(mod1) == "glmerMod" | class(mod1) == "lmerMod"){

                data <- data.frame(eval(mod1@call$data))
        }

        count <- plyr::count(data[, samplecol], 1) # count frequencies
        nobs <- max(count$freq) # count number of observations by subject
        a <- as.vector(count$x[count$freq < nobs]) # identify subjects with less than n observations

        data$ppt <- rep(seq(1, length(unique(data[,samplecol])), 1), each = nobs)

        if(length(a)>0){

                # if needed, remove subjects with less than n observations (nobs)
                for (i in 1:length(a) ) {
                        data <- data[!data[, samplecol]==as.numeric(a[i]), ]
                }

                warning("Different numbers of observation by subject. Subjects with less than max(nobs) have been removed.")
        }

        startrow <- min(which(as.numeric(as.character(data$ppt))==nmin)) # check the row number of the nmin
        # endrow <- max(which(as.numeric(as.character(data$ppt))==max(as.numeric(as.character(data[,samplecol])))))
        endrow <- length(data[,samplecol]) # check the row number of the last subject

        pb = txtProgressBar(min = 0, max = endrow, initial = 0, style = 3)

        for (i in seq(startrow, endrow, nobs) ) {

                maxrow <- i - 1 + nobs
                #DF <- data[1:maxrow,] # create a dataframe for each nmin + 1 sample

                if((class(mod1) == "glmerMod")){

                        mod1 <- glmer(formula(mod1), family = family(mod1)$family, data[1:maxrow,])
                        mod2 <- glmer(formula(mod2), family = family(mod2)$family, data[1:maxrow,])
                }

                if((class(mod1) == "lmerMod")){

                        mod1 <- lmer(formula(mod1), REML = FALSE, data[1:maxrow,])
                        mod2 <- lmer(formula(mod2), REML = FALSE, data[1:maxrow,])
                }

                if((class(mod1) == "lm")){

                        mod1 <- lm(formula(mod1), data[1:maxrow,])
                        mod2 <- lm(formula(mod2), data[1:maxrow,])
                }

                tabtab <- AICcmodavg::aictab(list(mod1,mod2), modnames = c("mod1","mod2"), sort = FALSE)
                tempER <- data.frame(cbind(tabtab$AICcWt[tabtab$Modnames=="mod2"] /
                                tabtab$AICcWt[tabtab$Modnames=="mod1"], data$ppt[i] ) )

                # if the merged dataset doesn't exist, create it, else, rbind it
                if (!exists("ER")) ER <- tempER else ER <- rbind(ER,tempER)

                rm(tempER)

                setTxtProgressBar(pb,i)

        }

        ER <- data.frame(ER[c(2,1)])
        colnames(ER) <- c("ppt","ER")

        class(ER) <- c("seqER", "data.frame")

        return(ER)

}

#' @S3method plot seqER
plot.seqER <- function(x, ...) {

        plot(x$ppt, x$ER, type = "l", xlab = expression(Sample~ ~size),
                ylab = expression(Evidence~ ~Ratio~ ~(ER[10])), bty = "n",
                log = "y", panel.first = grid (0, NULL, lty = 3))

        text(max(x$ppt), tail((x$ER), 1) * 1.1, as.character(round(tail(x$ER, 1), 2)))

}
