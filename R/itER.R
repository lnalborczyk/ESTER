#' Iterative evidence ratios.
#'
#' \code{itER} computes evidence ratios (ER) as a function of sample size, for a given data set.
#'
#' @param mod1 A mathematical model, of class "lm" or "lmerMod".
#' @param mod2 A mathematical model, of class "lm" or "lmerMod" (of the same class of mod1).
#' @param nmin Minimum sample size from which start to compute iterative evidence ratios (ER).
#' @param samplecol Name of the subject/observation column of your dataframe, as a character vector.
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from the model specification (mod1 and mod2).
#'
#' @importFrom stats aggregate family formula lm
#' @importFrom AICcmodavg aictab
#' @importFrom lme4 lmer glmer
#'
#' @examples
#' library(lme4)
#' data <- sleepstudy
#' mod1 <- lm(Reaction ~ 1, data)
#' mod2 <- lm(Reaction ~ Days, data)
#' itER(mod1, mod2, "Subject", 10, data)
#'
#' @export itER

itER <- function(mod1, mod2, samplecol, nmin, data) {

        if(!class(mod1)==class(mod2)){

                stop("Error: mod1 and mod2 have to be of the same class")

                }

        data <- data.frame(data)

        count <- plyr::count(data[, samplecol], 1) # count frequencies
        nobs <- max(count$freq) # count number of observations by subject
        a <- as.vector(count$x[count$freq < nobs]) # identify subjects with less than n observations

        data$ppt <- rep(seq(1, length(unique(data[,samplecol])), 1), each = nobs)

        if(length(a)>0){

                # if needed, remove subjects with less than n observations (nobs)
                for (i in 1:length(a) ) {
                        data <- data[!data[, samplecol]==as.numeric(a[i]), ]
                }
        }

        startrow <- min(which(as.numeric(as.character(data$ppt))==nmin)) # check the row number of the nmin
        # endrow <- max(which(as.numeric(as.character(data$ppt))==max(as.numeric(as.character(data[,samplecol])))))
        endrow <- length(data[,samplecol]) # check the row number of the last subject

        for (i in seq(startrow, endrow, nobs) ) {

                maxrow <- i - 1 + nobs
                DF <- data[1:maxrow,]

                if((class(mod1) == "glmerMod")){

                        mod1 <- lmer(formula(mod1), REML = FALSE, family = family(mod1)$family, DF)
                        mod2 <- lmer(formula(mod2), REML = FALSE, family = family(mod2)$family, DF)
                }

                if((class(mod1) == "lmerMod")){

                        mod1 <- lmer(formula(mod1), REML = FALSE, DF)
                        mod2 <- lmer(formula(mod2), REML = FALSE, DF)
                }

                if((class(mod1) == "lm")){

                        mod1 <- lm(formula(mod1), DF)
                        mod2 <- lm(formula(mod2), DF)
                }

                tabtab <- aictab(list(mod1,mod2), modnames = c("mod1","mod2"), sort = FALSE)
                tempER <- data.frame(cbind(tabtab$AICcWt[tabtab$Modnames=="mod2"] /
                                tabtab$AICcWt[tabtab$Modnames=="mod1"], data$ppt[i] ) )

                # if the merged dataset doesn't exist, create it, else, rbind it
                if (!exists("ER")) ER <- tempER else ER <- rbind(ER,tempER)

                rm(tempER)

        }

        ER <- data.frame(ER[c(2,1)])
        colnames(ER) <- c("ppt","ER")

        ER$ER <- log(ER$ER)

        return(ER)

}
