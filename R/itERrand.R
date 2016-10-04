#' Iterative evidence ratios of original sample and random rearrangments.
#'
#' \code{itERrand} computes evidence ratios (ER) as a function of sample size,
#' for a given data set and for N random rearrangments of this dataset.
#'
#' @inheritParams itER
#' @param order_nb Number of random rearrangments to evaluate.
#' @param data Data.
#'
#' @importFrom stats aggregate family formula lm
#' @importFrom AICcmodavg aictab
#' @importFrom lme4 lmer glmer
#' @importFrom Rmisc CI
#'
#' @examples
#' library(lme4)
#' data <- sleepstudy
#' mod1 <- lm(Reaction ~ 1, data)
#' mod2 <- lm(Reaction ~ Days, data)
#' itERrand(mod1, mod2, "Subject", 10, 10, data)
#'
#' @export itERrand

itERrand <- function(mod1, mod2, samplecol, order_nb, nmin = 10, data) {

        if(!class(mod1)==class(mod2)){stop("Error: mod1 and mod2 have to be of the same class")}

        order_nb <- order_nb + 1
        data1 <- data.frame(data)

        count <- plyr::count(data1[, samplecol], 1) # count frequencies
        nobs <- as.numeric(max(count$freq)) # count number of observations by subject
        a <- as.vector(count$x[count$freq < nobs]) # identify subjects with less than n observations

        data1$ppt <- rep(seq(1, length(unique(data1[, samplecol])), 1), each = nobs)

        if(length(a)>0){

                # if needed, remove subjects with less than n observations (nobs)
                for (i in 1:length(a) ) {
                        data <- data[!data[,samplecol]==as.numeric(a[i]),]
                }
        }

        for(i in (order_nb-order_nb+2):order_nb){

                assign(paste0("data", i), data1[sample(nrow(data1)),])

        }

        list = ls(pattern = "data*")
        list = list[-which((lapply(list, nchar)<5))]

        pair <- function(data){

                data <- data[order(factor(data$ppt, levels = unique(data$ppt))),]
                return(data)
        }

        for(i in 1:length(list)){

                assign(list[i], pair(get(list[i])))

        }

        # replace this loop by a lapply ? lapply(list, pair)

        nmin <- nmin * nobs

        randER <- function(data) {

                endrow <- as.numeric(nrow(data))

                for (i in seq(nmin, endrow, nobs) ) {

                        DF <- data[1:i,]

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

                        if (!exists("ER")) ER <- tempER else ER <- rbind(ER,tempER)

                        rm(tempER)

                }

                ER <- data.frame(ER[c(2,1)])
                colnames(ER) <- c("ppt","ER")

                ER$ER <- log(ER$ER)

                return(ER)

        }

        for(i in 1:order_nb){

                assign(paste0("ER", i), randER(get(paste0("data", i))) )

        }

        for(i in (order_nb-order_nb+2):order_nb){

                ERi <- rep(paste0(paste0("ER", i)), nrow(get(paste0("ER", i))))

                tempER <- cbind(get(paste0("ER", i)), ERi)
                colnames(tempER) <- c("ppt", "ER", "ERi")
                tempER <- tempER[,c(3,1,2)]

                if (!exists("ER")) ER <- tempER else ER <- rbind(ER,tempER)

        }

        agg_ER <- data.frame(as.matrix(aggregate(ER ~ ppt, ER, CI)))

        class(agg_ER) <- c("itERrand", "data.frame")

        return(agg_ER)

}
