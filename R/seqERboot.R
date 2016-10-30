#' Sequential evidence ratios of original sample and random rearrangments.
#'
#' \code{seqERboot} computes evidence ratios (ER) as a function of sample size,
#' for a given data set and for N random rearrangments of this dataset.
#'
#' @inheritParams seqER
#' @param order_nb Number of random rearrangments to evaluate.
#' @param replace If TRUE, corresponds to bootstrap with replacement.
#'
#' @importFrom stats aggregate family formula lm loess
#' @importFrom AICcmodavg aictab
#' @importFrom lme4 lmer glmer
#' @importFrom graphics lines
#'
#' @examples
#' library(lme4)
#' data <- sleepstudy
#' mod1 <- lm(Reaction ~ 1, data)
#' mod2 <- lm(Reaction ~ Days, data)
#' seqERboot(mod1, mod2, samplecol = "Subject", order_nb = 10, nmin = 10, replace = FALSE)
#'
#' @export seqERboot

seqERboot <- function(mod1, mod2, samplecol, order_nb, nmin = 10, replace = FALSE) {

        if(!class(mod1)==class(mod2)){stop("Error: mod1 and mod2 have to be of the same class")}

        if(class(mod1) == "lm"){

                data1 <- data.frame(eval(mod1$call[["data"]]))
        }

        if(class(mod1) == "lmerMod" | class(mod1) == "glmerMod"){

                data1 <- data.frame(eval(mod1@call$data))

        }

        order_nb <- order_nb + 1

        count <- plyr::count(data1[, samplecol], 1) # count frequencies
        nobs <- as.numeric(max(count$freq)) # count number of observations by subject
        a <- as.vector(count$x[count$freq < nobs]) # identify subjects with less than n observations

        data1$ppt <- rep(seq(1, length(unique(data1[, samplecol])), 1), each = nobs)

        if(length(a)>0){

                # if needed, remove subjects with less than n observations (nobs)
                for (i in 1:length(a) ) {
                        data1 <- data1[!data1[,samplecol]==as.numeric(a[i]), ]
                }
        }

        for(i in (order_nb-order_nb+2):order_nb){

                assign(paste0("data", i), data1[sample(nrow(data1), replace = replace), ])

        }

        list = ls(pattern = "data*")

        pair <- function(data){

                data <- data[order(factor(data$ppt, levels = unique(data$ppt))), ]
                return(data)
        }

        #list <- list(lapply(ls(pattern="data*"), function(x) get(x)))
        #rapply(list, pair, how = "replace")

        for(i in 1:length(list)){

                assign(list[i], pair(get(list[i])) )

        }

        startrow <- nmin * nobs

        randER <- function(data) {

                endrow <- as.numeric(nrow(data))

                for (i in seq(startrow, endrow, nobs) ) {

                        #DF <- data[1:i,]

                        if((class(mod1) == "glmerMod")){

                                mod1 <- lme4::lmer(formula(mod1), family = family(mod1)$family, data[1:i,])
                                mod2 <- lme4::lmer(formula(mod2), family = family(mod2)$family, data[1:i,])

                        }

                        if((class(mod1) == "lmerMod")){

                                mod1 <- lme4::lmer(formula(mod1), REML = FALSE, data[1:i,])
                                mod2 <- lme4::lmer(formula(mod2), REML = FALSE, data[1:i,])

                        }

                        if((class(mod1) == "lm")){

                                mod1 <- lm(formula(mod1), data[1:i,])
                                mod2 <- lm(formula(mod2), data[1:i,])

                        }


                        tabtab <- AICcmodavg::aictab(list(mod1,mod2), modnames = c("mod1","mod2"), sort = FALSE)
                        tempER <- data.frame(cbind(tabtab$AICcWt[tabtab$Modnames=="mod2"] /
                                        tabtab$AICcWt[tabtab$Modnames=="mod1"]) )

                        if (!exists("ER")) ER <- tempER else ER <- rbind(ER,tempER)

                        rm(tempER)

                }

                if(class(data1[,samplecol])=="factor") nbppt <- as.numeric(length(levels(data1[,samplecol])))
                if(class(data1[,samplecol])=="integer") nbppt <- as.numeric(length(unique(data1[,samplecol])))
                if(class(data1[,samplecol])=="numeric") nbppt <- as.numeric(length(unique(data1[,samplecol])))

                ER <- data.frame(cbind(ER, seq(nmin, nbppt, 1) ) )
                ER <- data.frame(ER[c(2,1)])
                colnames(ER) <- c("ppt","ER")

                return(ER)

        }

        pb = txtProgressBar(min = 0, max = order_nb, initial = 0, style = 3)

        for(i in 1:order_nb){

                assign(paste0("ER", i), randER(get(paste0("data", i))) )

                setTxtProgressBar(pb,i)

        }

        for(i in 1:order_nb){

                ERi <- rep(paste0(paste0("ER", i)), nrow(get(paste0("ER", i))))

                tempER <- cbind(get(paste0("ER", i)), ERi)
                colnames(tempER) <- c("ppt", "ER", "ERi")
                tempER <- tempER[,c(3,1,2)]

                if (!exists("ER")) ER <- tempER else ER <- rbind(ER,tempER)

        }

        agg_ER <- data.frame(as.matrix(aggregate(ER ~ ppt, ER, mean)))

        class(ER) <- c("itER", "data.frame")
        class(agg_ER) <- c("itERrand", "data.frame")

        return_list <- list("ER" = ER, "agg_ER" = agg_ER)
        class(return_list) <- c("ERlist", "list")

        return(return_list)

}

#' @S3method plot ERlist
plot.ERlist <- function(x, ...) {

        ylim <- c(min(x$ER$ER)/1.1, max(x$ER$ER)*1.1 )

        plot(x$ER$ppt[x$ER$ERi=="ER1"], x$ER$ER[x$ER$ERi=="ER1"], type = "l",
                lwd = 1.5, xlab = expression(Sample~ ~size),
                ylab = expression(Evidence~ ~Ratio~ ~(ER[10])),
                bty = "n", log = "y", ylim = ylim,
                panel.first = grid(0, NULL, lty = 3) )

        for(i in 2:nlevels(x$ER$ERi)){

                lines(loess( x$ER$ER[x$ER$ERi==as.character(paste0("ER", i))] ~
                                x$ER$ppt[x$ER$ERi==as.character(paste0("ER", i))]  ),
                        lwd = 0.8, col = "grey" )

        }

        options( digits = 3 )

        text(max(x$ER$ppt[x$ER$ERi=="ER1"]), tail((x$ER$ER[x$ER$ERi=="ER1"]), 1) * 1.1,
                as.character(round(tail(x$ER$ER[x$ER$ERi=="ER1"], 1), 2)))

}

