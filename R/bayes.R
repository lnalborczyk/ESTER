#' Computing likelihood and posterior probabilities of population effect sizes
#' from simulations ran with \code{simER}
#'
#' Computing the likelihood of observing an evidence ratio above or below
#' a certain threshold, given a certin effect size and a sample size
#' from simulations ran with \code{simER}. Then, using Bayes theorem, we
#' can compute the inverse probability (the posterior probability) of
#' the population effect size given a certain evidence ratio and
#' sample size (and prior probability).
#'
#' @param sim A \code{simER} object.
#' @param type Should we plot the likelihood or the posterior, or should we output a table ?
#' @param threshold A logical test (see examples)
#' @param prior Defining the prior on effect sizes (by default uniform prior)
#' @param palette Defining the color palette to be used for plots.
#'
#' @return Either a plot of the evolution of the likelihood or the posterior
#' probability, or a dataframe summarising these probabilities for each effect size
#' and sample size.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom lazyeval expr_text
#' @importFrom stats na.omit
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(ESTER)
#' # loading some simulated data
#' load(url("https://github.com/lnalborczyk/shiny/blob/master/sims_ERs.RData?raw=true") )
#' # plotting the posterior probability
#' bayes(sims, type = "posterior", threshold = ER > 10)
#'
#' # plotting the likelihood
#' bayes(sims, type = "likelihood", threshold = ER > 10)
#'
#' # changing priors
#' pr <- c(0.2, 0.2, 0.4, 0.2)
#' bayes(sims, type = "posterior", threshold = ER > 10, prior = pr)
#'
#' # obtaining the posterior probability of observing an ER > 10 if the
#' # population effect size is equal to 0 and n = 113 (with uniform priors)
#' bayes(sims, type = "table", threshold = ER > 10) %>%
#' filter(n == 113 & true.ES == 0)
#' }
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{simER}}
#'
#' @export

bayes <- function(sim, type = "posterior", threshold, prior = NULL, palette = "PuBu") {

    UseMethod("bayes")

}

#' @export

bayes.simER <- function(
    sim, type = "posterior", threshold, prior = NULL, palette = "PuBu") {

    if (!any(class(sim) == "simER") ) {

        stop("sim should be a simER object")

    }

    if (is.null(prior) ) {

        priors <-
            rep(
                1 / n_distinct(na.omit(sim$true.ES) ),
                n_distinct(na.omit(sim$true.ES) )
                )

    } else {

        priors <- prior

    }

    if (sum(priors) != 1) {

        stop("Priors on effect sizes should sum to 1")

    }

    ###################################################################
    # Using Bayes theorem to compute p(d|ER>10,n)
    #####################################################

    cond <- expr_text(threshold)

    bayesER <-
        sim %>%
        na.omit %>%
        group_by_("n", "true.ES") %>%
        dplyr::summarise_(
            # computing the likelihood of being above a particular boundary
            "likelihood" = "sum(eval(parse(text = cond) ) ) / length(ER) ") %>%
        data.frame %>%
        # there will be a problem with cases for which the likelihood is 0
        # (i.e., no simulation above the threshold) or 1 (i.e., all simulations
        # abova the threshold)... maybe replace it by very low and very high
        # probabilitites... ?
        mutate_(
            "likelihood" = "ifelse(likelihood == 0, 0.0001, likelihood)",
            "likelihood" = "ifelse(likelihood == 1, 0.9999, likelihood)",
            "true.ES" = "factor(true.ES)" ) %>%
        mutate_(
            "prior" = "rep_len(priors, nrow(.) )",
            "marginal" = "likelihood * prior") %>%
        group_by(n) %>%
        mutate_(
            "marginal" = "sum(marginal)" ) %>%
        ungroup %>%
        mutate_(
            "posterior" = "(likelihood * prior) / marginal") %>%
        data.frame

    ##############################
    # plotting
    #####################

    if (type == "likelihood") {

        ##################################################################
        # plotting the likelihood p(ER>threshold|n,d)
        ###################################################

        bayesER %>%
            mutate_(
                "true.ES" = "factor(true.ES, levels = rev(levels(true.ES) ) )"
                ) %>%
            ggplot(
                aes_string(
                    x = "n", y = "likelihood", colour = "true.ES", fill = "true.ES")
                ) +
            geom_area(size = 1, na.rm = TRUE, position = "stack") +
            scale_x_log10() +
            scale_fill_manual(
                values = rev(colorRampPalette(brewer.pal(9,
                    palette) )(n_distinct(sim$true.ES) ) ) ) +
            scale_colour_manual(
                values = rev(colorRampPalette(brewer.pal(9,
                    palette) )(n_distinct(sim$true.ES) ) ) ) +
            annotation_logticks(sides = "b") +
            theme_bw(base_size = 12) +
            xlab("Sample size") +
            ylab("p(ER,n|d)")

    } else if (type == "posterior") {

        ##################################################################
        # plotting posterior p(d|ER>threshold,n)
        ###################################################

        bayesER %>%
            mutate_(
                "true.ES" = "factor(true.ES, levels = rev(levels(true.ES) ) )"
                ) %>%
            ggplot(
                aes_string(
                    x = "n", y = "posterior", colour = "true.ES", fill = "true.ES")
                ) +
            geom_area(size = 1, na.rm = TRUE, position = "stack") +
            scale_x_log10() +
            scale_fill_manual(
                values = rev(colorRampPalette(brewer.pal(9,
                    palette) )(n_distinct(sim$true.ES) ) ) ) +
            scale_colour_manual(
                values = rev(colorRampPalette(brewer.pal(9,
                    palette) )(n_distinct(sim$true.ES) ) ) ) +
            annotation_logticks(sides = "b") +
            theme_bw(base_size = 12) +
            xlab("Sample size") +
            ylab("p(d|ER,n)")

    } else {

        bayesER

    }

}
