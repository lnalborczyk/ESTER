#' Simulates sequential testing procedures, using either evidence ratios computed
#' from WAIC or LOOIC with pseudo-BMA weights or Bayes Factors, for an independent
#' two-groups comparison
#'
#' Simulates sequential testing procedures, using either evidence ratios computed
#' from WAIC or LOOIC with pseudo-BMA weights or Bayes Factors, for an independent
#' two-groups comparison as a function of sample size and standardized mean
#' difference (Cohen's d). Bayes Factors are computed either via the Savage-Dickey
#' method (BF_SD) or via bridge sampling (BF_BS).
#'
#' @param cohensd Expected effect size
#' @param nmin Minimum sample size from which start computing ERs
#' @param nmax Maximum sample size
#' @param boundary The Evidence Ratio or the Bayes Factor (or its reciprocal)
#' to which the run is stopped as well
#' @param prior The prior that sould be placed on the slope. This will strongly
#' influence the Bayes Factor computed using the Savage-Dickey method (BF_SD).
#' @param B Number of bootstrap samples (should be dividable by cores)
#' @param cores Number of parallel processes. If cores is set to 1, no parallel framework is used.
#'
#' @return An object of class \code{data.frame}, which contains the value of the
#' evidence ratio (either WAIC_ER or LOO_ER) or the Bayes Factor (either DF_SD or BF_BS),
#' at a specific sample size n, for each simulation id.
#'
#' @importFrom stats lm rnorm update runif
#' @importFrom magrittr %>%
#' @importFrom tidyr gather_
#' @import doParallel
#' @import foreach
#' @import ggplot2
#' @import dplyr
#' @import brms
#'
#' @examples
#' \dontrun{
#' sim <- compER(cohensd = 0.8, nmin = 20, nmax = 100, boundary = 10, B = 20, cores = 4)
#' plot(sim)
#' }
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{simER}}, \code{\link{distER}}
#'
#' @export

compER <- function(
    cohensd = 0, nmin = 20, nmax = 100, boundary = 20, prior = "normal(0, 5)",
    B = 20, cores = 2) {

    if (nmin == 0) {

        stop("nmin should be a positive integer")

    }

    if (nmin > nmax) {

        stop("n should be superior to nmin")

    }

    if (nmin < 10) {

        warning("nmin should usually be set above 10...")

    }

    x <- cbind(rnorm(nmax / 2, 0, 1), rep(-0.5, nmax / 2) )
    y <- cbind(rnorm(nmax / 2, cohensd, 1), rep(0.5, nmax / 2) )

    df_pop <-
        rbind(y, x) %>%
        as.data.frame %>%
        set_names(c("value", "group") ) %>%
        sample_n(nrow(.) )

    prior1 <- set_prior("normal(0, 5)", class = "Intercept")

    prior2 <- c(
        set_prior("normal(0, 5)", class = "Intercept"),
        set_prior(prior, class = "b")
    )

    mod1 <-
        brm(value ~ 1, data = df_pop, prior = prior1,
            save_all_pars = TRUE, chains = 2, cores = cores,
            seed = sample(1e6, size = 1), refresh = 0)

    mod2.1 <-
        brm(value ~ group, data = df_pop, prior = prior2,
        save_all_pars = TRUE, chains = 2, cores = cores,
            seed = sample(1e6, size = 1), refresh = 0)

    mod2.2 <-
        brm(value ~ group, data = df_pop, prior = prior2,
        sample_prior = TRUE, chains = 2, cores = cores,
            seed = sample(1e6, size = 1), refresh = 0)

    registerDoParallel(cores = cores)

    ################################################
    # THE simulation
    ##############################

    sim <-
        foreach(
            batch = 1:getDoParWorkers(), .combine = rbind,
            .packages = c("brms", "Rcpp", "dplyr", "magrittr", "rstan", "stats")
            ) %dopar% {

        max_b <- round(B / getDoParWorkers() )
        res.counter <- 1

        # res saves the statistics at each step
        res <- matrix(NA, nrow = length(nmin:nmax) * max_b, ncol = 9,
            dimnames = list(NULL,
                c("id", "true.ES", "boundary", "prior", "n",
                    "WAIC_ER", "LOO_ER", "BF_SD", "BF_BS") ) )

        # run max_b iterations in each parallel worker
        for (b in 1:max_b) {

            # Draw a new maximum sample at each step
            x <- cbind(rnorm(nmax / 2, 0, 1), rep(-0.5, nmax / 2) )
            y <- cbind(rnorm(nmax / 2, cohensd, 1), rep(0.5, nmax / 2) )

            df_pop <-
                rbind(y, x) %>%
                as.data.frame %>%
                set_names(c("value", "group") ) %>%
                sample_n(nrow(.) )

            # res0 keeps the accumulating sample variables from this specific run
            res0 <- matrix(NA, nrow = length(nmin:nmax), ncol = ncol(res),
                dimnames = dimnames(res) )

            # increase sample size up to nmax
            for (i in nmin:nmax) {

                samp <- df_pop[1:i, ]

                mod1 <-
                    update(mod1, newdata = samp,
                        chains = 1, cores = 1, seed = sample(1e6, size = 1) )

                mod2.1 <-
                    update(mod2.1, newdata = samp,
                        chains = 1, cores = 1, seed = sample(1e6, size = 1) )

                mod2.2 <-
                    update(mod2.2, newdata = samp,
                        chains = 1, cores = 1, seed = sample(1e6, size = 1) )

                mods <- list(mod1 = mod1, mod2.1 = mod2.1)

                model_comp <- ictab(mods, WAIC)

                WAIC_ER <-
                    model_comp$ic_wt[model_comp$modnames == "mod2.1"] /
                    model_comp$ic_wt[model_comp$modnames == "mod1"]

                model_comp <- ictab(mods, LOO)

                LOO_ER <-
                    model_comp$ic_wt[model_comp$modnames == "mod2.1"] /
                    model_comp$ic_wt[model_comp$modnames == "mod1"]

                BF_SD <-
                    1 / hypothesis(mod2.2, "group = 0",
                        seed = sample(1e6, size = 1) )$hypothesis$Evid.Ratio

                # resetting the seed fixed by "bayes_factor"
                set.seed(NULL) # set.seed(Sys.time() )

                BF_BS <- bayes_factor(mod2.1, mod1, silent = TRUE)$bf

                # resetting the seed fixed by "bayes_factor"
                set.seed(NULL) # set.seed(Sys.time() )

                res0[which(nmin:nmax == i), ] <-
                    c(
                        # id is a unique id for each trajectory
                        id = batch * 10^(floor(log(max_b, base = 10) ) + 2) + b,
                        true.ES	= cohensd,
                        boundary = boundary,
                        prior = prior,
                        n = i,
                        WAIC_ER	= WAIC_ER,
                        LOO_ER = LOO_ER,
                        BF_SD = BF_SD,
                        BF_BS = BF_BS
                        )

            } # end of i

            bound_hit <- function(x) {

                if (any(x > boundary) ) {

                    first <- which(x > boundary)[1]
                    x[first:length(x)] <- boundary

                } else if (any(x < (1 / boundary) ) ){

                    first <- which(x < (1 / boundary) )[1]
                    x[first:length(x)] <- 1 / boundary

                } else{

                    x <- x

                }

                return(x)

            }

            res0[, 6:9] <-
                apply(res0[, 6:9], 2, function(x) bound_hit(x) )

            res[res.counter:(res.counter + nrow(res0) - 1), ] <- res0
            res.counter <- res.counter + nrow(res0)

        } # end of b's

        batch <- NULL

        return(res)

        } # end of %dopar%

    res <- data.frame(sim, stringsAsFactors = FALSE)
    class(res) <- c("compER", "data.frame")

    return(res)

}

#' @export

plot.compER <- function(x, log = TRUE, ... ) {

    boundary <- unique(x$boundary)
    logBoundary <- log(sort(c(boundary, 1 / boundary) ) )

    bound_hit <- function(x) {

        if (any(x > boundary) ) {

            first <- which(x > boundary)[1]
            x[first:length(x)] <- boundary

        } else if (any(x < (1 / boundary) ) ){

            first <- which(x < (1 / boundary) )[1]
            x[first:length(x)] <- 1 / boundary

        } else{

            x <- x

        }

        return(x)

    }

    bound_na <- function(x) {

        if (any(x == boundary) ) {

            first <- which(x == boundary)[1]

            if(first < length(x) ) {

                x[(first + 1):length(x)] <- NA

            }

        } else if (any(x == 1 / boundary) ) {

            first <- which(x == 1 / boundary)[1]

            if(first < length(x) ) {

                x[(first + 1):length(x)] <- NA

            }

        } else {

            x <- x

        }

        return(x)

    }

    y <-
        x %>%
        group_by(id) %>%
        mutate_(
            "WAIC_ER" = "bound_na(bound_hit(WAIC_ER))",
            "LOO_ER" = "bound_na(bound_hit(LOO_ER))",
            "BF_SD" = "bound_na(bound_hit(BF_SD))",
            "BF_BS" = "bound_na(bound_hit(BF_BS))" ) %>%
        ungroup()

    .dots <- list(~log_value == abs(logBoundary[1]) )

    final_point_boundary <-
        y %>%
        group_by(id) %>%
        gather_("index", "value", names(x)[6:9]) %>%
        mutate_("log_value" = "log(value)" ) %>%
        filter_(.dots = list(~log_value %in% logBoundary) )

    y %>%
        gather_("index", "value", names(x)[6:9]) %>%
        ggplot(aes_string(
            x = "n", y = "value",
            group = "interaction(id, index)", colour = "index") ) +
        geom_line(alpha = 0.4, size = 0.6, na.rm = TRUE) +
        geom_point(data = final_point_boundary,
            aes_string(
                x = "n", y = "value",
                group = "interaction(id, index)", colour = "index"),
            alpha = 0.6) +
        {if(log) scale_y_log10()} +
        {if(log) annotation_logticks(sides = "l")} +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor.x = element_blank(), legend.title = element_blank() ) +
        xlab("sample size") +
        ylab(expression(ER[10] ~ - ~ BF[10]) ) +
        scale_x_continuous(breaks = seq(min(y$n), max(y$n), 10 ) ) +
        annotate(
            "text",
            x = max(y$n), y = boundary,
            label = paste0(sum(final_point_boundary$value == boundary) /
                (length(unique(y$id) ) * 4) * 100, "%") ) +
        annotate(
            "text",
            x = max(y$n), y = (1 / boundary),
            label = paste0(sum(final_point_boundary$value == (1 / boundary) ) /
                    (length(unique(y$id) ) * 4) * 100, "%") )

}
