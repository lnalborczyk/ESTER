#' Analysing the results of simulations ran with \code{simER} or \code{compER}
#'
#' Analysing the results of simulations ran with \code{simER} or \code{compER}.
#' It computes the average sample number (ASN) at which the boundary is attained
#' (either the lower or the upper one), the percentage of hits of the lower
#' boundary as well as hits of the upper boundary, and the percentage of
#' trajectories that did not hit none of the boundaries.
#'
#' @param sim A \code{simER} or a \code{compER} object.
#'
#' @return An object of class \code{data.frame}, which contains the average
#' sample number (ASN) at which the boundary is attained (either the lower or
#' the upper one), the percentage of hits of the lower boundary as well as hits
#' of the upper boundary, and the percentage of trajectories that did not hit
#' none of the boundaries (and thus end at nmax).
#'
#' @importFrom stats sd na.omit
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(ESTER)
#' sim <- simER(cohensd = 0.8, nmin = 20, nmax = 100, boundary = 10, nsims = 100, ic = bic)
#' analysER(sim)
#' }
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{simER}}, \code{\link{compER}}
#'
#' @export

analysER <- function(sim) {

    UseMethod("analysER")

}

#' @export

analysER.simER <- function(sim) {

    if (!any(class(sim) %in% c("simER", "compER") ) ) {

        stop("sim should be a simER or a compER object")

    }

    boundary <- as.numeric(unique(sim$boundary) )
    logboundary <- log(sort(c(boundary, 1 / boundary) ) )

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

            if (first < length(x) ) {

                x[(first + 1):length(x)] <- NA

            }

        } else if (any(x == 1 / boundary) ) {

            first <- which(x == 1 / boundary)[1]

            if (first < length(x) ) {

                x[(first + 1):length(x)] <- NA

            }

        } else {

            x <- x

        }

        return(x)

    }

    sim2 <-
        sim %>%
        group_by(id) %>%
        mutate_at(
            .vars = vars( - (1:3) ),
            .funs = funs(bound_hit) ) %>%
        mutate_at(
            .vars = vars( - (1:3) ),
            .funs = funs(bound_na) ) %>%
        ungroup()

    sim2 <-
        sim %>%
        group_by(id) %>%
        mutate_("ER" = "bound_na(bound_hit(ER) )" ) %>%
        ungroup()

    nmax.hit <-
        sim2 %>%
        group_by(id) %>%
        gather_("index", "value", names(sim2[, - c(1:4)]) ) %>%
        mutate_("log_value" = "log(value)" ) %>%
        filter_(.dots = list(~n == max(n) ) ) %>%
        mutate(hitCondition = "nmax") %>%
        na.omit()

    boundary.hit <-
        sim2 %>%
        group_by(id) %>%
        gather_("index", "value", names(sim2[, - c(1:4)]) ) %>%
        mutate_("log_value" = "log(value)" ) %>%
        filter_(.dots = list(~log_value %in% logboundary) ) %>%
        mutate(hitCondition = "boundary") %>%
        na.omit()

    endpoint <- bind_rows(nmax.hit, boundary.hit)

    res <-
        endpoint %>%
        group_by_("index") %>%
        summarise_(
            "ASN" = "mean(n)",
            "ASN_sd" = "sd(n)",
            "Lower_hit" = "sum(value == 1 / boundary) / n_distinct(sim$id)",
            "Upper_hit" = "sum(value == boundary) / n_distinct(sim$id)",
            "Inconclusive" = "1 - (Lower_hit + Upper_hit)" ) %>%
        data.frame %>%
        mutate_at(.vars = 4:6, .funs = funs(paste0(. * 100, "%") ) )

    return(res)

}

#' @export

analysER.compER <- function(sim) {

    if (!any(class(sim) %in% c("simER", "compER") ) ) {

        stop("sim should be a simER or a compER object")

    }

    boundary <- as.numeric(unique(sim$boundary) )
    logboundary <- log(sort(c(boundary, 1 / boundary) ) )

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

            if (first < length(x) ) {

                x[(first + 1):length(x)] <- NA

            }

        } else if (any(x == 1 / boundary) ) {

            first <- which(x == 1 / boundary)[1]

            if (first < length(x) ) {

                x[(first + 1):length(x)] <- NA

            }

        } else {

            x <- x

        }

        return(x)

    }

    sim2 <-
        sim %>%
        group_by(id) %>%
        mutate_(
            "WAIC_ER" = "bound_na(bound_hit(WAIC_ER))",
            "LOO_ER" = "bound_na(bound_hit(LOO_ER))",
            "BF_SD" = "bound_na(bound_hit(BF_SD))",
            "BF_BS" = "bound_na(bound_hit(BF_BS))" ) %>%
        ungroup()

    nmax.hit <-
        sim2 %>%
        group_by(id) %>%
        gather_("index", "value", names(sim2)[6:9]) %>%
        mutate_("log_value" = "log(value)" ) %>%
        filter_(.dots = list(~n == max(n) ) ) %>%
        mutate(hitCondition = "nmax") %>%
        na.omit()

    boundary.hit <-
        sim2 %>%
        group_by(id) %>%
        gather_("index", "value", names(sim2)[6:9]) %>%
        mutate_("log_value" = "log(value)" ) %>%
        filter_(.dots = list(~log_value %in% logboundary) ) %>%
        mutate(hitCondition = "boundary") %>%
        na.omit()

    endpoint <- bind_rows(nmax.hit, boundary.hit)

    res <-
        endpoint %>%
        group_by_("index") %>%
        summarise_(
            "ASN" = "mean(n)",
            "ASN_sd" = "sd(n)",
            "Lower_hit" = "sum(value == 1 / boundary) / n_distinct(sim$id)",
            "Upper_hit" = "sum(value == boundary) / n_distinct(sim$id)",
            "Inconclusive" = "1 - (Lower_hit + Upper_hit)" ) %>%
        data.frame %>%
        mutate_at(.vars = 4:6, .funs = funs(paste0(. * 100, "%") ) )

    return(res)

}
