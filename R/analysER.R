#' Analyse the results of simulations ran with compER
#'
#' Analyse the results of simulations ran with compER. It computes the average
#' sample number (ASN) at which the boundary is attained (either the lower or
#' the upper one), the percentage of hits of the lower boundary as well as hits
#' of the upper boundary, and the percentage of trajectories that did not hit
#' none of the boundaries.
#'
#' @param sim A compER object.
#'
#' @return An object of class \code{data.frame}, which contains the average
#' sample number (ASN) at which the boundary is attained (either the lower or
#' the upper one), the percentage of hits of the lower boundary as well as hits
#' of the upper boundary, and the percentage of trajectories that did not hit
#' none of the boundaries.
#'
#' @importFrom stats sd na.omit
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(ESTER)
#' sim <- compER(cohensd = 0.8, nmin = 20, nmax = 100, boundary = 10, B = 20, cores = 4)
#' analysER(sim)
#' }
#'
#' @author Ladislas Nalborczyk <\email{ladislas.nalborczyk@@gmail.com}>
#'
#' @seealso \code{\link{compER}}
#'
#' @export

analysER <- function(sim) {

    if (!any(class(sim) == "compER") ) {

        stop("sim should be a compER object")

    }

    boundary <- unique(sim$boundary)
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
        gather_("index", "value", names(sim2)[5:8]) %>%
        mutate_("log_value" = "log(value)" ) %>%
        #filter(n == max(n),
        #    max(log_value, na.rm = TRUE) <= logBoundary[2] &
        #        min(log_value, na.rm = TRUE) >= logBoundary[1]) %>%
        filter_(.dots = list(~n == max(n) ) ) %>%
            #~(max(log_value, na.rm = TRUE) <= logBoundary[2] &
            #    min(log_value, na.rm = TRUE) >= logBoundary[1]) ) ) %>%
        mutate(hitCondition = "nmax") %>%
        na.omit()

    boundary.hit <-
        sim2 %>%
        group_by(id) %>%
        gather_("index", "value", names(sim2)[5:8]) %>%
        mutate_("log_value" = "log(value)" ) %>%
        filter_(.dots = list(~log_value %in% logBoundary) ) %>%
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
        data.frame

    class(res) <- c("analysER", "data.frame")

    return(res)

}

#' @export

print.analysER <- function(x, digits = 2, ... ) {

    x %>%
        lapply(function(y) if (is.numeric(y) ) round(y, digits) else y) %>%
        data.frame %>%
        print

}
