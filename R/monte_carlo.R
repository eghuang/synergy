#' @title Runs the Monte Carlo method on a baseline no-synergy/antagonism
#'        mixture DER.
#'
#' @description
#'
#' @param n Numeric integer of the number of samples to be drawn.
#' @param dose Numeric vector of all total dose values to be evaluated.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of dose proportions applied on component DERs.
#' @param E
#' @param models
#' @param vcov Boolean for assessing inter-parameter correlation.
#' @param interval_length Numeric double of the confidence interval width.
#' @param seed Numeric value for pseudorandom generators.
#' @param ... Optional arguments to DER functions in E.

#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'
#' @return Named list representing lower and upper bounds for a Monte Carlo
#'         confidence interval for a mixture DER over an interval of doses.
#'
#' @examples
#'
#' @author Yimin Lin, Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

monte_carlo <- function(dose, LET, ratios, E, models, n = 200, vcov = TRUE,
                        interval_length = 0.95, seed = 100, ...) {
  if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  }
  # Set the pseudorandom seed
  set.seed(seed)
  # Generate N randomly generated samples of parameters of HZE model.
  curve_list <- .generate_samples(n, dose, LET, ratios, E, models, vcov, ...)
  monte_carlo_ci <- matrix(nrow = 2, ncol = length(dose))

  # Calculate CI for each dose point
  for (i in 1:length(dose)) { # EGH: Possible vectorization opportunity
    monte_carlo_ci[, i] <- .generate_ci(n, i, curve_list, interval_length)
  }
  return(list(monte_carlo = monte_carlo_ci))
}


# Monte Carlo hidden functions =================================================

# Sampling ---------------------------------------------------------------------

#' @title Generates mixture no-synergy/antagonism DER samples with
#'        parameters drawn from a Gaussian distribution.
#'
#' @description
#'
#' @param n Numeric integer of the number of samples to be drawn.
#' @param dose Numeric vector of all total dose values to be evaluated.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of dose proportions applied on component DERs.
#' @param E
#' @param models
#' @param vcov Boolean for assessing inter-parameter correlation.
#' @param ... Optional arguments to DER functions in E.
#'
#' @return Numeric vector of sample mixture baseline DERs evaluated at the
#'         given doses.

.generate_samples <- function(n, dose, LET, ratios, E, models, vcov, ...) {
  params <- list()
  if (vcov) {
    for (i in 1:length(models)) {
      params[i] <- rmvnorm(n, mean = coef(models[i]), sigma = vcov(models[i]))
    }
  } else {
    for (i in 1:length(models)) {
      params[i] <- mapply(stats::rnorm, rep(n, length(coef(models[i]))),
                          coef(models[i]),
                          summary(models[i])$coefficients[, "Std. Error"])
    }
  }
  curve_list <- list()
  for (i in 1:n) {
    curve_list[[i]] <- iea(dose, LET, ratios, E, coef = params[[, ]][i, ], ...)
    cat(paste("  Currently at Monte Carlo step:", toString(i), "of",
              toString(n)), sprintf('\r'))
  }
  return(curve_list)
}


# Interval construction --------------------------------------------------------

#' @title Generates confidence intervals for mixture baseline DER samples
#'              at a dose.
#'
#' @description
#'
#' @param n Numeric integer of the number of samples to be drawn.
#' @param index Numeric integer of doses
#' @param sample_curves Numeric list of sampled mixture DER values.
#' @param interval_length Numeric double of the confidence interval width.
#'
#' @return Numeric length-two vector of an upper and lower bound for the
#'         confidence interval of a dose.

.generate_ci <- function(n, index, sample_curves, interval_length) {
  # For each sample curve, evalute them at input dose, and sort.
  sample_values <- sort(sapply(sample_curves, function(x) x[, 2][index]))
  # Returning resulting CI
  return(c(sample_values[ceiling((1 - interval_length) / 2 * n)],
           sample_values[(interval_length + (1 - interval_length) / 2) * n]))
}
