#' @title Runs the Monte Carlo method on a baseline no-synergy/antagonism
#'        mixture DER.
#'
#' @description
#'
#' @param n Numeric integer of the number of samples to be drawn.
#' @param dose Numeric vector of all total dose values to be evaluated.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of dose proportions applied on component DERs.
#' @param E Vector of DER functions
#' @param models List of nls models.
#' @param vcov Boolean for assessing inter-parameter correlation.
#' @param interval_length Numeric double of the confidence interval width.
#' @param seed Numeric value for pseudorandom generators.
#' @param ... Optional arguments to DER functions in E. Should be constant
#'            across all supplied E.

#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'
#' @return Named list representing lower and upper bounds for a Monte Carlo
#'         confidence interval for a mixture DER over an interval of doses.
#'
#' @examples
#'
#' ion_data <- load_ion_data("one_ion.csv")
#' HZE_nte_der <- make_der(ion_data, TRUE, TRUE)
#' check_der(HZE_nte_der)
#' # Iron and silicon two-ion mixture.
#' mc(dose = 0:100, LET = c(193, 70), ratios = c(1/3, 2/3),
#'    E = rep(c(HZE_nte_der), 2),
#'    models = rep(list(HZE_nte_model), 2), n = 20, vcov = TRUE)
#'
#' @author Yimin Lin, Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

mc <- function(dose, LET, ratios, E, models, dE = NULL, n = 200, vcov = TRUE,
               interval_length = 0.95, seed = 100, check = TRUE, ...) {
  if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  } else if (check) {
    for (DER in E) {
      if (!check_der(DER, ...)) {
        stop("DER has invalid properties.")
      }
    }
  }
  if (is.numeric(seed)) {
    set.seed(seed) # Set the pseudorandom seed
  }
  # Generate N randomly generated samples of parameters of HZE model.
  curve_list <- .generate_samples(n, dose, LET, ratios, E, dE,
                                  models, vcov, ...)
  ci <- matrix(nrow = length(dose), ncol = 2)

  # Calculate CI for each dose point
  for (i in 1:length(dose)) { # EGH: Possible vectorization opportunity
    ci[i, ] <- .generate_ci(length(curve_list), i, curve_list, interval_length)
  }
  colnames(ci) <- c("bottom", "top")
  return(data.frame(dose, ci))
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
#' @param dE
#' @param models
#' @param vcov Boolean for assessing inter-parameter correlation.
#' @param ... Optional arguments to DER functions in E.
#'
#' @return Numeric vector of sample mixture baseline DERs evaluated at the
#'         given doses.

.generate_samples <- function(n, dose, LET, ratios, E, dE, models, vcov, ...) {
  params <- curve_list <- list()
  if (vcov) {
    for (i in 1:length(models)) {
      params[[i]] <- mvtnorm::rmvnorm(n, mean  = stats::coef(models[[i]]),
                                         sigma = stats::vcov(models[[i]]))
    }
  } else {
    for (i in 1:length(models)) {
      params[[i]] <- sapply(rep(n, length(stats::coef(models[[i]]))),
                                          stats::rnorm,
                       mean = stats::coef(models[[i]]),
                         sd = summary(models[[i]])$coefficients[, "Std. Error"])
    }
  }
  start <- proc.time() # Starting time
  fail_count <- 0
  for (i in 1:n) {
    result <- tryCatch( # Counter implementation
      { # Prints step information
        out <- iea(dose, LET, ratios, E, dE,
                   coeff = sapply(params, function(x) x[i, ],
                   simplify = FALSE), check = FALSE, ...)
        curve_list[[length(curve_list) + 1]] <- out # Appends successful results
        message("Currently at Monte Carlo step: ", i, " of ", n,
                ". Elapsed time: ", (proc.time() - start)[["elapsed"]],
                " seconds.", rep(" ", 15), "\r", appendLF = FALSE)
        0 # Increments failure count by zero
      },
      error = function(cond) {
        message("Failure [", fail_count + 1, "] at step [", i,
                "] after [", (proc.time() - start)[["elapsed"]],
                "] seconds.", rep(" ", 50))
        message("\n", cond) # Prints original error message
        return(1) # Increments by one
      },
      warning = function(cond) {
        message("Warning at step ", i, " of ", n, " after ",
                (proc.time() - start)[["elapsed"]], " seconds.", rep(" ", 50))
        message("\n", cond) # Prints original warning message
        return(0) # Increments by zero
      },
      finally = { }
    )
    fail_count <- fail_count + result # Applies the increment
  }
  message("\n", "Total number of failures: ", fail_count)
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
