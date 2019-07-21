#' @title Applies IEA to get a baseline no-synergy/antagonism mixture DER.
#'
#' @description
#'
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param E Vector of dose effect relationship functions.
#' @param dE Optional vector of functions corresponding to the derivatives of E.
#' @param coeff
#' @param check Optional boolean whether to check DERs in each run.
#' @param ... Optional arguments to DER functions in E.
#'
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'
#' @return Numeric vector representing the estimated Harderian Gland tumor
#'         prevalence from an IEA mixture DER constructed from the given
#'         one-ion DERs parameters.
#'
#' @examples
#' ion_data <- load_ion_data("one_ion.csv")
#' HZE_nte_der <- make_der(ion_data, TRUE, TRUE)
#' check_der(HZE_nte_der)
#' # Iron and silicon two-ion mixture.
#' iea(0:100, c(70, 30), c(1/3, 2/3), rep(c(HZE_nte_der), 2))
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

iea <- function(dose, LET, ratios, E,
                dE = NULL, coeff = NULL, check = FALSE, ...) {
  ll <- length(LET)
  lr <- length(ratios)
  le <- length(E)
 if (ll != 1 && lr != 1 && ll != lr) {
    stop("Length of LET and ratio arguments do not match.")
  } else if (le != 1 && le != min(ll, lr)) {
    return("Length of DERs and components do not match.")
  } else if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  } else if (!is.numeric(dose)) {
    stop("Dose argument is not numeric.")
  } else if (any(dose < 0)) {
    stop("Dose argument has negative elements.")
  } else if (!is.null(dE)) { # Derivative argument specified
    if (le != length(dE)) {
      stop("Length of DERs and derivatives do not match.")
    }
  }
  if (check) {
    for (DER in E) {
      if (!check_der(DER, ...)) {
        stop("DER has invalid properties.")
      }
    }
  } # End error handling.

  # Set up
  D <- list()
  for (m in 1:length(E)) {
    if (!is.null(coeff)) { # Should only be used in Monte Carlo calls
      D[[m]] <- function(y) E[[m]](y, LET[m], coeff = coeff[[m]], ...)
    } else { # Stores the LET information
      D[[m]] <- function(y) E[[m]](y, LET[m], ...)
    }

  }
  if (is.null(dE)) { # Finds numeric derivative
    dE <- list()
    for (k in 1:length(E)) {
      dE[[k]] <- function(y) numDeriv::grad(D[[k]], y)
    }
  }
  H <- list()
  for (n in 1:length(E)) {
    H[[n]] <- .inverse(D[[n]]) # Finds numerical inverse
  }
  # Add capability to cycle single DER or LET value

  dI <- function(t, y, parms) { # I'(d) = \sum_{j = 1}^{N} r_j * E'_j(D_j(I))
    with(as.list(c(y, parms)), { # Required syntax for ode call
      i <- 0
      for (j in 1:length(LET)) { # Loops over the summation
        slope <- function(y) dE[[j]](y) # Must require only one argument
        i <- i + ratios[j] * slope(H[[j]](I)) # Increments the summation
      }
      return(list(i)) # I'(d)
    })
  }

  pars  <- c() # Setting up arguments to ode call
  yini  <- c(I = 0)
  times <- dose
  out   <- data.frame(deSolve::ode(yini, times, dI, pars, method = "bdf_d")) # "adams" method recommended for nonstiff
  colnames(out) <- c("dose", "effect")
  return(out)
}


# Hidden functions =============================================================

.inverse <- function(f, lower = 10 ^ (-2), upper = 10 ^ 4) {
  function(y) stats::uniroot(function(x) f(x) - y, lower = lower, upper = upper,
                             extendInt = "yes", maxiter = 10 ^ 4)$root
}
