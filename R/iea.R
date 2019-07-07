#' @title Applies IEA to get a baseline no-synergy/antagonism mixture DER.
#'
#' @description
#'
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param E Vector of dose effect relationship functions.
#' @param dE Optional vector of functions corresponding to the derivatives of E.
#' @param showSummary Boolean
#'
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'
#' @return Numeric vector representing the estimated Harderian Gland tumor
#'         prevalence from an IEA mixture DER constructed from the given
#'         one-ion DERs parameters.
#'
#' @examples
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

calculate_IEA <- function(dose, LET, ratios, model = "NTE",
                         coef = list(NTE = HZE_nte_model_coef,
                                     TE = HZE_te_model_coef,
                                     lowLET = low_LET_model_coef),
                         ders = list(NTE = calibrated_HZE_nte_der,
                                     TE = calibrated_HZE_te_der,
                                     lowLET = calibrated_low_LET_der),
                         calculate_dI = c(NTE = .calculate_dI_nte,
                                          TE = .calculate_dI_te),
                         phi = 2000) {
  dE <- function(yini, state, pars) { # Constructing an ODE from the DERS.
    with(as.list(c(state, pars)), {

      # Screen out low LET values.
      lowLET_total <- lowLET_ratio <- 0
      remove <- c()
      for (i in 1:length(LET)) {
        if (LET[i] <= 3) { # Ion is low-LET.
          lowLET_total <- lowLET_total + LET[i]
          lowLET_ratio <- lowLET_ratio + ratios[i]
          remove <- unique(c(remove, LET[i]))
          ratios[i] <- 0
        }
      }
      LET <- c(setdiff(LET, remove))
      ratios <- ratios[! ratios == 0]

      # Begin calculating dI values.
      aa <- u <- dI <- vector(length = length(LET))
      if (length(LET) > 0) {
        for (i in 1:length(LET)) {
          aa[i] <- pars[1] * LET[i] * exp( - pars[2] * LET[i])
          u[i] <- uniroot(function(dose) HZE_der(dose, LET[i], pars) - I,
                          interval = c(0, 20000),
                          extendInt = "yes",
                          tol = 10 ^ - 10)$root
          dI[i] <- ratios[i] * calc_dI(aa[i], u[i], pars[3])
        }
      }
      if (lowLET_ratio > 0) {
        # If low-LET DER is present then include it at the end of the dI vector. RKS: make it first as asked for above? last as asked for here? any location at all, which seems to work?
        u[length(LET) + 1] <- uniroot(function(dose)
          low_der(dose, LET = lowLET_total,
                  alph_low = coef[["lowLET"]]) - I,
          interval = c(0, 20000),
          extendInt = "yes",
          tol = 10 ^ - 10)$root
        dI[length(LET) + 1] <- lowLET_ratio * low_LET_slope(u[length(LET) + 1],
                                                            LET = lowLET_total)
      }
      return(list(sum(dI)))
    })
  }
  p <- list(pars = coef[[model]],
            HZE_der = ders[[model]],
            low_der = ders[["lowLET"]],
            calc_dI = calculate_dI[[model]])
  return(ode(c(I = 0), times = dose, dE, parms = p, method = "radau"))
}

#===============================================================================
iea <- function(dose, LET, ratios, E, dE = NULL, check = FALSE) {
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
    for (f in E) {
      if (!check_DER(f)) {
        stop("DER has invalid properties.")
      }
    }
  } # End error handling.

  dI <- function(t, y, parms) { # I'(d) = \sum_{j = 1}^{N} r_j * E'_j(D_j(I))
    with(as.list(c(y, parms)), { # Required syntax for ode call
      i <- 0
      for (j in 1:length(LET)) { # Loops over the summation
        dj <- function(y) E[[j]](y, LET[j]) # Stores the LET information
        if (!is.null(dE)) { # Finds numeric derivative
          slope <- function(y) dE[[j]](y, LET[j]) # Must require only one argument
        } else {
          slope <- function(y) numDeriv::grad(dj, y) # Slope function is accurate, inflated in past scripts
        }
        invrs <- .inverse(dj) # Finds numerical inverse
        i <- i + ratios[j] * slope(invrs(I)) # Increments the summation
        print(invrs(I))
      }
      return(list(i)) # I'(d)
    })
  }

  pars  <- c() # Setting up arguments to ode call
  yini  <- c(I = 0)
  times <- dose
  out   <- ode(yini, times, dI, pars) # method = "adams" recommended for nonstiff
  return(out)
}


# HIDDEN FUNCTIONS -------------------------------------------------------------
.calculate_dI_nte <- function(aa = 0.0000830213, u, kk1 = 0.0314282300) {
  return((aa + exp( - phi * u) * kk1 * phi) *
           exp( - (aa * u + (1 -exp( - phi * u)) * kk1)))
}

.calculate_dI_te <- function(aa, u, pars = NULL) {
  return(aa * exp(- aa * u))
}

.inverse <- function(f, lower = 10 ^ (-2), upper = 10 ^ 3) {
  function(y) uniroot(function(x) f(x) - y, lower = lower, upper = upper,
                      extendInt = "yes", maxiter = 10 ^ 4)$root
}
