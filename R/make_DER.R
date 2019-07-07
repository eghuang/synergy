#' @title Makes a dose effect relationship function for a specified dataset.
#'
#' @description Outputs a dose effect relationship (DER) for a single ion
#'              with parameters dose and LET. The ion may be a low-LET or HZE
#'              ion and non-targeted effects may be specified.
#'
#' @param data Data.frame of single ion data.
#' @param phi Numeric constant
#' @param y_0 Numeric constant corresponding to background effect
#' @param HZE Boolean for presence of HZE ions
#' @param NTE Boolean for presence of non-targeted effects
#'
#'
#' @details data must have a numeric column "dose" and a numeric weight column
#'          "NWeight". If a paper uses dose in Gy care is needed in preceding
#'          and following lines to rescale from cGy.
#'
#' @return Function
#'
#' @examples
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

make_DER <- function(data, HZE, NTE = TRUE, phi = 2000, y_0 = 0.046404) {
  # Error handling in make_model
  model <- make_model(data, HZE, NTE, phi, y_0)
  model_coef <- coef(model)

  # Swift light ================================================================
  if (!HZE) {
    # Calibrated Low LET model. Use L = 0, but maybe later will use small L > 0.
    DER <- function(dose, LET, alph_low = model_coef) {
      return(1 - exp( - alph_low * dose))
    }

  # HZE models =================================================================
  } else {

    if (NTE) { # HZE NTE -------------------------------------------------------
      # Calibrated hazard function.
      HZE_NTE_hazard_func <- function(dose, LET, coef) {
        return(coef[1] * LET * dose * exp( - coef[2] * LET)
               + (1 - exp( - phi * dose)) * coef[3])
      }
      # Calibrated HZE NTE DER.
      DER <- function(dose, LET, coef = model_coef) {
        return(1 - exp( - HZE_NTE_hazard_func(dose, LET, coef)))
      }

    } else if (!NTE) { # HZE TE ------------------------------------------------
      # Calibrated hazard function.
      HZE_TE_hazard_func <- function(dose, LET, coef) {
        return(coef[1] * LET * dose * exp( - coef[2] * LET))
      }
      # Calibrated HZE TE one-ion DER.
      DER <- function(dose, LET, coef = model_coef) {
        return(1 - exp( - HZE_TE_hazard_func(dose, LET, coef)))
      }
    }
  }
  print(summary(model, correlation = TRUE))
  return(DER)
}


# Hidden functions =============================================================

.low_LET_slope <- function(dose, LET) { # Slope dE/dd of the low LET, low Z model.
  low_LET_model_coef * exp( - low_LET_model_coef * dose)
}

.HZE_NTE_slope <- function(aa, u, kk1) {
  return((aa + exp( - phi * u) * kk1 * phi) *
           exp( - (aa * u + (1 -exp( - phi * u)) * kk1)))
}

.HZE_TE_slope <- function(aa, u, pars = NULL) {
  return(aa * exp(- aa * u))
}
