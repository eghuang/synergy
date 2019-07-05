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
#'          "NWeight".
#'
#' @return Function
#'
#' @examples
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

make_DER <- function(data, HZE, NTE, phi = 2000, y_0 = 0.04604) {

  # Error handling in make_model
  model <- make_model(data, HZE, NTE, phi = 2000, y_0 = 0.04604)
  model_coef <- coef(model)

  # Swift light ================================================================
  if (!HZE) {
    # Calibrated Low LET model. Use L = 0, but maybe later will use small L > 0.
    DER <- function(dose, LET, alph_low = model_coef) {
      return(1 - exp( - alph_low * dose))
    }
    # .low_LET_slope <- function(dose, LET) { # Slope dE/dd of the low LET, low Z model.
    #   low_LET_model_coef * exp( - low_LET_model_coef * dose)
    # }

  # HZE models =================================================================
  } else {

    if (NTE) { # HZE NTE -------------------------------------------------------
      # Calibrated hazard function.
      calibrated_nte_hazard_func <- function(dose, LET, coef) {
        return(coef[1] * LET * dose * exp( - coef[2] * LET)
               + (1 - exp( - phi * dose)) * coef[3])
      }
      # Calibrated HZE NTE DER.
      DER <- function(dose, LET, coef = model_coef) {
        return(1 - exp( - calibrated_nte_hazard_func(dose, LET, coef)))
      }

    } else if (!NTE) { # HZE TE ------------------------------------------------
      # Calibrated hazard function.
      calibrated_te_hazard_func <- function(dose, LET, coef) {
        return(coef[1] * LET * dose * exp( - coef[2] * LET))
      }
      # Calibrated HZE TE one-ion DER.
      DER <- function(dose, LET, coef = model_coef) {
        return(1 - exp( - calibrated_te_hazard_func(dose, LET, coef)))
      }

    }
  }
  print(summary(model, correlation = TRUE))
  return(DER)
}

# Model generation =============================================================

#' @title Makes a dose effect relationship model for a specified dataset.
#'
#' @description Outputs a dose effect relationship (DER) model for a single ion
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
#'          "NWeight".
#'
#' @return nls model
#'
#' @examples
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

make_model <- function(data, HZE, NTE, phi = 2000, y_0 = 0.04604) {
  # Error handling =============================================================
  if (!is.data.frame(data)) {
    stop("Argument given for data is not a data.frame.")
  } else if (class(HZE) != "logical") {
    stop("Argument given for HZE is not boolean.")
  } else if (class(NTE) != "logical") {
    stop("Argument given for NTE is not boolean.")
  } else if (is.na(as.numeric(as.character(phi)))) {
    stop("Argument given for phi is not numeric.")
  } else if (is.na(as.numeric(as.character(y_0)))) {
    stop("Argument given for y_0 is not numeric.")
  }
  # Swift light ================================================================
  if (!HZE) {
    # Swift light ions: here protons and alpha particles.
    low_LET_data <- dplyr::select(filter(data, Z < 4), 1:length(data))

    model <- nls(
      Prev ~ y_0 + 1 - exp( - alpha_low * dose), # alpha is used throughout radioiology for dose coefficients.
      data = low_LET_data,
      weights = NWeight, # Must have weight column called NWeight
      start = list(alpha_low = .005))

  # HZE models ========================================================-========
  } else {
    HZE_data <- dplyr::select(filter(data, Z > 3), 1:length(data)) # Includes 1-ion data iff Z > 3

    if (NTE) { # HZE NTE -------------------------------------------------------
      model <- nls(  # Calibrate params in model modifying 17Cuc. hazard function NTE models. RKS: Models include Y_0, DERs do not.
        Prev ~ y_0 + (1 - exp ( - (aa1 * LET * dose * exp( - aa2 * LET)
                                   + (1 - exp( - phi * dose)) * kk1))),
        data = HZE_data,
        weights = NWeight, # Must have weight column called NWeight
        start = list(aa1 = .00009, aa2 = .001, kk1 = .06))

    } else if (!NTE) { # HZE TE ------------------------------------------------
      model <- nls( # Calibrating parameters in a TE only model.
        Prev ~ y_0 + (1 - exp ( - (aate1 * LET * dose * exp( - aate2 * LET)))),
        data = HZE_data,
        weights = NWeight, # Must have weight column called NWeight
        start = list(aate1 = .00009, aate2 = .01))

      # If a paper uses dose in Gy care is needed in preceding and following lines to rescale from cGy.
    }
  }
  return(model)
}

# Hidden functions =============================================================

