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

  # Error handling
  if (!is.data.frame(data)) {
    stop("Argument given for phi is not a data.frame.")
  } else if (is.na(as.numeric(as.character(phi)))) {
    stop("Argument given for phi is not numeric.")
  } else if (is.na(as.numeric(as.character(y_0)))) {
    stop("Argument given for y_0 is not numeric.")
  } else if (class(HZE) != "logical") {
    stop("Argument given for HZE is not boolean.")
  } else if (class(NTE) != "logical") {
    stop("Argument given for NTE is not boolean.")
  }

  # Swift light ================================================================
  if (!HZE) {
    low_LET_data <- dplyr::select(filter(data, Z < 4), 1:length(data))  # Swift light ions: here protons and alpha particles.

    low_LET_model <- nls(
      Prev ~ y_0 + 1 - exp( - alpha_low * dose), # alpha is used throughout radioiology for dose coefficients.
      data = low_LET_data,
      weights = NWeight, # Must have weight column called NWeight
      start = list(alpha_low = .005))

    low_LET_model_coef <- coef(low_LET_model) # Calibrated central values of the parameter.

    # Calibrated Low LET model. Use L = 0, but maybe later will use small L > 0.
    calibrated_low_LET_der <- function(dose, LET, alph_low = low_LET_model_coef[1]) {
      return(1 - exp( - alph_low * dose))
    }

    # .low_LET_slope <- function(dose, LET) { # Slope dE/dd of the low LET, low Z model.
    #   low_LET_model_coef * exp( - low_LET_model_coef * dose)
    # }
    model <- low_LET_model
    DER <- calibrated_low_LET_der
  }

  else { # HZE models ==========================================================

    HZE_data <- dplyr::select(filter(data, Z > 3), 1:length(data)) # Includes 1-ion data iff Z > 3

    # HZE NTE ------------------------------------------------------------------
    if (!NTE) {

      HZE_nte_model <- nls(  # Calibrate params in model modifying 17Cuc. hazard function NTE models. RKS: Models include Y_0, DERs do not.
        Prev ~ y_0 + (1 - exp ( - (aa1 * LET * dose * exp( - aa2 * LET)
                                   + (1 - exp( - phi * dose)) * kk1))),
        data = HZE_data,
        weights = NWeight, # Must have weight column called NWeight
        start = list(aa1 = .00009, aa2 = .001, kk1 = .06))

      # If a paper uses dose in Gy care is needed in preceding and following lines to rescale from cGy.

      HZE_nte_model_coef <- coef(HZE_nte_model) # Calibrated central values of the 3 parameters.
      aa1 <- HZE_nte_model_coef['aa1']
      aa2 <- HZE_nte_model_coef['aa2']
      kk1 <- HZE_nte_model_coef['kk1']

      # The DER, = 0 at dose 0.
      calibrated_nte_hazard_func <- function(dose, LET, coef) { # Calibrated hazard function.
        return(coef[1] * LET * dose * exp( - coef[2] * LET)
               + (1 - exp( - phi * dose)) * coef[3])
      }

      calibrated_HZE_nte_der <- function(dose, LET, coef = HZE_nte_model_coef) { # Calibrated HZE NTE DER.
        return(1 - exp( - calibrated_nte_hazard_func(dose, LET, coef)))
      }

      model <- HZE_nte_model
      DER <- calibrated_HZE_nte_der
    }


    # HZE TE -------------------------------------------------------------------
    if (!NTE) {

      HZE_te_model <- nls( # Calibrating parameters in a TE only model.
        Prev ~ y_0 + (1 - exp ( - (aate1 * LET * dose * exp( - aate2 * LET)))),
        data = HZE_data,
        weights = NWeight, # Must have weight column called NWeight
        start = list(aate1 = .00009, aate2 = .01))

      HZE_te_model_coef <- coef(HZE_te_model) # Calibrated central values of the 2 parameters.

      # The DER, = 0 at dose 0.
      calibrated_te_hazard_func <- function(dose, LET, coef) { # Calibrated hazard function.
        return(coef[1] * LET * dose * exp( - coef[2] * LET))
      }

      calibrated_HZE_te_der <- function(dose, LET, coef = HZE_te_model_coef) {
        return(1 - exp( - calibrated_te_hazard_func(dose, LET, coef))) # Calibrated HZE TE one-ion DER.
      }

      model <- HZE_te_model
      DER <- calibrated_HZE_te_der
    }

  }
  print(summary(model, correlation = TRUE))
  return(DER)
}

# Hidden functions =============================================================

# Swift light


# HZE NTE


# HZE TE

# Developer functions ----------------------------------------------------------
.test_runtime <- function(f, ...) { # Naive runtime check
  start_time <- Sys.time()
  f(...)
  Sys.time() - start_time
}
