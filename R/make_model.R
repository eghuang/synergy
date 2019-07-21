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
#' ion_data <- load_ion_data("one_ion.csv")
#'
#' LLET_model    <- make_model(ion_data, FALSE)
#' HZE_nte_model <- make_model(ion_data, TRUE, TRUE)
#' HZE_te_model  <- make_model(ion_data, TRUE, FALSE)
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

make_model <- function(data, HZE, NTE = TRUE, phi = 2000, y_0 = 0.046404) {
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
  if (!HZE) { # Swift light ions: here protons and alpha particles.
    low_LET_data <- dplyr::select(dplyr::filter(data, Z < 4), 1:length(data))

    model <- stats::nls(
      Prev ~ y_0 + 1 - exp( - alpha_low * dose), # alpha is used throughout radioiology for dose coefficients.
      data    = low_LET_data,
      weights = NWeight, # Must have weight column called NWeight
      start   = list(alpha_low = .005))

  # HZE models =================================================================
  } else {
    HZE_data <- dplyr::select(dplyr::filter(data, Z > 3), 1:length(data)) # Includes 1-ion data iff Z > 3

    if (NTE) { # HZE NTE -------------------------------------------------------
      model <- stats::nls(  # Calibrate params in model modifying 17Cuc. hazard function NTE models.
        Prev ~ y_0 + (1 - exp ( - (aa1 * LET * dose * exp( - aa2 * LET)
                                   + (1 - exp( - phi * dose)) * kk1))),
        data    = HZE_data,
        weights = NWeight, # Must have weight column called NWeight
        start   = list(aa1 = .00009, aa2 = .001, kk1 = .06))

    } else if (!NTE) { # HZE TE ------------------------------------------------
      model <- stats::nls( # Calibrating parameters in a TE only model.
        Prev ~ y_0 + (1 - exp ( - (aate1 * LET * dose * exp( - aate2 * LET)))),
        data    = HZE_data,
        weights = NWeight, # Must have weight column called NWeight
        start   = list(aate1 = .00009, aate2 = .01))
    }
  }
  return(model)
}
