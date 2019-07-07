#' @title Applies Simple Effect Additivity to get a baseline mixture DER.
#'
#' @description
#'
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param E Vector of DER functions taking dose and LET as arguments.
#' @param ... Optional arguments to DER functions in E.
#'
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'
#' @return Numeric vector representing the estimated Harderian Gland
#'         prevalence from a SEA mixture DER constructed from the given DER
#'         parameters.
#'
#' @examples
#' sea(.01 * 0:40, c(70, 195), c(1/2, 1/2), n = 2)
#' sea(.01 * 0:70, c(0.4, 195), c(4/7, 3/7))
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

sea <- function(dose, LET, ratios, E, ...) {
  if (length(LET) != 1 && length(ratios) != 1
      && length(LET) != length(ratios)) {
    stop("Length of LET and ratio arguments do not match.")
  } else if (length(E) != 1 &&
             length(E) != min(length(LET), length(ratios))) {
    return("Length of DERs and components do not match.")
  } else if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  }
  for (DER in E) {
    if (!check_DER(DER, ...)) {
      stop("DER has invalid properties.")
    }
  }
  # End error handling.
  total <- 0
  for (i in 1:length(ratios)) { # Iterate over HZE ions in the mixture.
    total <- total + E[[i]](dose * ratios[i], LET[i], ...)
  }
  return(total)
}
