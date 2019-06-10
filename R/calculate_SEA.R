#' @title Applies Simple Effect Additivity to get a baseline mixture DER.
#'
#' @description
#'
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param lowLET Boolean of whether a low LET DER should be included in the mixture DER.
#' @param n Number of DERs, optional argument used to check parameter validity.
#'
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'
#' @return Numeric vector representing the estimated Harderian Gland
#'         prevalence from a SEA mixture DER constructed from the given DER
#'         parameters.
#'
#' @examples
#' calculate_SEA(.01 * 0:40, c(70, 195), c(1/2, 1/2), n = 2)
#' calculate_SEA(.01 * 0:70, c(0.4, 195), c(4/7, 3/7))
#' # RKS please explain n=NULL and lowLET = false a bit more. Also, what happens if we add more models in the cross validation?
#' @export

calculate_SEA <- function(dose, LET, ratios, lowLET = FALSE, n = NULL) {
  if (!is.null(n) && (n != length(ratios) | n != length(LET))) {
    stop("Length of arguments do not match.")
  } else if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  } #  End error handling
  total <- 0
  i <- 1
  if (lowLET == TRUE) {
    # First elements of ratios and LET should be the low-LET DER.
    #RKS is this still true, or is low LET now identified by LET, not location in a vector?
    total <- total + calibrated_low_LET_der(dose * ratios[i], LET[i])
    i <- i + 1
  }
  while (i < length(ratios) + 1) { # Iterate over HZE ions in the mixture.
    total <- total + calibrated_HZE_nte_der(dose * ratios[i], LET[i])
    i <- i + 1
  }
  return(total)
}
