#' @title Applies Simple Effect Additivity to get a baseline mixture DER.
#'
#' @description
#'
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param DER DER function taking dose and LET as arguments.
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
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

calculate_SEA <- function(dose, LET, ratios, DER, n = NULL) {
  if (!is.null(n) && (n != length(ratios) | n != length(LET))) {
    stop("Length of arguments do not match.")
  } else if (sum(ratios) != 1) {
    stop("Sum of ratios do not add up to one.")
  }
  # else if (!check_DER(DER)) {
  # stop("DER has invalid properties".)
  # }
  #  End error handling
  total <- 0
  i <- 0
  while (i < length(ratios)) { # Iterate over HZE ions in the mixture.
    total <- total + DER(dose * ratios[i], LET[i])
    i <- i + 1
  }
  return(total)
}
