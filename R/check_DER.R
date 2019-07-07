#' @title Checks whether a DER satisfies the desired properties
#'
#' @description Checks whether the input DER is nonnegative when dose is
#'              nonnegative.
#'
#' @param DER Dose effect relationship function taking dose and LET as
#'            arguments.
#'
#' @details Only tests LET values in [1, 500].
#'
#' @return Boolean.
#'
#' @examples
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

check_DER <- function(DER, upper = 10 ^ 3, upperLET = 500) {
  if (any(DER(0, 1:upperLET) < 0)) { # Nonnegative when dose >= 0.
    return(FALSE)
  }
  for (i in 1:upperLET) {
    values <- DER(0:upper, i) # Check up to upper
    if (any(values > 1) || any(values < 0)) { # Prevalence does not exceed 1.
      return(FALSE)
    }
  }
  # Numerical inverse exists.
  message("DER verified up to ", upper, " cGy and ", upperLET, " LET.")
  return(TRUE) # All tests pass.
}
