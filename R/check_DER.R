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

check_DER <- function(DER) {
  if (sum(DER(0, 1:500) < 0) < 0) { # Nonnegative when dose >= 0.
    return(FALSE)
  } else if (integrate(DER,
                       lower = 0, upper = .Machine$integer.max, LET = 1:500) !=
             integrate(function(x, LET) abs(DER(x, LET)),
                       lower = 0, upper = .Machine$integer.max, LET = 1:500)) {
    return(FALSE)
  }
  else if (DER(0:.Machine$integer.max, 1:500) > 1) { # Prevalence does not exceed 1.
    return(FALSE)
  }
  # Numerical inverse exists.
  return(TRUE) # All tests pass.
}
