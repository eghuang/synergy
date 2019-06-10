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
  # Nonnegative when dose >= 0
  if (sum(DER(0, 1:500) < 0) < 0) {
    return(FALSE)
  } else if (integrate(DER,
                       lower = 0, upper = .Machine$integer.max, LET = 1:500) !=
             integrate(function(x, LET) abs(DER(x, LET)),
                       lower = 0, upper = .Machine$integer.max, LET = 1:500)) {
    return(FALSE)
  }

  # All tests pass
  return(TRUE)
}
