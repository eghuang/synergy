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
#' ion_data <- load_ion_data("one_ion.csv")
#' HZE_nte_der <- make_der(ion_data, TRUE, TRUE)
#' check_der(HZE_nte_der)
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

check_der <- function(DER, upper = 10 ^ 3, upperLET = 500, coeff = NULL, ...) {
  if (!is.null(coeff)) {
    if (any(DER(0, 1:upperLET, coeff, ...) < 0)) { # Nonnegative when dose >= 0.
      return(FALSE)
    }
    for (i in 1:upperLET) {
      values <- DER(0:upper, i, coeff, ...) # Check up to upper
      if (any(values > 1) || any(values < 0)) { # Prevalence does not exceed 1.
        return(FALSE)
      }
    }
    message(deparse(substitute(DER)), " is verified up to ", upper,
            " cGy and ", upperLET, " LET for the supplied parameters.")
  } else {
    if (any(DER(0, 1:upperLET, ...) < 0)) { # Nonnegative when dose >= 0.
      return(FALSE)
    }
    for (i in 1:upperLET) {
      values <- DER(0:upper, i, ...) # Check up to upper
      if (any(values > 1) || any(values < 0)) { # Prevalence does not exceed 1.
        return(FALSE)
      }
    }
    message(deparse(substitute(DER)), " is verified up to ", upper,
            " cGy and ", upperLET, " LET for the default parameters.")
  }
  # Numerical inverse exists.
  invisible(TRUE) # All tests pass.
}
