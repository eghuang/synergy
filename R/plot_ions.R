#' @title Plots a mixed field.
#'
#' @description
#'
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param E Vector of DER functions taking dose and LET as arguments.
#' @param models
#' @param dE
#' @param n
#' @param ... Optional arguments to DER functions in E.
#'
#' @details
#'
#' @return
#'
#' @examples
#' ion_data      <- load_ion_data("one_ion.csv")
#' HZE_nte_der   <- make_der(ion_data, TRUE, TRUE)
#' HZE_nte_model <- make_model(ion_data, TRUE, TRUE)
#' plot_ions(0:100, c(193, 70), c(1/3, 2/3), E = c(HZE_nte_der, HZE_nte_der),
#'           models = list(HZE_nte_model, HZE_nte_model))
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

plot_ions <- function(dose, LET, ratios, E, models,
                      dE = NULL, n = 20, ...) {
  graphics.off() # Wipe plots

  # Contruct points for wide and narrow ribbons
  wci <- mc(dose, LET, ratios, E, models, dE, n = n, vcov = FALSE, ...)
  nci <- mc(dose, LET, ratios, E, models, dE, n = n, vcov = TRUE,
            check = FALSE, ...)

  # Plot window
  plot(c(dose[1], dose[length(dose)]), c(0, 1),
       col = "white", bty = 'L', ann = FALSE)

  # Plot ribbons
  polygon(x = c(wci$dose, rev(wci$dose)), y = c(wci$bottom, rev(wci$top)),
          xpd = -1, col = "yellow", lwd = .4, border = "orange") # Wide
  polygon(x = c(nci$dose, rev(nci$dose)), y = c(nci$bottom, rev(nci$top)),
          xpd = -1, col = "orange", lwd = .4, border = "orange") # Narrow

  # Plot DERs
  for (i in 1:length(LET)) {
    lines(dose, E[[i]](dose, LET[i], ...))
  }

  # Plot I(d)
  lines(dose, iea(dose, LET, ratios, E, dE, ...)$effect, col = "red")
}

# ion_data <- load_ion_data("one_ion.csv")
#
# LLET_der    <- make_der(ion_data, FALSE)
# HZE_nte_der <- make_der(ion_data, TRUE, TRUE)
# HZE_te_der  <- make_der(ion_data, TRUE, FALSE)
#
# LLET_model    <- make_model(ion_data, FALSE)
# HZE_nte_model <- make_model(ion_data, TRUE, TRUE)
# HZE_te_model  <- make_model(ion_data, TRUE, FALSE)
#
# check_der(LLET_der)
# check_der(HZE_nte_der)
# check_der(HZE_te_der)
#
# ci <- mc(dose = 0:100, LET = c(193, 70), ratios = c(1/3, 2/3),
#          E = rep(c(HZE_nte_der), 2),
#          models = rep(list(HZE_nte_model), 2), n = 20, vcov = TRUE)
#
# plot(c(0, 100), c(0, 0.6), col = "white", bty = 'L', ann = FALSE) # Set plot area
# polygon(x = c(ci$dose, rev(ci$dose)), y = c(ci$bottom, rev(ci$top)),
#         xpd = -1, col = "yellow", lwd = .4, border = "orange") # Narrow CI ribbon
# lines(0:100, HZE_nte_der(0:100, 193))
# lines(0:100, HZE_nte_der(0:100, 70))
# lines(0:100, iea(0:100, c(193, 70), c(1/3, 2/3), rep(c(HZE_nte_der), 2))$effect)
#
# plot_ions(0:100, c(193, 70), c(1/3, 2/3),
#           E      = c(HZE_nte_der, HZE_nte_der),
#           models = list(HZE_nte_model, HZE_nte_model))


