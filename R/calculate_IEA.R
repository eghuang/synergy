#' @title Applies IEA to get a baseline no-synergy/antagonism mixture DER.
#'
#' @description
#'
#' @param dose Numeric vector corresponding to the sum dose in cGy.
#' @param LET Numeric vector of all LET values, must be length n.
#' @param ratios Numeric vector of all dose ratios, must be length n.
#' @param model String value corresponding to the model to be used, either "NTE" or "TE".
#' @param coef Named list of numeric vectors containing coefficients for one-ion DERs.
#' @param ders Named list of functions containing relevant one-ion DER models.
#' @param calculate_dI Named vector of functions to calculate dI depending on
#'                     the selected model.
#' @param phi Numeric value, used in NTE models.
#'
#' @details Corresponding elements of ratios, LET should be associated with the
#'          same DER.
#'
#' @return Numeric vector representing the estimated Harderian Gland
#'         prevalence from an IEA mixture DER constructed from the given
#'         one-ion DERs parameters.
#'
#' @examples
#' calculate_IEA(.01 * 0:40, c(70, 195), c(1/2, 1/2))
#' calculate_IEA(.01 * 0:70, c(.4, 195), c(4/7, 3/7), model = "TE")
#' @export

calculate_IEA <- function(dose, LET, ratios, model = "NTE",
                         coef = list(NTE = HZE_nte_model_coef,
                                     TE = HZE_te_model_coef,
                                     lowLET = low_LET_model_coef),
                         ders = list(NTE = calibrated_HZE_nte_der,
                                     TE = calibrated_HZE_te_der,
                                     lowLET = calibrated_low_LET_der),
                         calculate_dI = c(NTE = .calculate_dI_nte,
                                          TE = .calculate_dI_te),
                         phi = 2000) {
  dE <- function(yini, state, pars) { # Constructing an ODE from the DERS.
    with(as.list(c(state, pars)), {

      # Screen out low LET values.
      lowLET_total <- lowLET_ratio <- 0
      remove <- c()
      for (i in 1:length(LET)) {
        if (LET[i] <= 3) { # Ion is low-LET.
          lowLET_total <- lowLET_total + LET[i]
          lowLET_ratio <- lowLET_ratio + ratios[i]
          remove <- unique(c(remove, LET[i]))
          ratios[i] <- 0
        }
      }
      LET <- c(setdiff(LET, remove))
      ratios <- ratios[! ratios == 0]

      # Begin calculating dI values.
      aa <- u <- dI <- vector(length = length(LET))
      if (length(LET) > 0) {
        for (i in 1:length(LET)) {
          aa[i] <- pars[1] * LET[i] * exp( - pars[2] * LET[i])
          u[i] <- uniroot(function(dose) HZE_der(dose, LET[i], pars) - I,
                          interval = c(0, 20000),
                          extendInt = "yes",
                          tol = 10 ^ - 10)$root
          dI[i] <- ratios[i] * calc_dI(aa[i], u[i], pars[3])
        }
      }
      if (lowLET_ratio > 0) {
        # If low-LET DER is present then include it at the end of the dI vector. RKS: make it first as asked for above? last as asked for here? any location at all, which seems to work?
        u[length(LET) + 1] <- uniroot(function(dose)
          low_der(dose, LET = lowLET_total,
                  alph_low = coef[["lowLET"]]) - I,
          interval = c(0, 20000),
          extendInt = "yes",
          tol = 10 ^ - 10)$root
        dI[length(LET) + 1] <- lowLET_ratio * low_LET_slope(u[length(LET) + 1],
                                                            LET = lowLET_total)
      }
      return(list(sum(dI)))
    })
  }
  p <- list(pars = coef[[model]],
            HZE_der = ders[[model]],
            low_der = ders[["lowLET"]],
            calc_dI = calculate_dI[[model]])
  return(ode(c(I = 0), times = dose, dE, parms = p, method = "radau"))
}

#============= dI HIDDEN FUNCTIONS =============#
.calculate_dI_nte <- function(aa, u, kk1) {
  return((aa + exp( - phi * u) * kk1 * phi) *
           exp( - (aa * u + (1 -exp( - phi * u)) * kk1)))
}

.calculate_dI_te <- function(aa, u, pars = NULL) {
  return(aa * exp(- aa * u))
}
