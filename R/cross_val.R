#' @title Applies cross validation to a mixture of ions with respect to
#'        a synergy theory.
#'
#' @description
#'
#' @param ions List of dataframes corresponding to ion.
#' @param model String, "NTE" or "TE" for non-targeted or targeted effects.
#' @param prev Numeric vector of observed prevalence values.
#' @param w Numeric vector of experimental weights.
#'
#' @details Weight vector elements should correspond to dataframe element order.
#'          i.e. w[n] = ions[length(ions[:, 1]) / n][length(ions[:, 1]) % n]
#'
#' @return Numeric vector representing the estimated Harderian Gland
#'         prevalence from a SEA mixture DER constructed from the given DER
#'         parameters. RKS 5/17/2019. Again: why SEA?
#'
#' @details Tested for the two examples below. General correctness is
#'          unverified and should be treated at some future data.
#'
#' @examples
#' ##=================== Cross validation ====================#
#' # Seperate Data into 8 blocks, i.e. test/training sets:
#' data_len <- 1:length(HZE_data)
#' O_350 <- dplyr::select(filter(HZE_data, Beam == "O"), data_len) #RKS added 3 oxygen points 5/17/2019)
#' Ne_670 <- dplyr::select(filter(HZE_data, Beam == "Ne"), data_len)
#' Si_260 <- dplyr::select(filter(HZE_data, Beam == "Si"), data_len)
#' Ti_1000 <- dplyr::select(filter(HZE_data, Beam == "Ti"), data_len)
#' Fe_600 <- dplyr::select(filter(HZE_data, LET == 193), data_len)
#' Fe_350 <- dplyr::select(filter(HZE_data, LET == 253), data_len)
#' Nb_600 <- dplyr::select(filter(HZE_data, Beam == "Nb"), data_len)
#' La_593 <- dplyr::select(filter(HZE_data, Beam == "La"), data_len)
#' set_list <- list(O_350, Ne_670, Si_260, Ti_1000,
#'                  Fe_600, Fe_350, Nb_600, La_593)
#' actual_prev <- HZE_data$Prev
#'
#' NTE_cv <- cross_val(set_list,"NTE", HZE_data$Prev, HZE_data$NWeight)
#' TE_cv <- cross_val(set_list, "TE", HZE_data$Prev, HZE_data$NWeight)
#'
#' @export

cross_val <- function(ions, model, prev, w) {
  theoretical <- vector()
  for (i in 1:length(ions)) {
    test <- ions[[i]]
    excluded_list <- ions[-i]
    train <- excluded_list[[1]]
    for (j in 2:length(excluded_list)) {
      train <- rbind(train, excluded_list[[j]])
    }
    if (model == "NTE") {
      f <- Prev ~ Y_0 + (1 - exp ( -(aa1 * LET * dose * exp( -aa2 * LET)
                                     + (1 - exp( - phi * dose)) * kk1)))
      s <- list(aa1 = .00009, aa2 = .001, kk1 = .06)
    } else {
      f <- Prev ~ Y_0 + (1 - exp ( -(aate1 * LET * dose * exp( -aate2 * LET))))
      s <- list(aate1 = .00009, aate2 = .01)
    }
    m <- nls(f, data = train, weights = NWeight, start = s)
    pred <- predict(m, test)
    theoretical <- c(theoretical, pred)
  }
  errors <- (theoretical - prev) ^ 2
  return(weighted.mean(errors, w))
}

