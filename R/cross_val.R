#' @title Applies leave-one-out cross validation to a mixture of ions with
#'        respect to the NTE or TE models used from Huang 2019.
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
#' @return Numeric weighted mean errors.
#'
#' @details Tested for the two examples below.
#'
#' @examples
#' ##=================== Cross validation ====================#
#' # Seperate Data into 8 blocks, i.e. test/training sets:
#' HZE_data <- select(filter(ion_data, Z > 3), 1:length(ion_data))
#' data_len <- 1:length(HZE_data)
#' O_350 <- select(filter(HZE_data, LET == 20), data_len)
#' Ne_670 <- select(filter(HZE_data, LET == 25), data_len)
#' Si_260 <- select(filter(HZE_data, LET == 70), data_len)
#' Ti_1000 <- select(filter(HZE_data, LET == 100), data_len)
#' Fe_600 <- select(filter(HZE_data, LET == 193), data_len)
#' Fe_350 <- select(filter(HZE_data, LET == 253), data_len)
#' Nb_600 <- select(filter(HZE_data, LET == 464), data_len)
#' La_593 <- select(filter(HZE_data, LET == 953), data_len)
#' set_list <- list(O_350, Ne_670, Si_260, Ti_1000,
#'                  Fe_600, Fe_350, Nb_600, La_593)
#' actual_prev <- HZE_data$Prev
#'
#' NTE_cv <- cross_val(set_list,"NTE", HZE_data$Prev, HZE_data$NWeight)
#' TE_cv <- cross_val(set_list, "TE", HZE_data$Prev, HZE_data$NWeight)
#' cv_table <- cbind(NTE_cv, TE_cv)
#' cv_table
#'
#' @author Edward Greg Huang <eghuang@@berkeley.edu>
#' @export

loo_cv <- function(ions, prev, wgts, model = c("NTE", "TE"),
                      y_0 = 0.046404, phi = 2000) { # Leave one out
  pred <- c()
  for (i in 1:length(ions)) {
    test <- ions[[i]]
    train <- Reduce(rbind, ions[-i]) # Constructs excluded dataframe
    if (model == "NTE") { # Fits NTE model
      f <- Prev ~ y_0 + (1 - exp ( -(aa1 * LET * dose * exp( -aa2 * LET)
                                     + (1 - exp( - phi * dose)) * kk1)))
      init <- list(aa1 = .00009, aa2 = .001, kk1 = .06)
    } else if (model == "TE") { # Fits TE model
      f <- Prev ~ y_0 + (1 - exp ( -(aate1 * LET * dose * exp( -aate2 * LET))))
      init <- list(aate1 = .00009, aate2 = .01)
    } else { # Error handling
      stop("Invalid model argument specified.")
    }
    modl <- nls(f, data = train, weights = NWeight, start = init)
    pred <- c(pred, predict(modl, test)) # Appends predictions
  }
  errors <- (pred - prev) ^ 2
  return(weighted.mean(errors, wgts)) # Weighted mean squared error
}
