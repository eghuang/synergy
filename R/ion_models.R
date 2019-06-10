# Copyright:    (C) 2017-2019 Sachs Undergraduate Research Apprentice Program
#               This program and its accompanying materials are distributed
#               under the terms of the GNU General Public License v3.
# Filename:     synergy.R
# Purpose:      Concerns radiogenic mouse Harderian gland tumorigenesis.
#               Contains relevant synergy theory models, information coefficient
#               calculations and useful objects. It is part of the
#               source code for the Chang 2019 HG project.
# Contact:      Rainer K. Sachs
# Website:      https://github.com/rainersachs/mouseHG_Chang_2019plus
# Mod history:  04 Apr 2019
# Details:      See data_info.R for further licensing, attribution,
#               references, and abbreviation information.

# source("data_info.R") # Load in the data. Remark: dose is in units of cGy;
# LET usually in keV/micron; prevalence Prev always < 1
# (i.e. not in %, which would mean prevalence < 100 but is strongly deprecated).

# library(deSolve) # Solving differential equations.
# library(minpack.lm) # Package for non-linear regression
# library(mvtnorm) # Package for calculating confidence intervals by Monte Carlo
# simulation based on variance-covariance matrices
# library(Hmisc) # Plotting, e.g. ribbons
# library(dplyr) # Helps manipulate data frames

#========================= MISC. OBJECTS & VARIABLES ==========================#
# In next line phi controls how fast NTE build up from zero; not really needed
# during calibration since phi * Dose >> 1 at every observed Dose !=0. phi is
# needed for later synergy calculations.
# d_0 = 1 / phi = 5 x 10-4 cGy = 5 x 10^-6 Gy.

Y_0 <- 0.04604 # HG tumor prevalence for sham irradiated controls.
# Y_0 SD = 0.01 but this value is not used.

phi <- 2000 # even larger phi should give the same final results,
# but might cause extra problems with R.

#================================ DER MODELS ==================================#

#=============== Subsets of 1-ion dataframe =================#
# (HZE = high charge and energy;
HZE_data <- dplyr::select(filter(ion_data, Z > 3), 1:length(ion_data)) # Includes 1-ion data iff Z > 3
low_LET_data <- dplyr::select(filter(ion_data, Z < 4), 1:length(ion_data))  # Swift light ions: here protons and alpha particles.

#=============== HZE/NTE MODEL =================#
#  NTE = non-targeted effects are included (in addition to TE)
HZE_nte_model <- nls(  # Calibrate params in model modifying 17Cuc. hazard function NTE models. RKS: Models include Y_0, DERs do not.
  Prev ~ Y_0 + (1 - exp ( - (aa1 * LET * dose * exp( - aa2 * LET)
                             + (1 - exp( - phi * dose)) * kk1))),
  data = HZE_data, weights = NWeight,
  start = list(aa1 = .00009, aa2 = .001, kk1 = .06))

summary(HZE_nte_model, correlation = TRUE) # Parameter values & accuracy.

# If a paper uses dose in Gy care is needed in preceding and following lines to rescale from cGy.
vcov(HZE_nte_model) # Variance-covariance matrix.
HZE_nte_model_coef <- coef(HZE_nte_model) # Calibrated central values of the 3 parameters.
aa1 <- HZE_nte_model_coef['aa1']
aa2 <- HZE_nte_model_coef['aa2']
kk1 <- HZE_nte_model_coef['kk1']

# The DER, = 0 at dose 0.
calibrated_nte_hazard_func <- function(dose, LET, coef) { # Calibrated hazard function.
  return(coef[1] * LET * dose * exp( - coef[2] * LET)
         + (1 - exp( - phi * dose)) * coef[3])
}

calibrated_HZE_nte_der <- function(dose, LET, coef = HZE_nte_model_coef) { # Calibrated HZE NTE DER.
  return(1 - exp( - calibrated_nte_hazard_func(dose, LET, coef)))
}


#================ HZE/TE MODEL =================#
# (TE = targeted effects only). This chunk runs and gives good results.
HZE_te_model <- nls( # Calibrating parameters in a TE only model.
  Prev ~ Y_0 + (1 - exp ( - (aate1 * LET * dose * exp( - aate2 * LET)))),
  data = HZE_data,
  weights = NWeight,
  start = list(aate1 = .00009, aate2 = .01))

summary(HZE_te_model, correlation = TRUE) # Parameter values & accuracy.
vcov(HZE_te_model) # Variance-covariance matrix.
HZE_te_model_coef <- coef(HZE_te_model) # Calibrated central values of the 2 parameters.

# The DER, = 0 at dose 0.
calibrated_te_hazard_func <- function(dose, LET, coef) { # Calibrated hazard function.
  return(coef[1] * LET * dose * exp( - coef[2] * LET))
}

calibrated_HZE_te_der <- function(dose, LET, coef = HZE_te_model_coef) {
  return(1 - exp( - calibrated_te_hazard_func(dose, LET, coef))) # Calibrated HZE TE one-ion DER.
}

#==== LIGHT ION, LOW Z (<= 3), LOW LET MODEL ===#

low_LET_model <- nls(
  Prev ~ Y_0 + 1 - exp( - alpha_low * dose), # alpha is used throughout radioiology for dose coefficients.
  data = low_LET_data,
  weights = NWeight,
  start = list(alpha_low = .005))

summary(low_LET_model, correlation = TRUE)
low_LET_model_coef <- coef(low_LET_model) # Calibrated central values of the parameter.

# Calibrated Low LET model. Use L = 0, but maybe later will use small L > 0.
calibrated_low_LET_der <- function(dose, LET, alph_low = low_LET_model_coef[1]) {
  return(1 - exp( - alph_low * dose))
}

# Slope dE/dd of the low LET, low Z model.
low_LET_slope <- function(dose, LET) {
  low_LET_model_coef * exp( - low_LET_model_coef * dose)
}


#=========================== INFORMATION CRITERION ============================#
info_crit_table <- cbind(AIC(HZE_te_model, HZE_nte_model),
                         BIC(HZE_te_model, HZE_nte_model))
print(info_crit_table)

##=================== Cross validation ====================#
# Seperate Data into 8 blocks, i.e. test/training sets:
data_len <- 1:length(HZE_data)
O_350 <- dplyr::select(filter(HZE_data, Beam == "O"), data_len) #RKS added 3 oxygen points 5/17/2019)
Ne_670 <- dplyr::select(filter(HZE_data, Beam == "Ne"), data_len)
Si_260 <- dplyr::select(filter(HZE_data, Beam == "Si"), data_len)
Ti_1000 <- dplyr::select(filter(HZE_data, Beam == "Ti"), data_len)
Fe_600 <- dplyr::select(filter(HZE_data, LET == 193), data_len)
Fe_350 <- dplyr::select(filter(HZE_data, LET == 253), data_len)
Nb_600 <- dplyr::select(filter(HZE_data, Beam == "Nb"), data_len)
La_593 <- dplyr::select(filter(HZE_data, Beam == "La"), data_len)
set_list <- list(O_350, Ne_670, Si_260, Ti_1000,
                 Fe_600, Fe_350, Nb_600, La_593)
actual_prev <- HZE_data$Prev

#============================ DEVELOPER FUNCTIONS =============================#
.test_runtime <- function(f, ...) { # Naive runtime check
  start_time <- Sys.time()
  f(...)
  Sys.time() - start_time
}
