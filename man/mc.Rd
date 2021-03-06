% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc.R
\name{mc}
\alias{mc}
\title{Runs the Monte Carlo method on a baseline no-synergy/antagonism
       mixture DER.}
\usage{
mc(dose, LET, ratios, E, models, dE = NULL, n = 200, vcov = TRUE,
  interval_length = 0.95, seed = 100, check = TRUE, ...)
}
\arguments{
\item{dose}{Numeric vector of all total dose values to be evaluated.}

\item{LET}{Numeric vector of all LET values, must be length n.}

\item{ratios}{Numeric vector of dose proportions applied on component DERs.}

\item{E}{Vector of DER functions}

\item{models}{List of nls models.}

\item{n}{Numeric integer of the number of samples to be drawn.}

\item{vcov}{Boolean for assessing inter-parameter correlation.}

\item{interval_length}{Numeric double of the confidence interval width.}

\item{seed}{Numeric value for pseudorandom generators.}

\item{...}{Optional arguments to DER functions in E. Should be constant
across all supplied E.}
}
\value{
Named list representing lower and upper bounds for a Monte Carlo
        confidence interval for a mixture DER over an interval of doses.
}
\description{

}
\details{
Corresponding elements of ratios, LET should be associated with the
         same DER.
}
\examples{

ion_data <- load_ion_data("one_ion.csv")
HZE_nte_der <- make_der(ion_data, TRUE, TRUE)
check_der(HZE_nte_der)
# Iron and silicon two-ion mixture.
mc(dose = 0:100, LET = c(193, 70), ratios = c(1/3, 2/3),
   E = rep(c(HZE_nte_der), 2),
   models = rep(list(HZE_nte_model), 2), n = 20, vcov = TRUE)

}
\author{
Yimin Lin, Edward Greg Huang <eghuang@berkeley.edu>
}
