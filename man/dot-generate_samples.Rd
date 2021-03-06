% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc.R
\name{.generate_samples}
\alias{.generate_samples}
\title{Generates mixture no-synergy/antagonism DER samples with
       parameters drawn from a Gaussian distribution.}
\usage{
.generate_samples(n, dose, LET, ratios, E, dE, models, vcov, ...)
}
\arguments{
\item{n}{Numeric integer of the number of samples to be drawn.}

\item{dose}{Numeric vector of all total dose values to be evaluated.}

\item{LET}{Numeric vector of all LET values, must be length n.}

\item{ratios}{Numeric vector of dose proportions applied on component DERs.}

\item{vcov}{Boolean for assessing inter-parameter correlation.}

\item{...}{Optional arguments to DER functions in E.}
}
\value{
Numeric vector of sample mixture baseline DERs evaluated at the
        given doses.
}
\description{

}
