% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_binary.R
\name{fit_binary}
\alias{fit_binary}
\title{Fit binary outcome model in TMB}
\usage{
fit_binary(
  binom_df,
  cluster = "cluster",
  y = "y",
  Ntrials = "N",
  region = "admin1",
  hiv_adj = NA,
  natl = NULL,
  natl_sd = NULL,
  pop_weights = NULL,
  Q_struct_space = NULL,
  intercept_pri = c(0, 31.62278),
  nsamp = 1000,
  benched = "unbenched",
  family = "binomial"
)
}
\arguments{
\item{binom_df}{a dataframe containing binomial counts, and the following columns:
\itemize{
\item \code{cluster}: cluster ID
\item \code{y}: deaths
\item \code{Ntrials}: total number of person-months in this \code{age}, \code{cluster}, and \code{region}
\item \code{region}: area
}}

\item{cluster}{the column in \code{binom_df} corresponding to cluster id}

\item{y}{the column in \code{binom_df} corresponding to deaths}

\item{Ntrials}{the column in \code{binom_df} corresponding to total number of person-months}

\item{region}{the column in \code{binom_df} corresponding to region/area. Must be numeric and
start from 1.}

\item{hiv_adj}{An optional log offset in time to include in the linear predictor. Defaults to
no log offset.}

\item{natl}{a vector of national level estimates, arranged in order of time}

\item{natl_sd}{a vector of standard deviations for national level estimates, arranged in order of
time}

\item{pop_weights}{a vector of population weights for use in the benchmarking constraint. Must
sum to one at each time point, and be in order arrange(region)}

\item{Q_struct_space}{An ICAR precision matrix. Should be unscaled, as scaling will happen
internally.}

\item{intercept_pri}{Prior specification for the intercept. Defaults to c(0, 31.62278), corresponding
to the default prior for the intercept in INLA, with mean 0 and precision 0.001. Must be
a vector of length 2, with specificaiton c(mean, sd) for a Normal distribution. Currently only
an option for unbenchmarked models.}

\item{nsamp}{Number of posterior samples to take from joint posterior. Defaults to 1000}

\item{benched}{A string, either \code{"benched"}, \code{"unbenched"}, or \code{"both"}, determining
whether to fit a benchmarked model, and unbenchmarked model, or both. Defaults to \code{"unbenched"}.}

\item{family}{A string, either \code{"binomial"} or \code{"betabinomial"}, specifying which likelihood
to use for the clustered observations. If \code{"binomial"} is specified, does the lono-binomial
correction, as described in Dong and Wakefield (2021). Defaults to \code{"binomial"}.}
}
\value{
A list containing:
\itemize{
\item fitted_mat: a matrix of posterior samples of fitted values in order arrange(region, time)
\item re_list: a list contaning matrices of posterior samples for each random effect term
\item param_list: a list containing matrices of posterior samples for fixed effects and hyperparameters
\item runtime: the time it took to fit the model in TMB and get samples from the joint posterior
}
If \code{benched = "both"}, a list of two will be returned containing the above list for both
benchmarked and unbenchmarked models.
}
\description{
Fit binary outcome model in TMB. Currently can only fit a space-only model
at a single time point. Linear predictor includes a single intercept and
a BYM2 sptial random effect (and a IID cluster-level random effect if a lono-binomial
likelihood is chosen). Can fit benchmarked models, unbenchmarked models, or both.
Simultaneous benchmarking is performed via a second likelihood for the national level estimate.
Produces posterior samples for fitted values, hyperparameters,
random effects, and fixed effects. All models use either a lono-binomial likelihood, described
in Dong and Wakefield (2021) or a betabinomial likelihood at the cluster level.
}
\author{
Taylor Okonek
}
