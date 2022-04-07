#' Fit binary outcome model in TMB
#' 
#' Fit binary outcome model in TMB. Currently can only fit a space-only model
#' at a single time point. Linear predictor includes a single intercept and
#' a BYM2 sptial random effect (and a IID cluster-level random effect if a lono-binomial
#' likelihood is chosen). Can fit benchmarked models, unbenchmarked models, or both.
#' Simultaneous benchmarking is performed via a second likelihood for the national level estimate. 
#' Produces posterior samples for fitted values, hyperparameters, 
#' random effects, and fixed effects. All models use either a lono-binomial likelihood, described
#' in Dong and Wakefield (2021) or a betabinomial likelihood at the cluster level.
#' 
#' @param binom_df a dataframe containing binomial counts, and the following columns:
#' \itemize{
#'  \item \code{cluster}: cluster ID
#'  \item \code{y}: deaths
#'  \item \code{Ntrials}: total number of person-months in this \code{age}, \code{cluster}, and \code{region}
#'  \item \code{region}: area
#' }
#' @param cluster the column in \code{binom_df} corresponding to cluster id
#' @param y the column in \code{binom_df} corresponding to deaths
#' @param Ntrials the column in \code{binom_df} corresponding to total number of person-months
#' @param region the column in \code{binom_df} corresponding to region/area. Must be numeric and
#' start from 1.
#' @param hiv_adj An optional log offset in time to include in the linear predictor. Defaults to 
#' no log offset. 
#' @param natl a vector of national level estimates, arranged in order of time
#' @param natl_sd a vector of standard deviations for national level estimates, arranged in order of
#' time
#' @param pop_weights a vector of population weights for use in the benchmarking constraint. Must
#' sum to one at each time point, and be in order arrange(region)
#' @param Q_struct_space An ICAR precision matrix. Should be unscaled, as scaling will happen 
#' internally.
#' @param nsamp Number of posterior samples to take from joint posterior. Defaults to 1000
#' @param benched A string, either \code{"benched"}, \code{"unbenched"}, or \code{"both"}, determining
#' whether to fit a benchmarked model, and unbenchmarked model, or both. Defaults to \code{"unbenched"}.
#' @param family A string, either \code{"binomial"} or \code{"betabinomial"}, specifying which likelihood
#' to use for the clustered observations. If \code{"binomial"} is specified, does the lono-binomial 
#' correction, as described in Dong and Wakefield (2021). Defaults to \code{"binomial"}.
#' @return A list containing: 
#' \itemize{
#' \item fitted_mat: a matrix of posterior samples of fitted values in order arrange(region, time)
#' \item re_list: a list contaning matrices of posterior samples for each random effect term
#' \item param_list: a list containing matrices of posterior samples for fixed effects and hyperparameters
#' \item runtime: the time it took to fit the model in TMB and get samples from the joint posterior
#' } 
#' If \code{benched = "both"}, a list of two will be returned containing the above list for both
#' benchmarked and unbenchmarked models.
#' 
#' @author Taylor Okonek
#' @export fit_binary
fit_binary <- function(binom_df, 
                       cluster = "cluster",
                       y = "y", 
                       Ntrials = "N", 
                       region = "admin1", 
                       hiv_adj = NA,
                       natl = NULL,
                       natl_sd = NULL,
                       pop_weights = NULL,
                       Q_struct_space = NULL, 
                       nsamp = 1000,
                       benched = "unbenched",
                       family = "binomial") {
  # error handling
  
  # betabinomial currently unsupported
  if (!(family %in% c("binomial", "betabinomial"))) {
    stop("family not currently supported")
  }
  
  # region, time, age, age_group, cluster, y, Ntrials must all be in binom_df
  if (!(region %in% colnames(binom_df))) {
    stop("region must be a column in binom_df")
  }
  if (!(cluster %in% colnames(binom_df))) {
    stop("cluster must be a column in binom_df")
  }
  if (!(y %in% colnames(binom_df))) {
    stop("y must be a column in binom_df")
  }
  if (!(Ntrials %in% colnames(binom_df))) {
    stop("Ntrials must be a column in binom_df")
  }
  
  # region, time, age, age_group, cluster must all be numeric and start from 1
  if (!(is.numeric(binom_df[,region]) & (min(binom_df[,region]) == 1))) {
    stop("region must be numeric and start from 1")
  }
  if (!(is.numeric(binom_df[,cluster]) & (min(binom_df[,cluster]) == 1))) {
    stop("cluster must be numeric and start from 1")
  }
  if (!(is.numeric(binom_df[,y]))) {
    stop("y must be numeric")
  }
  if (!(is.numeric(binom_df[,Ntrials]))) {
    stop("Ntrials must be numeric")
  }
  
  # relabel clusters from 1:length(unique(cluster))
  binom_df$cluster <- binom_df[,cluster]
  suppressMessages(binom_df <- data.frame(cluster = sort(unique(binom_df[,cluster])), cluster_new = 1:length(unique(binom_df[,cluster]))) %>%
                     right_join(binom_df))
  cluster <- "cluster_new"
  
  # create (scaled) Q_struct_space
  
  Q_struct_space <- INLA::inla.scale.model(Q_struct_space, 
                                           constr = list(A = matrix(1, nrow = 1, ncol = nrow(Q_struct_space)), e = 0))
  
  # determine which model to fit
  
  if (family == "binomial") {
    if (benched == "unbenched") {
      message("Unbenched Binomial model")
      unbenched_res <- tmb_binary_intercepts_bym2(binom_df,
                                                  cluster,
                                                  y,
                                                  Ntrials,
                                                  region,
                                                  hiv_adj,
                                                  Q_struct_space,
                                                  nsamp)
    } else if (benched == "benched") {
      message("Benched Binomial model")
      benched_res <- tmb_binary_intercepts_bym2_benched(binom_df,
                                                        cluster,
                                                        y,
                                                        Ntrials,
                                                        region,
                                                        hiv_adj,
                                                        natl,
                                                        natl_sd,
                                                        pop_weights,
                                                        Q_struct_space,
                                                        nsamp)
    } else {
      message("Unenched Binomial model")
      unbenched_res <- tmb_binary_intercepts_bym2(binom_df,
                                                  cluster,
                                                  y,
                                                  Ntrials,
                                                  region,
                                                  hiv_adj,
                                                  Q_struct_space,
                                                  nsamp)
      
      message("Benched Binomial model")
      benched_res <- tmb_binary_intercepts_bym2_benched(binom_df,
                                                        cluster,
                                                        y,
                                                        Ntrials,
                                                        region,
                                                        hiv_adj,
                                                        natl,
                                                        natl_sd,
                                                        pop_weights,
                                                        Q_struct_space,
                                                        nsamp)
      
    }
  } else {
    stop("betabinomial models not currently implemented")
  }
  
  # return appropriate results
  if (benched == "unbenched") {
    return(unbenched_res)
  } else if (benched == "benched") {
    return(benched_res)
  } else {
    return(list(benched_res = benched_res,
                unbenched_res = unbenched_res))
  }
  
}





