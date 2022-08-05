#' Fit normal outcome model in TMB
#' 
#' Fit normal outcome model in TMB. Currently can only fit a space-only model
#' at a single time point. Linear predictor includes a single intercept and
#' a BYM2 sptial random effect, and an IID cluster-level random effect is included
#' if unit-level data is input (otherwise, no cluster-level random effect). 
#' Can fit benchmarked models, unbenchmarked models, or both.
#' Simultaneous benchmarking is performed via a second likelihood for the national level estimate. 
#' Produces posterior samples for fitted values, hyperparameters, 
#' random effects, and fixed effects. 
#' 
#' @param df a dataframe containing normal outcome data, and the following columns:
#' \itemize{
#'  \item \code{cluster}: cluster ID
#'  \item \code{y}: normally distributed outcome
#'  \item \code{region}: area
#' }
#' @param cluster the column in \code{binom_df} corresponding to cluster id. If data is area-level
#' then cluster should be set to NA. By default, cluster will be set to NA.
#' @param y the column in \code{binom_df} corresponding to the outcome
#' @param region the column in \code{binom_df} corresponding to region/area. Must be numeric and
#' start from 1.
#' @param SE_i a vector of fixed standard errors for use in the normal likelihood. Length of
#' \code{SE_i} must be equal to the number of regions.
#' @param hiv_adj An optional log offset in time to include in the linear predictor. Defaults to 
#' no log offset. 
#' @param natl a vector of national level estimates, arranged in order of time
#' @param natl_sd a vector of standard deviations for national level estimates, arranged in order of
#' time
#' @param pop_weights a vector of population weights for use in the benchmarking constraint. Must
#' sum to one at each time point, and be in order arrange(region)
#' @param Q_struct_space An ICAR precision matrix. Should be unscaled, as scaling will happen 
#' internally.
#' @param intercept_pri Prior specification for the intercept. Defaults to c(0, 31.62278), corresponding
#' to the default prior for the intercept in INLA, with mean 0 and precision 0.001. Must be
#' a vector of length 2, with specificaiton c(mean, sd) for a Normal distribution. Currently only 
#' an option for the area-level, unbenchmarked model.
#' @param nsamp Number of posterior samples to take from joint posterior. Defaults to 1000
#' @param benched A string, either \code{"benched"}, \code{"unbenched"}, or \code{"both"}, determining
#' whether to fit a benchmarked model, and unbenchmarked model, or both. Defaults to \code{"unbenched"}.
#' @param expit_outcome A boolean for whether or not the normally distributed outcome needs to be
#' expit-ed in order to be on the same scale as national estimates. This parameter will only affect
#' benchmarked models. If TRUE, linear predictors will be expit-ed before being multiplied by population
#' weights in the benchmarking constraint. This parameter must be specified by the user if 
#' \code{benched = "benched"}.
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
#' @export fit_normal
fit_normal <- function(df, 
                       cluster = NA,
                       y = "y", 
                       region = "admin1", 
                       SE_i,
                       hiv_adj = NA,
                       natl = NULL,
                       natl_sd = NULL,
                       pop_weights = NULL,
                       Q_struct_space = NULL, 
                       intercept_pri = c(0, 31.62278),
                       nsamp = 1000,
                       benched = "unbenched",
                       expit_outcome = NULL) {
  
  # quick fix - taylor, fix this later
  binom_df <- df
  
  # error handling
  
  # alpha_pri must be length(2) and numeric
  alpha_pri <- intercept_pri
  if (!is.numeric(alpha_pri)) {
    stop("prior for intercept must be numeric")
  }
  if (length(alpha_pri) != 2) {
    stop("prior for intercept must be length(2), in form (mean, sd)")
  }
  
  # region, time, age, age_group, cluster, y, Ntrials must all be in binom_df
  if (!(region %in% colnames(binom_df))) {
    stop("region must be a column in binom_df")
  }
  if (!is.na(cluster)) {
    if (!(cluster %in% colnames(binom_df))) {
      stop("cluster must be a column in binom_df")
    }
  }
  if (!(y %in% colnames(binom_df))) {
    stop("y must be a column in binom_df")
  }
  
  # region, time, age, age_group, cluster must all be numeric and start from 1
  if (!(is.numeric(binom_df[,region]) & (min(binom_df[,region]) == 1))) {
    stop("region must be numeric and start from 1")
  }
  if (!is.na(cluster)) {
    if (!(is.numeric(binom_df[,cluster]) & (min(binom_df[,cluster]) == 1))) {
      stop("cluster must be numeric and start from 1")
    }
  }
  if (!(is.numeric(binom_df[,y]))) {
    stop("y must be numeric")
  }
  
  # make sure length(SE_i) is the same as length(unique(binom_df[,region]))
  if (length(SE_i) != length(unique(binom_df[,region]))) {
    stop("length of SE_i must equal the unique number of regions in df")
  }
  
  # if benched = "benched", make sure that expit_outcome is specified
  if (benched == "benched" | benched == "both") {
    if (is.null(expit_outcome)) {
      stop("expit_outcome must be specified by the user for benchmarked models")
    }
  } else {
    if (!is.null(expit_outcome)) {
      message("expit_outcome parameter does not affect unbenchmarked models")
    }
  }
  
  # check that expit_outcome is boolean
  if (!is.null(expit_outcome)) {
    if (!is.logical(expit_outcome)) {
      stop("expit_outcome must be TRUE or FALSE")
    }
  }
  
  # relabel clusters from 1:length(unique(cluster))
  if (!is.na(cluster)) {
    binom_df$cluster <- binom_df[,cluster]
    suppressMessages(binom_df <- data.frame(cluster = sort(unique(binom_df[,cluster])), cluster_new = 1:length(unique(binom_df[,cluster]))) %>%
                       right_join(binom_df))
    cluster <- "cluster_new"
  }
  
  # create (scaled) Q_struct_space
  
  Q_struct_space <- INLA::inla.scale.model(Q_struct_space, 
                                           constr = list(A = matrix(1, nrow = 1, ncol = nrow(Q_struct_space)), e = 0))
  
  # determine which model to fit
  
  # unit-level space-only model
  if (!is.na(cluster)) {
    if (benched == "unbenched") {
      message("Unbenched Normal model")
      unbenched_res <- tmb_normal_intercepts_bym2(binom_df = binom_df,
                                                  cluster = cluster,
                                                  y = y,
                                                  region = region,
                                                  SE_i = SE_i,
                                                  hiv_adj = hiv_adj,
                                                  Q_struct_space = Q_struct_space,
                                                  nsamp = nsamp)
    } else if (benched == "benched") {
      message("Benched Normal model")
      benched_res <- tmb_normal_intercepts_bym2_benched(binom_df = binom_df,
                                                        cluster = cluster,
                                                        y = y,
                                                        region = region,
                                                        SE_i = SE_i,
                                                        hiv_adj = hiv_adj,
                                                        natl = natl,
                                                        natl_sd = natl_sd,
                                                        pop_weights = pop_weights,
                                                        Q_struct_space = Q_struct_space,
                                                        expit_outcome = expit_outcome,
                                                        nsamp = nsamp)
    } else {
      message("Unbenched Normal model")
      unbenched_res <- tmb_normal_intercepts_bym2(binom_df = binom_df,
                                                  cluster = cluster,
                                                  y = y,
                                                  region = region,
                                                  SE_i = SE_i,
                                                  hiv_adj = hiv_adj,
                                                  Q_struct_space = Q_struct_space,
                                                  nsamp = nsamp)
      
      message("Benched Normal model")
      benched_res <- tmb_normal_intercepts_bym2_benched(binom_df = binom_df,
                                                        cluster = cluster,
                                                        y = y,
                                                        region = region,
                                                        SE_i = SE_i,
                                                        hiv_adj = hiv_adj,
                                                        natl = natl,
                                                        natl_sd = natl_sd,
                                                        pop_weights = pop_weights,
                                                        Q_struct_space = Q_struct_space,
                                                        expit_outcome = expit_outcome,
                                                        nsamp = nsamp)
      
    }
  
  # area-level space-only model
  } else {
    if (benched == "unbenched") {
      message("Unbenched Normal model")
      unbenched_res <- tmb_normal_arealevel_intercepts_bym2(binom_df = binom_df,
                                                            y = y,
                                                            region = region,
                                                            SE_i = SE_i,
                                                            hiv_adj = hiv_adj,
                                                            Q_struct_space = Q_struct_space,
                                                            alpha_pri = alpha_pri,
                                                            nsamp = nsamp)
    } else if (benched == "benched") {
      message("Benched Normal model")
      benched_res <- tmb_normal_arealevel_intercepts_bym2_benched(binom_df = binom_df,
                                                                  y = y,
                                                                  region = region,
                                                                  SE_i = SE_i,
                                                                  hiv_adj = hiv_adj,
                                                                  natl = natl,
                                                                  natl_sd = natl_sd,
                                                                  pop_weights = pop_weights,
                                                                  Q_struct_space = Q_struct_space,
                                                                  expit_outcome = expit_outcome,
                                                                  nsamp = nsamp)
    } else {
      message("Unenched Normal model")
      unbenched_res <- tmb_normal_arealevel_intercepts_bym2(binom_df = binom_df,
                                                            y = y,
                                                            region = region,
                                                            SE_i = SE_i,
                                                            hiv_adj = hiv_adj,
                                                            Q_struct_space = Q_struct_space,
                                                            alpha_pri = alpha_pri,
                                                            nsamp = nsamp)
      
      message("Benched Normal model")
      benched_res <- tmb_normal_arealevel_intercepts_bym2_benched(binom_df = binom_df,
                                                                  y = y,
                                                                  region = region,
                                                                  SE_i = SE_i,
                                                                  hiv_adj = hiv_adj,
                                                                  natl = natl,
                                                                  natl_sd = natl_sd,
                                                                  pop_weights = pop_weights,
                                                                  Q_struct_space = Q_struct_space,
                                                                  expit_outcome = expit_outcome,
                                                                  nsamp = nsamp)
      
    }
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





