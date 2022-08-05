#' Fit Benchmarked binary outcome model in TMB
#' 
#' Fit benchmarked model with BYM2 in space, single intercept. Only allows for a single
#' benchmark, and therefore assumes one time point.Produces posterior samples for fitted values, 
#' hyperparameters, random effects, and fixed effects. Simultaneous benchmarking is performed 
#' via a second likelihood for the national level estimate.
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
#' @param region the column in \code{binom_df} corresponding to region/area
#' @param hiv_adj An optional log offset in time to include in the linear predictor. Defaults to 
#' no log offset. 
#' @param natl a national level estimate
#' @param natl_sd the standard deviation for the national level estimate
#' @param pop_weights a vector of population weights for use in the benchmarking constraint. Must
#' sum to one, and be in order arrange(region)
#' @param Q_struct_space An ICAR precision matrix. Should be scaled.
#' @param nsamp Number of posterior samples to take from joint posterior. Defaults to 1000
#' @return A list containing: 
#' \itemize{
#' \item fitted_mat: a matrix of posterior samples of fitted values in order arrange(region, time)
#' \item re_list: a list contaning matrices of posterior samples for each random effect term
#' \item param_list: a list containing matrices of posterior samples for fixed effects and hyperparameters
#' \item runtime: the time it took to fit the model in TMB and get samples from the joint posterior
#' } 
#' 
#' @author Taylor Okonek
#' @keywords internal
#' @noRd
tmb_binary_intercepts_bym2_benched <- function(binom_df, 
                                               cluster = "cluster",
                                               y = "y", 
                                               Ntrials = "N", 
                                               region = "admin1", 
                                               hiv_adj = NA,
                                               natl,
                                               natl_sd,
                                               pop_weights,
                                               Q_struct_space,
                                               nsamp = 1000) {
  
  N <- dim(Q_struct_space)[1]
  
  # if hiv_adj is NA, make it 1
  if (is.na(hiv_adj)) {
    hiv_adj <- 1
  }
   
  # define data object
  t.data <- list(model = "binary_intercepts_bym2_benched",
                 N = N,
                 nclusts = length(unique(binom_df[,cluster])),
                 nobs = nrow(binom_df),
                 y_ic = binom_df[,y],
                 N_ic = binom_df[,Ntrials],
                 hiv_adj = hiv_adj,
                 region_id = binom_df[,region],
                 cluster_id = binom_df[,cluster],
                 natl = natl,
                 natl_sd = natl_sd,
                 pop_weights = pop_weights,
                 Q_struc_space = Q_struct_space,
                 risk_region_id = 1:length(unique(binom_df[,region])))
  
  # starting vals
  t.params <- list(alpha = 0,
                   space_logit_phi = 0,
                   space_log_tau = 0,
                   cluster_log_tau = 0,
                   Epsilon_s = matrix(0, ncol = 1, nrow = 2 * N),
                   cluster_effect = matrix(0, ncol = 1, nrow = length(unique(binom_df[,cluster]))))
  
  # rand effs
  t.rand <- c('Epsilon_s', 
              "cluster_effect")
  
  ## and lastly, we can NULL out params that are in the c++ template but
  ## which won't be used in this run. this allows us to build up a more
  ## complicated template that can take different options. named items
  ## in the list which as set to factor(NA) will be left out
  ADmap <- list()
  
  ## fit in TMB
  
  # make AD functional
  obj <- TMB::MakeADFun(data  = t.data,
                        parameters = t.params,
                        random     = t.rand,
                        map        = ADmap,
                        hessian    = TRUE,
                        DLL        = "stbench_TMBExports")
  
  message("Model fitting...")
  
  begin_time <- Sys.time()
  
  # call the optimizer
  opt0 <- try(do.call("nlminb",
                      list(start     = obj$par,
                           objective = obj$fn,
                           gradient  = obj$gr,
                           lower     = rep(-10, length(obj$par)),
                           upper     = rep( 10, length(obj$par)),
                           control   = list(trace=1))))
  
  # report the estimates
  # calculates standard deviations of all model parameters
  SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                       getReportCovariance = TRUE,
                       bias.correct = TRUE,
                       bias.correct.control = list(sd = TRUE))
  
  end_time <- Sys.time()
  
  message("Post-processing...")
  
  # get all values for things with REPORT(things); in the c++ code
  R0 <- obj$report()
  
  ####### Organize results from TMB
  
  ## take TMB draws
  
  # first obtain point estimates (means) for fixed and random effects
  mu <- c(SD0$par.fixed, SD0$par.random)
  
  # first, try to get the cholesky decomp of the joint precision (fixed and random effects)
  # note that length(mu) = nrow(L) = ncol(L)
  # L <- try(suppressWarnings(Cholesky(SD0$jointPrecision, super = T)), silent = TRUE)
  
  # take the draws
  # this calls on the function rmvnorm_prec, which samples from a multivariate normal 
  # using a cholesky decomposition - NOTE: look into this more
  # returns a list containing two matrices: one with unconstrained draws, one with constrained draws
  # get ids for which values of mu correspond to the structured bym2 components in space and time
  struct_space_idx <- which(names(mu) == "Epsilon_s") %>% tail(N)
  struct_space_both_idx <- which(names(mu) == "Epsilon_s")
  cluster_effect_idx <- which(names(mu) == "cluster_effect")
  
  # create list of constraint matrices for these terms
  A.mat.list <- list()
  A.mat.list[[1]] <- matrix(1, nrow = 1, ncol = N) # space
  
  t.draws <- multiconstr_prec(mu = mu,
                              prec = SD0$jointPrecision,
                              n.sims = nsamp,
                              #constrain = FALSE,
                              constrain.idx.list = list(struct_space_idx),
                              A.mat.list = A.mat.list)
  
  # take the constrained draws
  t.draws <- t.draws$x.c
  
  # get list of draws for each random effect
  re_list <- list()
  re_list$struct_space <- t.draws[struct_space_both_idx, ]
  t.parnames <- c(names(SD0$par.fixed), names(SD0$par.random))
  re_list$cluster_effect <- t.draws[cluster_effect_idx, ]
  
  ## separate out the tmb_draws
  # NOTE: t.parnames will be the same as names(mu)
  t.param.idx <- grep("alpha", t.parnames)
  t.space.phi.idx <- grep("space_logit_phi", t.parnames)
  t.space.tau.idx <- grep("space_log_tau", t.parnames)
  t.cluster.tau.idx <- grep("cluster_log_tau", t.parnames)
  
  # make a list containing all draws from parameters
  param_list <- list()
  param_list$intercepts <- t.draws[t.param.idx, ]
  param_list$space_phi <- expit(t.draws[t.space.phi.idx,])
  param_list$space_tau <- exp(t.draws[t.space.tau.idx,])
  param_list$cluster_tau <- exp(t.draws[t.cluster.tau.idx,])
  
  # estimates of tmb field
  # below is the total bym2 term for space, not just the structured part
  t.space.total.idx <- head(grep("Epsilon_s", t.parnames), N)
  # id's for each age intercept is t.param.idx
  
  # combine draws for linear predictor to get space/time/age risks
  # this bit adds in the total for space, intercepts for each age group
  e.fitted <- t.draws[rep(t.param.idx, length(t.space.total.idx)),] +
    t.draws[t.space.total.idx,] 
  
  # lono correction
  h <- 16 * sqrt(3) / (15 * pi)
  cluster_vars <- 1/param_list$cluster_tau
  denom_samps <- sqrt(1 + h^2 * cluster_vars)

  # correct fitted samps
  for (i in 1:nrow(e.fitted)) {
    e.fitted[i,] <- e.fitted[i,] / denom_samps
  }
  
  # matrix of samples of fitted values on probability scale
  fitted_mat <- expit(e.fitted)
  
  # return fitted mat, re_list, param_list
  return(list(fitted_mat = fitted_mat,
              re_list = re_list,
              param_list = param_list,
              runtime = end_time - begin_time))
}






