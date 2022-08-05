#' Fit Unbenchmarked normal outcome model in TMB
#' 
#' Fit unbenchmarked model with BYM2 in space, single intercept. Produces posterior 
#' samples for fitted values, hyperparameters, random effects, and fixed effects. 
#' 
#' @param binom_df a dataframe containing the following columns: (note to self, change this parameter
#' name since we're not working with binomial data here)
#' \itemize{
#'  \item \code{y}: normally distributed outcome
#'  \item \code{region}: area
#' }
#' @param y the column in \code{binom_df} corresponding to the outcome
#' @param region the column in \code{binom_df} corresponding to region/area
#' @param SE_i a vector of fixed standard errors for use in the normal likelihood. Length of
#' \code{SE_i} must be equal to the number of regions.
#' @param hiv_adj An optional log offset in time to include in the linear predictor. Defaults to 
#' no log offset. 
#' @param Q_struct_space An ICAR precision matrix. Should be scaled.
#' @param alpha_pri Prior specification for the intercept. Defaults to c(0, 31.62278), corresponding
#' to the default prior for the intercept in INLA, with mean 0 and precision 0.001. Must be
#' a vector of length 2, with specificaiton c(mean, sd) for a Normal distribution.
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
tmb_normal_arealevel_intercepts_bym2 <- function(binom_df, 
                                       y = "y", 
                                       region = "admin1", 
                                       SE_i,
                                       hiv_adj = NA,
                                       Q_struct_space,
                                       alpha_pri = c(0, 31.62278),
                                       nsamp = 1000) {
  
  N <- dim(Q_struct_space)[1]
  
  # if hiv_adj is NA, make it 1
  if (is.na(hiv_adj)) {
    hiv_adj <- 1
  }
  
  # get eigenvalues of Q_struct_space
  # gamma <- eigen(Q_struct_space)$values
  # gamma_tilde <- 1/gamma
  # gamma_tilde[length(gamma_tilde)] <- 0
  
  # define data object
  t.data <- list(model = "normal_arealevel_intercepts_bym2",
                 N = N,
                 y_ic = binom_df[,y],
                 SE_i = SE_i,
                 hiv_adj = hiv_adj,
                 region_id = binom_df[,region],
                 Q_struc_space = Q_struct_space,
                 alpha_pri = alpha_pri,
                 risk_region_id = 1:length(unique(binom_df[,region])))
  
  # starting vals
  t.params <- list(alpha = 0,
                   space_logit_phi = 0,
                   space_log_tau = 0,
                   Epsilon_s = matrix(0, ncol = 1, nrow = 2 * N))
  
  # rand effs
  t.rand <- c('Epsilon_s')
  
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
  # cluster_effect_idx <- which(names(mu) == "cluster_effect")
  
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
  # re_list$cluster_effect <- t.draws[cluster_effect_idx, ]
  
  ## separate out the tmb_draws
  # NOTE: t.parnames will be the same as names(mu)
  t.param.idx <- grep("alpha", t.parnames)
  t.space.phi.idx <- grep("space_logit_phi", t.parnames)
  t.space.tau.idx <- grep("space_log_tau", t.parnames)
  # t.cluster.tau.idx <- grep("cluster_log_tau", t.parnames)
  
  # make a list containing all draws from parameters
  param_list <- list()
  param_list$intercepts <- t.draws[t.param.idx, ]
  param_list$space_phi <- expit(t.draws[t.space.phi.idx,])
  param_list$space_tau <- exp(t.draws[t.space.tau.idx,])
  # param_list$cluster_tau <- exp(t.draws[t.cluster.tau.idx,])
  
  # estimates of tmb field
  # below is the total bym2 term for space, not just the structured part
  t.space.total.idx <- head(grep("Epsilon_s", t.parnames), N)
  # id's for each age intercept is t.param.idx
  
  # combine draws for linear predictor to get space/time/age risks
  # this bit adds in the total for space, intercepts for each age group
  e.fitted <- t.draws[rep(t.param.idx, length(t.space.total.idx)),] +
    t.draws[t.space.total.idx,] 
  
  # matrix of samples of fitted values on probability scale
  fitted_mat <- e.fitted
  
  # return fitted mat, re_list, param_list
  return(list(fitted_mat = fitted_mat,
              re_list = re_list,
              param_list = param_list,
              runtime = end_time - begin_time))
}






