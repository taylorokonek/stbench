#' Fit U5MR model in TMB
#' 
#' Fit U5MR model with age-specific intercepts, BYM2 in space.
#' Produces posterior samples for fitted values, hyperparameters, 
#' random effects, and fixed effects.
#' 
#' @param binom_df a dataframe containing binomial counts, and the following columns:
#' \itemize{
#'  \item \code{cluster}: cluster ID
#'  \item \code{y}: deaths
#'  \item \code{Ntrials}: total number of person-months in this \code{age}, \code{cluster}, and \code{region}
#'  \item \code{region}: area
#'  \item \code{age}: id for each age group (must be numeric) you want intercepts for
#'  \item \code{age_group}: id for each age group you want RW2s for 
#' }
#' @param cluster the column in \code{binom_df} corresponding to cluster id
#' @param y the column in \code{binom_df} corresponding to deaths
#' @param Ntrials the column in \code{binom_df} corresponding to total number of person-months
#' @param region the column in \code{binom_df} corresponding to region/area
#' @param age the column in \code{binom_df} corresponding to age id, used for age-specific intercepts
#' @param age_group the column in \code{binom_df} corresponding to age group id, used for
#' age-group-specific RWs in time
#' @param age_n the number of months in each age group specified by the \code{age} column in 
#' \code{binom_df}. The length of \code{age_n} must be equal to the number of unique age groups
#' specified in \code{age}.
#' @param hiv_adj An optional log offset in time to include in the linear predictor. Defaults to 
#' no log offset. If included, the length of \code{hiv_adj} must be equal to the number of unique
#' time points specified in \code{time}.
#' @param Q_struct_space An ICAR precision matrix. Should be scaled.
#' @param nsamp Number of posterior samples to take from joint posterior. Defaults to 1000
#' @return A list containing: 
#' \itemize{
#' \item fitted_mat: a matrix of posterior samples of fitted values in order arrange(region)
#' \item re_list: a list containing matrices of posterior samples for each random effect term
#' \item param_list: a list containing matrices of posterior samples for fixed effects and hyperparameters
#' \item runtime: the time it took to fit the model in TMB and get samples from the joint posterior
#' } 
#' 
#' @author Taylor Okonek
#' @keywords internal
#' @noRd
tmb_u5mr_intercepts_bym2 <- function(binom_df, 
                                     cluster = "cluster",
                                     y = "y", 
                                     Ntrials = "N", 
                                     region = "admin1", 
                                     age = "age_id",
                                     age_n = c(1,11,12,12,12,12),
                                     age_group = "age_group_id",
                                     hiv_adj = NA,
                                     Q_struct_space, 
                                     nsamp) {
  
  N <- dim(Q_struct_space)[1]
  
  # if hiv_adj is NA, make it 1
  if (is.na(hiv_adj[1])) {
    hiv_adj <- 1
  }
  
  # create spaceage_id variable: 1:(N*A_m), in order arrange(region, age)
  spaceage_id_df <- expand.grid(region = 1:N,
                                    age = 1:length(unique(binom_df[,age]))) %>%
    arrange(region, age)
  spaceage_id_df$spaceage_id <- 1:nrow(spaceage_id_df)
  # left join with binom_df to get full spacetime_id variable
  binom_df$region <- binom_df[,region]
  binom_df$age <- binom_df[,age]
  binom_df <- binom_df %>% arrange(region, age)
  suppressMessages(spaceage_id <- left_join(binom_df, spaceage_id_df)$spaceage_id)
  
  # define data object
  t.data <- list(model = "u5mr_intercepts_bym2",
                 N = N,
                 A_m = length(unique(binom_df[,age])),
                 nclusts = length(unique(binom_df[,cluster])),
                 nobs = nrow(binom_df),
                 y_iac = binom_df[,y],
                 N_iac = binom_df[,Ntrials],
                 hiv_adj = hiv_adj,
                 region_id = binom_df[,region],
                 age_id = binom_df[,age],
                 spaceage_id = spaceage_id,
                 cluster_id = binom_df[,cluster],
                 Q_struc_space = Q_struct_space)
  
  # starting vals
  t.params <- list(alpha = matrix(0, ncol = 1, nrow = length(unique(binom_df[,age]))),
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
  # unstruct_time_idx <- which(names(mu) == "Epsilon_t") 
  struct_space_idx <- which(names(mu) == "Epsilon_s") %>% tail(N)
  struct_space_both_idx <- which(names(mu) == "Epsilon_s")
  # struct_time_idx1 <- which(names(mu) == "age_rw2_1")
  # struct_time_idx2 <- which(names(mu) == "age_rw2_2")
  # struct_time_idx3 <- which(names(mu) == "age_rw2_3")
  # spacetime_idx <- which(names(mu) == "Epsilon_st")
  cluster_effect_idx <- which(names(mu) == "cluster_effect")
  
  # create list of constraint matrices for these terms
  A.mat.list <- list()
  A.mat.list[[1]] <- matrix(1, nrow = 1, ncol = N) # space
  # A.mat.list[[2]] <- matrix(1, nrow = 1, ncol = S) # time - RW2
  # A.mat.list[[3]] <- matrix(1, nrow = 1, ncol = S) # time - RW2
  # A.mat.list[[4]] <- matrix(1, nrow = 1, ncol = S) # time - RW2
  # A.mat.list[[2]] <- rbind(matrix(1, nrow = 1, ncol = length(unique(binom_df[,time]))),
  #                          matrix(1:length(unique(binom_df[,time])), nrow = 1, ncol = length(unique(binom_df[,time])))) # time - RW2
  # A.mat.list[[3]] <- rbind(matrix(1, nrow = 1, ncol = length(unique(binom_df[,time]))),
  #                          matrix(1:length(unique(binom_df[,time])), nrow = 1, ncol = length(unique(binom_df[,time])))) # time - RW2
  # A.mat.list[[4]] <- rbind(matrix(1, nrow = 1, ncol = length(unique(binom_df[,time]))),
  #                          matrix(1:length(unique(binom_df[,time])), nrow = 1, ncol = length(unique(binom_df[,time])))) # time - RW2
  # A.mat.list[[5]] <- A_st
  
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
  # re_list$unstruct_time <- t.draws[unstruct_time_idx, ]
  # re_list$struct_time_1 <- t.draws[struct_time_idx1, ]
  # re_list$struct_time_2 <- t.draws[struct_time_idx2, ]
  # re_list$struct_time_3 <- t.draws[struct_time_idx3, ]
  re_list$struct_space <- t.draws[struct_space_both_idx, ]
  # re_list$spacetime <- t.draws[spacetime_idx, ]
  t.parnames <- c(names(SD0$par.fixed), names(SD0$par.random))
  # t.param.beta.idx <- which(t.parnames == "beta")
  # re_list$betas <- t.draws[t.param.beta.idx, ]
  re_list$cluster_effect <- t.draws[cluster_effect_idx, ]
  
  ## separate out the tmb_draws
  # NOTE: t.parnames will be the same as names(mu)
  t.param.idx <- grep("alpha", t.parnames)
  t.space.phi.idx <- grep("space_logit_phi", t.parnames)
  t.space.tau.idx <- grep("space_log_tau", t.parnames)
  # t.time.tau.idx <- which("time_log_tau" == t.parnames)
  # t.age.rw2.tau.idx <- grep("age_rw2_log_tau", t.parnames)
  # t.spacetime.tau.idx <- grep("spacetime_log_tau", t.parnames)
  # t.beta.tau.idx <- grep("beta_log_tau", t.parnames)
  t.cluster.tau.idx <- grep("cluster_log_tau", t.parnames)
  
  # make a list containing all draws from parameters
  param_list <- list()
  param_list$intercepts <- t.draws[t.param.idx, ]
  param_list$space_phi <- expit(t.draws[t.space.phi.idx,])
  param_list$space_tau <- exp(t.draws[t.space.tau.idx,])
  # param_list$time_unstruct_tau <- exp(t.draws[t.time.tau.idx,])
  # param_list$time_struct_tau <- exp(t.draws[t.age.rw2.tau.idx,])
  # param_list$spacetime_tau <- exp(t.draws[t.spacetime.tau.idx,])
  # param_list$beta_tau <- exp(t.draws[t.beta.tau.idx,])
  param_list$cluster_tau <- exp(t.draws[t.cluster.tau.idx,])
  
  # estimates of tmb field
  # below is the total bym2 term for space, not just the structured part
  t.space.total.idx <- head(grep("Epsilon_s", t.parnames), N)
  # t.time.unstruct.idx <- grep("Epsilon_t", t.parnames)
  # t.time.age1.idx <- grep("age_rw2_1", t.parnames)
  # t.time.age2.idx <- grep("age_rw2_2", t.parnames)
  # t.time.age3.idx <- grep("age_rw2_3", t.parnames)
  # t.spacetime.idx <- grep("Epsilon_st", t.parnames)
  # id's for each age intercept is t.param.idx
  
  # now combine space and time and age appropriately to get estimates through space and time in order
  # of arrange(region, time, age)
  linpred_ids <- expand.grid(region = t.space.total.idx, 
                             # time = t.time.unstruct.idx, 
                             age = t.param.idx) %>% 
    arrange(region, age)
  # add in ids for t.time.age1.idx, t.time.age2.idx, t.time.age3.idx
  # diff1 <- t.time.age1.idx[1] - t.time.unstruct.idx[1]
  # diff2 <- t.time.age2.idx[1] - t.time.unstruct.idx[1]
  # diff3 <- t.time.age3.idx[1] - t.time.unstruct.idx[1]
  # linpred_ids$age1 <- linpred_ids$time + diff1
  # linpred_ids$age2 <- linpred_ids$time + diff2
  # linpred_ids$age3 <- linpred_ids$time + diff3
  # # get spacetime column based on region and time variables
  # linpred_ids_st <- expand.grid(region = t.space.total.idx,
  #                               time = t.time.unstruct.idx) %>% 
  #   arrange(region, time) %>%
  #   mutate(spacetime = t.spacetime.idx)
  # suppressMessages(linpred_ids <- linpred_ids %>% left_join(linpred_ids_st))
  
  # combine draws for linear predictor to get space/time/age risks
  # this bit adds in the total for space, intercepts for each age group
  e.fitted <- t.draws[linpred_ids$age,] + 
    # t.draws[linpred_ids$time,] +
    t.draws[linpred_ids$region,] 
    # t.draws[linpred_ids$spacetime,]
  
  # add in age_group column to linpred_ids
  # age_to_agegroup_df <- binom_df[,c(age, age_group)] %>% distinct(.keep_all = TRUE)
  # colnames(age_to_agegroup_df) <- c("age","age_group")
  # suppressMessages(linpred_ids <- linpred_ids %>% left_join(age_to_agegroup_df))
  
  # make region_id and time_id variables that will help us with the random slopes in time
  # these will be 1:length(unique(region)) and 1:length(unique(time))
  # diff_reg <- (linpred_ids$region %>% min()) - 1
  # diff_time <- (linpred_ids$time %>% min()) - 1
  # linpred_ids$region_id <- linpred_ids$region - diff_reg
  # linpred_ids$time_id <- linpred_ids$time - diff_time
  
  # get centered and scaled version of time for the slopes
  # time_slope <- 1:S
  # center <- (N + 1)/2 + 1e-05
  # time_slope <- (time_slope - center)/(N - 1)
  
  # add in rw2 for time in each age group now, as well as random slopes in time for each area
  # for (i in 1:nrow(e.fitted)) {
  #   which_agegroup <- linpred_ids$age_group[i]
  #   which_reg <- linpred_ids$region_id[i]
  #   which_time <- linpred_ids$time_id[i]
  #   # add in appropriate age group
  #   if (which_agegroup == 1) {
  #     e.fitted[i,] <- e.fitted[i,] + t.draws[linpred_ids$age1[i],] +
  #       t.draws[t.param.beta.idx[which_reg],] * time_slope[which_time]
  #   } else if (which_agegroup == 2) {
  #     e.fitted[i,] <- e.fitted[i,] + t.draws[linpred_ids$age2[i],] +
  #       t.draws[t.param.beta.idx[which_reg],] * time_slope[which_time]
  #   } else {
  #     e.fitted[i,] <- e.fitted[i,] + t.draws[linpred_ids$age3[i],] +
  #       t.draws[t.param.beta.idx[which_reg],] * time_slope[which_time]
  #   }
  # }
  
  # lono correction
  h <- 16 * sqrt(3) / (15 * pi)
  cluster_vars <- 1/param_list$cluster_tau
  denom_samps <- sqrt(1 + h^2 * cluster_vars)
  
  # correct fitted samps
  for (i in 1:nrow(e.fitted)) {
    e.fitted[i,] <- e.fitted[i,] / denom_samps
  }
  
  # combine age groups using hazards
  
  # create fitted mat for space to put results in
  e.fitted.st <- matrix(0, nrow = N, ncol = nsamp)
  for (i in 1:N) {
      # get rows of e.fitted that correspond to this age/time
      row_ids <- (spaceage_id_df %>% filter(region == i))$spaceage_id
      # do hazards computation
      temp_samps <- e.fitted[row_ids,]
      for (m in 1:length(row_ids)) {
        #temp_samps[m,] <- (1 - expit(e.fitted[row_ids[m],]))^age_n[m]
        temp_samps[m,] <- (1 - expit(temp_samps[m,]))^age_n[m]
      }
      samps <- 1 - apply(temp_samps,2,prod)
      #samps <- 1 - apply(1 - expit(e.fitted[row_ids,]), 2, prod)
      # assign to appropriate row of e.fitted.st
      e.fitted.st[i,] <- samps
  }
  
  # matrix of samples of fitted values on probability scale
  fitted_mat <- e.fitted.st
  
  # return fitted mat, re_list, param_list
  return(list(fitted_mat = fitted_mat,
              sta_mat = e.fitted,
              re_list = re_list,
              param_list = param_list,
              runtime = end_time - begin_time))
}






