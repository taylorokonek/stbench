#' Fit U5MR model in TMB
#' 
#' Fit U5MR model in TMB. Optionally include age-specific intercepts, IID random effects in time,
#' BYM2 spatial random effects, area-specific random slopes in time, RWs in time for each age group,
#' and Type IV Knorr-Held interactions. Can fit benchmarked models, unbenchmarked models, or both.
#' Simultaneous benchmarking is performed via a second likelihood for national level estimates 
#' at each time point. Produces posterior samples for fitted values, hyperparameters, 
#' random effects, and fixed effects. All models use either a lono-binomial likelihood, described
#' in Dong and Wakefield (2021) or a betabinomial likelihood at the cluster level.
#' 
#' @param binom_df a dataframe containing binomial counts, and the following columns:
#' \itemize{
#'  \item \code{cluster}: cluster ID
#'  \item \code{y}: deaths
#'  \item \code{Ntrials}: total number of person-months in this \code{age}, \code{cluster}, and \code{region}
#'  \item \code{region}: area
#'  \item \code{time}: year
#'  \item \code{survey}: survey id. Only needs to be in binom_df if you want a survey 
#'  random effect in your model
#'  \item \code{age}: id for each age group (must be numeric) you want intercepts for
#'  \item \code{age_group}: id for each age group you want RW2s for 
#'  \item \code{age_n}: the number of months in each \code{age}
#' }
#' @param terms a character vector containing the names of terms to include in the linear predictor.
#' Options include:
#' \itemize{
#'  \item \code{"intercepts"}: age-specific intercepts. If this is not included in \code{terms},
#'  a model will be fit with a global intercept. A model with no intercepts is not currently implemented.
#'  \item \code{"iid_time"}: an iid random effect in time
#'  \item \code{"bym2_space"}: a BYM2 random effect in space. This is the only option for a spatial
#'  random effect at this time.
#'  \item \code{"area_slopes"}: area-specific random slopes in time
#'  \item \code{"rw2s_time"}: rw2s in time for each age group specified in \code{age_group}
#'  \item \code{"rw1s_time"}: To Do
#'  \item \code{"typei"}: To Do. A Type I Knorr-Held space-time interaction, with RW1 in time
#'  \item \code{"typeii"}: To Do. A Type II Knorr-Held space-time interaction, with RW1 in time
#'  \item \code{"typeiii"}: To Do. A Type III Knorr-Held space-time interaction, with RW1 in time
#'  \item \code{"typeiv"}: A Type IV Knorr-Held space-time interaction, with RW1 in time
#' }
#' Current default is \code{terms = c("intercepts", "iid_time", "bym2_space", "area_slopes", 
#' "rw2s_time", "typeiv")}.
#' @param cluster the column in \code{binom_df} corresponding to cluster id
#' @param y the column in \code{binom_df} corresponding to deaths
#' @param Ntrials the column in \code{binom_df} corresponding to total number of person-months
#' @param region the column in \code{binom_df} corresponding to region/area. Must be numeric and
#' start from 1.
#' @param time the column in \code{binom_df} corresponding to year. Must be numeric and 
#' start from 1.
#' @param survey the column in \code{binom_df} corresponding to survey id. Must be numeric and
#' start from 1.
#' @param years_predict a numeric sequence containing the years in which we want to make 
#' predictions. Must start at 1, and the max year must be greater than or equal to the 
#' maximum value in \code{time}. Defaults to 1 through the maximum value in \code{time}.
#' @param age the column in \code{binom_df} corresponding to age id, used for age-specific 
#' intercepts. Must be numeric and start from 1.
#' @param age_group the column in \code{binom_df} corresponding to age group id, used for
#' age-group-specific RWs in time. Must be numeric and start from 1.
#' @param age_n the number of months in each age group specified by the \code{age} column in 
#' \code{binom_df}. The length of \code{age_n} must be equal to the number of unique age groups
#' specified in \code{age}.
#' @param hiv_adj An optional log offset in time to include in the linear predictor. Defaults to 
#' no log offset. If included, the length of \code{hiv_adj} must be equal to the number of unique
#' time points specified in \code{years_predict} times the number of surveys present in your dataset, 
#' in order arrange(survey, time).
#' @param natl a vector of national level estimates, arranged in order of time
#' @param natl_sd a vector of standard deviations for national level estimates, arranged in order of
#' time
#' @param pop_weights a vector of population weights for use in the benchmarking constraint. Must
#' sum to one at each time point, and be in order arrange(region, time)
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
#' @export fit_u5mr
fit_u5mr <- function(binom_df, 
                     terms = c("intercepts","iid_time","bym2_space","area_slopes","rw2s_time","typeiv"),
                     cluster = "cluster",
                     y = "y", 
                     Ntrials = "N", 
                     region = "admin1", 
                     time = NULL,
                     survey = NULL,
                     years_predict = NULL,
                     age = "age_id",
                     age_n = c(1,11,12,12,12,12),
                     age_group = "age_group_id",
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
  if (!is.null(time)) {
    if (!(time %in% colnames(binom_df))) {
      stop("time must be a column in binom_df")
    }
  }
  if (!(age %in% colnames(binom_df))) {
    stop("age must be a column in binom_df")
  }
  if (!(age_group %in% colnames(binom_df))) {
    stop("age_group must be a column in binom_df")
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
  if (!is.null(survey)) {
    if (!(survey %in% colnames(binom_df))) {
      stop("survey must be a column in binom_df")
    }
  }
  
  # region, time, age, age_group, cluster must all be numeric and start from 1
  if (!(is.numeric(binom_df[,region]) & (min(binom_df[,region]) == 1))) {
    stop("region must be numeric and start from 1")
  }
  if (!is.null(time)) {
    if (!(is.numeric(binom_df[,time]) & (min(binom_df[,time]) == 1))) {
      stop("time must be numeric and start from 1")
    }
  }
  if (!(is.numeric(binom_df[,age]) & (min(binom_df[,age]) == 1))) {
    stop("age must be numeric and start from 1")
  }
  if (!(is.numeric(binom_df[,age_group]) & (min(binom_df[,age_group]) == 1))) {
    stop("age_group must be numeric and start from 1")
  }
  if (!(is.numeric(binom_df[,cluster]))) {
    stop("cluster must be numeric")
  }
  if (!(is.numeric(binom_df[,y]))) {
    stop("y must be numeric")
  }
  if (!(is.numeric(binom_df[,Ntrials]))) {
    stop("Ntrials must be numeric")
  }
  if (!is.null(survey)) {
    if (!(is.numeric(binom_df[,survey]) & (min(binom_df[,survey]) == 1))) {
      stop("survey must be numeric and start from 1")
    }
  }
  
  # relabel clusters from 1:length(unique(cluster))
  binom_df$cluster <- binom_df[,cluster]
  suppressMessages(binom_df <- data.frame(cluster = sort(unique(binom_df[,cluster])), cluster_new = 1:length(unique(binom_df[,cluster]))) %>%
                     right_join(binom_df))
  cluster <- "cluster_new"
  
  # if years_predict is NULL, make it equal to the min in binom_df[,time] to the max in binom_df[,time]
  if (!is.null(years_predict)) {
    # check that year_label starts at 1
    if (min(years_predict) != 1) {
      stop("years_predict must start at 1")
    }
    if (max(years_predict) < max(binom_df[,time])) {
      stop("maximum value in years_predict must be greater than or equal to maximum value in time")
    }
    # addtl error check that year_label is a sequence from its min to max value - TO DO
    
  } else {
    if (!is.null(time)) {
      min_year <- min(binom_df[,time])
      max_year <- max(binom_df[,time])
      years_predict = min_year:max_year
    }
  }
  
  # create (scaled) Q_struct_space, Q_struct_time, Q_struct_time_interaction
  
  if ("bym2_space" %in% terms) {
    if (is.null(Q_struct_space)) {
      stop("If 'bym2_space' is included in terms, Q_struct_space must be specified")
    } else {
      Q_struct_space <- INLA::inla.scale.model(Q_struct_space, 
                                               constr = list(A = matrix(1, nrow = 1, ncol = nrow(Q_struct_space)), e = 0))
    }
  }
  
  if ("rw2s_time" %in% terms) {
    if ("rw1s_time" %in% terms) {
      stop("Must specify only one of 'rw1s_time', 'rw2s_time' in terms")
    } else {
      Q_struct_time <- INLA:::inla.rw(n = length(years_predict), order = 2, scale.model = TRUE, sparse = TRUE)
    }
  }
  
  if ("rw1s_time" %in% terms) {
    Q_struct_time <- INLA:::inla.rw(n = length(years_predict), order = 1, scale.model = TRUE, sparse = TRUE)
  }
  
  if ("typeiv" %in% terms) {
    Q_struct_time_interaction <- INLA:::inla.rw(n = length(years_predict), order = 1, scale.model = TRUE, sparse = TRUE)
  }
  if (any(c("typei","typeii","typeiii") %in% terms)) {
    stop("typei, typeii, typeiii are not currently implemented terms")
  }
  
  # determine which model to fit
  # space-only
  if (identical(sort(terms), c("bym2_space","intercepts"))) {
    if (is.null(survey)) {
      if (family == "binomial") {
        if (benched == "unbenched") {
          message("Unbenched Binomial model - single survey, space only")
          unbenched_res <- tmb_u5mr_intercepts_bym2(binom_df = binom_df,
                                                    cluster = cluster,
                                                    y = y,
                                                    Ntrials = Ntrials,
                                                    region = region,
                                                    age = age,
                                                    age_n = age_n,
                                                    age_group = age_group,
                                                    hiv_adj = hiv_adj,
                                                    Q_struct_space = Q_struct_space,
                                                    nsamp = nsamp)
        } else if (benched == "benched") {
          message("Benched Binomial model - single survey, space only")
          benched_res <- tmb_u5mr_intercepts_bym2_benched(binom_df = binom_df,
                                                          cluster = cluster,
                                                          y = y,
                                                          Ntrials = Ntrials,
                                                          region = region,
                                                          age = age,
                                                          age_n = age_n,
                                                          age_group = age_group,
                                                          hiv_adj = hiv_adj,
                                                          natl = natl,
                                                          natl_sd = natl_sd,
                                                          pop_weights = pop_weights,
                                                          Q_struct_space = Q_struct_space,
                                                          nsamp = nsamp)
        } else {
          message("Unbenched Binomial model - single survey, space only")
          unbenched_res <- tmb_u5mr_intercepts_bym2(binom_df = binom_df,
                                                    cluster = cluster,
                                                    y = y,
                                                    Ntrials = Ntrials,
                                                    region = region,
                                                    age = age,
                                                    age_n = age_n,
                                                    age_group = age_group,
                                                    hiv_adj = hiv_adj,
                                                    Q_struct_space = Q_struct_space,
                                                    nsamp = nsamp)
          message("Benched Binomial model - single survey, space only")
          benched_res <- tmb_u5mr_intercepts_bym2_benched(binom_df = binom_df,
                                                          cluster = cluster,
                                                          y = y,
                                                          Ntrials = Ntrials,
                                                          region = region,
                                                          age = age,
                                                          age_n = age_n,
                                                          age_group = age_group,
                                                          hiv_adj = hiv_adj,
                                                          natl = natl,
                                                          natl_sd = natl_sd,
                                                          pop_weights = pop_weights,
                                                          Q_struct_space = Q_struct_space,
                                                          nsamp = nsamp)
        }
      } else {
        stop("betabinomial space-only single survey model not currently implemented")
      }
    } else {
      stop("space-only multi-survey models not currently implemented")
    }
    
    # time-only model
    
  } else if (identical(sort(terms),c("iid_time", "intercepts", "rw2s_time"))){
    if (is.null(survey)) {
      if (family == "binomial") {
        if (benched == "unbenched") {
          message("Unbenched Binomial model - single survey, time only")
          unbenched_res <- tmb_u5mr_intercepts_iidtime_rw2s(binom_df = binom_df,
                                                            cluster = cluster,
                                                            y = y,
                                                            Ntrials = Ntrials,
                                                            time = time,
                                                            age = age,
                                                            age_n = age_n,
                                                            age_group = age_group,
                                                            hiv_adj = hiv_adj,
                                                            Q_struct_time = Q_struct_time,
                                                            nsamp = nsamp)
        } else if (benched == "benched") {
          message("Benched Binomial model - single survey, time only")
          benched_res <- tmb_u5mr_intercepts_iidtime_rw2s_benched(binom_df = binom_df,
                                                                  cluster = cluster,
                                                                  y = y,
                                                                  Ntrials = Ntrials,
                                                                  time = time,
                                                                  age = age,
                                                                  age_n = age_n,
                                                                  age_group = age_group,
                                                                  hiv_adj = hiv_adj,
                                                                  natl = natl,
                                                                  natl_sd = natl_sd,
                                                                  Q_struct_time = Q_struct_time,
                                                                  nsamp = nsamp)
        } else {
          message("Unbenched Binomial model - single survey, time only")
          unbenched_res <- tmb_u5mr_intercepts_iidtime_rw2s(binom_df = binom_df,
                                                            cluster = cluster,
                                                            y = y,
                                                            Ntrials = Ntrials,
                                                            time = time,
                                                            age = age,
                                                            age_n = age_n,
                                                            age_group = age_group,
                                                            hiv_adj = hiv_adj,
                                                            Q_struct_time = Q_struct_time,
                                                            nsamp = nsamp)
          message("Benched Binomial model - single survey, time only")
          benched_res <- tmb_u5mr_intercepts_iidtime_rw2s_benched(binom_df = binom_df,
                                                                  cluster = cluster,
                                                                  y = y,
                                                                  Ntrials = Ntrials,
                                                                  time = time,
                                                                  age = age,
                                                                  age_n = age_n,
                                                                  age_group = age_group,
                                                                  hiv_adj = hiv_adj,
                                                                  natl = natl,
                                                                  natl_sd = natl_sd,
                                                                  Q_struct_time = Q_struct_time,
                                                                  nsamp = nsamp)
        }
      } else {
        stop("betabinomial time-only single survey model not currently implemented")
      }
    } else {
      stop("time-only multi-survey models not currently implemented")
    }
    
    # full space-time model
  } else if (identical(sort(terms),c("area_slopes","bym2_space","iid_time","intercepts","rw2s_time","typeiv"))) {
    if (is.null(survey)) {
      if (family == "binomial") {
        if (benched == "unbenched") {
          message("Unbenched Binomial model - single survey")
          unbenched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV(binom_df = binom_df,
                                                                               cluster = cluster,
                                                                               y = y,
                                                                               Ntrials = Ntrials,
                                                                               region = region,
                                                                               time = time,
                                                                               age = age,
                                                                               age_n = age_n,
                                                                               age_group = age_group,
                                                                               hiv_adj = hiv_adj,
                                                                               Q_struct_space = Q_struct_space,
                                                                               Q_struct_time = Q_struct_time,
                                                                               Q_struct_time_interaction = Q_struct_time_interaction,
                                                                               nsamp = nsamp)
        } else if (benched == "benched") {
          message("Benched Binomial model - single survey")
          benched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV_benched(binom_df = binom_df,
                                                                                     cluster = cluster,
                                                                                     y = y,
                                                                                     Ntrials = Ntrials,
                                                                                     region = region,
                                                                                     time = time,
                                                                                     age = age,
                                                                                     age_n = age_n,
                                                                                     age_group = age_group,
                                                                                     hiv_adj = hiv_adj,
                                                                                     natl = natl,
                                                                                     natl_sd = natl_sd,
                                                                                     pop_weights = pop_weights,
                                                                                     Q_struct_space = Q_struct_space,
                                                                                     Q_struct_time = Q_struct_time,
                                                                                     Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                     nsamp = nsamp)
        } else {
          message("Unenched Binomial model - single survey")
          unbenched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV(binom_df = binom_df,
                                                                               cluster = cluster,
                                                                               y = y,
                                                                               Ntrials = Ntrials,
                                                                               region = region,
                                                                               time = time,
                                                                               age = age,
                                                                               age_n = age_n,
                                                                               age_group = age_group,
                                                                               hiv_adj = hiv_adj,
                                                                               Q_struct_space = Q_struct_space,
                                                                               Q_struct_time = Q_struct_time,
                                                                               Q_struct_time_interaction = Q_struct_time_interaction,
                                                                               nsamp = nsamp)
          
          message("Benched Binomial model - single survey")
          benched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV_benched(binom_df = binom_df,
                                                                                     cluster = cluster,
                                                                                     y = y,
                                                                                     Ntrials = Ntrials,
                                                                                     region = region,
                                                                                     time = time,
                                                                                     age = age,
                                                                                     age_n = age_n,
                                                                                     age_group = age_group,
                                                                                     hiv_adj = hiv_adj,
                                                                                     natl = natl,
                                                                                     natl_sd = natl_sd,
                                                                                     pop_weights = pop_weights,
                                                                                     Q_struct_space = Q_struct_space,
                                                                                     Q_struct_time = Q_struct_time,
                                                                                     Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                     nsamp = nsamp)
          
        }
      } else {
        if (benched == "unbenched") {
          message("Unbenched Betabinomial model - single survey")
          unbenched_res <- tmb_u5mr_betabinom_intercepts_iidtime_bym2_slopes_rw2s_typeIV(binom_df = binom_df,
                                                                                         cluster = cluster,
                                                                                         y = y,
                                                                                         Ntrials = Ntrials,
                                                                                         region = region,
                                                                                         time = time,
                                                                                         age = age,
                                                                                         age_n = age_n,
                                                                                         age_group = age_group,
                                                                                         hiv_adj = hiv_adj,
                                                                                         Q_struct_space = Q_struct_space,
                                                                                         Q_struct_time = Q_struct_time,
                                                                                         Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                         nsamp = nsamp)
        } else if (benched == "benched") {
          message("Benched Betabinomial model - single survey")
          benched_res <- tmb_u5mr_betabinom_intercepts_iidtime_bym2_slopes_rw2s_typeIV_benched(binom_df = binom_df,
                                                                                               cluster = cluster,
                                                                                               y = y,
                                                                                               Ntrials = Ntrials,
                                                                                               region = region,
                                                                                               time = time,
                                                                                               age = age,
                                                                                               age_n = age_n,
                                                                                               age_group = age_group,
                                                                                               hiv_adj = hiv_adj,
                                                                                               natl = natl,
                                                                                               natl_sd = natl_sd,
                                                                                               pop_weights = pop_weights,
                                                                                               Q_struct_space = Q_struct_space,
                                                                                               Q_struct_time = Q_struct_time,
                                                                                               Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                               nsamp = nsamp)
        } else {
          message("Unbenched Betabinomial model - single survey")
          unbenched_res <- tmb_u5mr_betabinom_intercepts_iidtime_bym2_slopes_rw2s_typeIV(binom_df = binom_df,
                                                                                         cluster = cluster,
                                                                                         y = y,
                                                                                         Ntrials = Ntrials,
                                                                                         region = region,
                                                                                         time = time,
                                                                                         age = age,
                                                                                         age_n = age_n,
                                                                                         age_group = age_group,
                                                                                         hiv_adj = hiv_adj,
                                                                                         Q_struct_space = Q_struct_space,
                                                                                         Q_struct_time = Q_struct_time,
                                                                                         Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                         nsamp = nsamp)
          
          message("Benched Betabinomial model - single survey")
          benched_res <- tmb_u5mr_betabinom_intercepts_iidtime_bym2_slopes_rw2s_typeIV_benched(binom_df = binom_df,
                                                                                               cluster = cluster,
                                                                                               y = y,
                                                                                               Ntrials = Ntrials,
                                                                                               region = region,
                                                                                               time = time,
                                                                                               age = age,
                                                                                               age_n = age_n,
                                                                                               age_group = age_group,
                                                                                               hiv_adj = hiv_adj,
                                                                                               natl = natl,
                                                                                               natl_sd = natl_sd,
                                                                                               pop_weights = pop_weights,
                                                                                               Q_struct_space = Q_struct_space,
                                                                                               Q_struct_time = Q_struct_time,
                                                                                               Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                               nsamp = nsamp)
          
        }
      }
    } else {
      if (family == "binomial") {
        if (benched == "unbenched") {
          message("Unbenched Binomial model - multi-survey")
          unbenched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV_multisurvey(binom_df = binom_df,
                                                                                           cluster = cluster,
                                                                                           y = y,
                                                                                           Ntrials = Ntrials,
                                                                                           region = region,
                                                                                           time = time,
                                                                                           survey = survey,
                                                                                           age = age,
                                                                                           age_n = age_n,
                                                                                           age_group = age_group,
                                                                                           hiv_adj = hiv_adj,
                                                                                           Q_struct_space = Q_struct_space,
                                                                                           Q_struct_time = Q_struct_time,
                                                                                           Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                           nsamp = nsamp)
        } else if (benched == "benched") {
          message("Benched Binomial model - multi-survey")
          benched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV_multisurvey_benched(binom_df = binom_df,
                                                                                                 cluster = cluster,
                                                                                                 y = y,
                                                                                                 Ntrials = Ntrials,
                                                                                                 region = region,
                                                                                                 time = time,
                                                                                                 survey = survey,
                                                                                                 age = age,
                                                                                                 age_n = age_n,
                                                                                                 age_group = age_group,
                                                                                                 hiv_adj = hiv_adj,
                                                                                                 natl = natl,
                                                                                                 natl_sd = natl_sd,
                                                                                                 pop_weights = pop_weights,
                                                                                                 Q_struct_space = Q_struct_space,
                                                                                                 Q_struct_time = Q_struct_time,
                                                                                                 Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                                 nsamp = nsamp)
        } else {
          message("Unenched Binomial model - multi-survey")
          unbenched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV_multisurvey(binom_df = binom_df,
                                                                                           cluster = cluster,
                                                                                           y = y,
                                                                                           Ntrials = Ntrials,
                                                                                           region = region,
                                                                                           time = time,
                                                                                           survey = survey,
                                                                                           age = age,
                                                                                           age_n = age_n,
                                                                                           age_group = age_group,
                                                                                           hiv_adj = hiv_adj,
                                                                                           Q_struct_space = Q_struct_space,
                                                                                           Q_struct_time = Q_struct_time,
                                                                                           Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                           nsamp = nsamp)
          
          message("Benched Binomial model - multi-survey")
          benched_res <- tmb_u5mr_intercepts_iidtime_bym2_slopes_rw2s_typeIV_multisurvey_benched(binom_df = binom_df,
                                                                                                 cluster = cluster,
                                                                                                 y = y,
                                                                                                 Ntrials = Ntrials,
                                                                                                 region = region,
                                                                                                 time = time,
                                                                                                 survey = survey,
                                                                                                 age = age,
                                                                                                 age_n = age_n,
                                                                                                 age_group = age_group,
                                                                                                 hiv_adj = hiv_adj,
                                                                                                 natl = natl,
                                                                                                 natl_sd = natl_sd,
                                                                                                 pop_weights = pop_weights,
                                                                                                 Q_struct_space = Q_struct_space,
                                                                                                 Q_struct_time = Q_struct_time,
                                                                                                 Q_struct_time_interaction = Q_struct_time_interaction,
                                                                                                 nsamp = nsamp)
          
        }
      } else {
        if (benched == "unbenched") {
          message("Unbenched Betabinomial model - multi-survey")
          stop("not currently implemented")
        } else if (benched == "benched") {
          message("Benched Betabinomial model - multi-survey")
          stop("not currently implemented")
        } else {
          message("Unbenched Betabinomial model - multi-survey")
          stop("not currently implemented")
          message("Benched Betabinomial model - multi-survey")
        }
      }
    }
  } else {
    stop("model with terms specified not currently implemented")
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





