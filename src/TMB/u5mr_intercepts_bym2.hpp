// U5MR model: 
//     Age-specific intercepts, 
//     BYM2 for space

#ifndef u5mr_intercepts_bym2_hpp
#define u5mr_intercepts_bym2_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// include libraries
// #include <TMB.hpp>
#include "stbench/addtl_densities.hpp"
#include "stbench/helper_funcs.hpp"
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;




///////////////////////////
// our main function     //
// to calculate the jnll //
///////////////////////////
template<class Type>
Type u5mr_intercepts_bym2(objective_function<Type>* obj)
{

  // Data and parameters passed in from R

  // Total number of space/time/age/survey/cluster 
  DATA_INTEGER( N );   // Number of spatial regions
  // DATA_INTEGER( S );   // Number of time points
  DATA_INTEGER( A_m ); // Number of age groups
  // DATA_INTEGER( K );   // Number of surveys
  DATA_INTEGER( nclusts ); // Number of clusters 
  DATA_INTEGER( nobs ); // Number of observations 

  // Data
  DATA_VECTOR( y_iac );   // Number of deaths in space/time/age/cluster
  DATA_VECTOR( N_iac );   // Ntrials for binomial in space/time/age/cluster
  DATA_VECTOR( hiv_adj );  // HIV adjustment ratios - log offsets for each time point. If no offset wanted, hiv_adj can be a vector of 1's
  
  // Ids for space/time/age/survey
  DATA_IVECTOR( region_id ); // Numeric indicator for which region an observation belongs to
  // DATA_IVECTOR( time_id ); // Numeric indicator for which time point an observation belongs to
  DATA_IVECTOR( age_id ); // Numeric indicator for which age group an observation belongs to
  // DATA_IVECTOR( spacetimeage_id ); // Numeric indicator for which space/time/age point an observation belongs to
  // DATA_IVECTOR( spacetime_id ); // Numeric indicator for which space/time point an observation belongs to
  DATA_IVECTOR( spaceage_id ); // Numeric indicator for which space/age point an observation belongs to
  // DATA_IVECTOR( age_group_id ); // Numeric indicator for which age group an observation belongs to
  DATA_IVECTOR( cluster_id );
  // DATA_IVECTOR( survey_id ); // Numeric indivator for which survey an observation belongs to

  // Benchmarking things
  // DATA_VECTOR( natl ); // National-level estimates we want to benchmark to, in order of time
  // DATA_VECTOR( natl_sd ); // SDs for national-level estimates, in order of time
  // DATA_VECTOR( pop_weights ); // population weights for benchmarking constraint. Must sum to 1 at each time point and be in order arrange(region, time)
  // DATA_IVECTOR( age_n ); // number of agemonths in each age group, needed for getting risk aggregated to space/time for benchmarking constraint
  
  // Precision matrices for space (structured), time (unstructured and structured), survey
  DATA_SPARSE_MATRIX( Q_struc_space ); // structured precision matrices must be scaled
  // DATA_SPARSE_MATRIX( Q_struc_time ); // structured precision matrices must be scaled
  // DATA_SPARSE_MATRIX( Q_struc_spacetime ); // structured precision matrices must be scaled
  // DATA_SPARSE_MATRIX( Q_unstruc_survey ); // should be an identity matrix (K x K) MULTIPLIED BY a fixed precision parameter

  // Constraint matrices
  // DATA_MATRIX( A_st ); // constraint matrix for the space-time-interaction term

  // Prior specifications - currently unused. 
  // *_tau_pri = c(u, prec), prior on log(tau) with sigma ~ N(u, 1/sqrt(prec))
  // DATA_VECTOR( bym2_alpha_pri ); // bym2_alpha_pri = c(u, prec) with alpha ~ Norm(u, 1/sqrt(prec)). same prior for each age group
  // DATA_VECTOR( bym2_space_phi_pri );   // bym2_space_phi_pri = c(a, b), prior on logit(phi) with phi~Beta(a,b)
  // DATA_VECTOR( bym2_space_tau_pri );   
  // DATA_VECTOR( time_tau_pri );   
  // DATA_VECTOR( age_rw2_tau_pri );  
  // DATA_VECTOR( spacetime_tau_pri ); 
  // DATA_VECTOR( d_pri );          // overdispersion parameter if we're fitting a betabinomial model
  // DATA_VECTOR( beta_pri ); //specifying u and prec for beta_i                     with beta_i~Norm(u,prec)

  // Fixed effects & hyperpars
  PARAMETER_ARRAY( alpha );           // Intercept. An array because there are different intercepts for each age group, a A_m x 1 array
  PARAMETER( space_logit_phi );          // logit of BYM2 mixing parameter
  PARAMETER( space_log_tau );            // Log of tau param (precision of space covariance matrix)
  // PARAMETER( time_log_tau );            // Log of tau param (precision of temporal covariance matrix - the iid bit)
  // PARAMETER( age_rw2_log_tau );          // Log of tau param (precision on temporal covariance matrix - the rw2s for age groups)
  // PARAMETER( spacetime_log_tau );        // Log of tau param (precision on spatio-temporal covariance matrix - the type IV Knorr-Held interaction)
  // PARAMETER_ARRAY( survey_coefs );     //  coefficients for the survey fixed effect K x 1
  // PARAMETER( beta_log_tau );
  PARAMETER( cluster_log_tau );          // Log of tau param for cluster random effect
  // PARAMETER( logit_d ); // overdispersion parameter for betabinomial distribution

  // Random effects
  PARAMETER_ARRAY( Epsilon_s );    // effects for each spatial region (length 2N). first half is total, second half is structured
  // PARAMETER_ARRAY( Epsilon_t );    // effects for each time point.
  // PARAMETER_ARRAY( Epsilon_st );   // effects for space-time interaction (length S*N)
  // PARAMETER_ARRAY( age_rw2_1 );    // rw2s for 3 different age groups
  // PARAMETER_ARRAY( age_rw2_2 );
  // PARAMETER_ARRAY( age_rw2_3 );
  // PARAMETER_ARRAY( beta );     // random slopes in time for each area, beta has dimension N x 1
  PARAMETER_ARRAY( cluster_effect );



  /////////////////////////////////
  // Define internal objects     //
  /////////////////////////////////


  // Compute time_slope, which is a centered and scaled version of 1:S, to use for our random slopes in time for each area
  // vector<Type> time_slope(S);
  // Type center = (N + 1)/2 + Type(0.00001);
  // for (int i = 0; i < S; i++) {
  //   time_slope(i) = (i + 1 - center) / (N - 1);
  // }

  // Multiply certain precision matrices by their precision parameters
  // Not needed for anything BYM2 since that is done in the bym2_Q() function
  // Not needed for Q_unstruc_survey since we fix the precision here and have it multiplied in already in R (prior to tmb)
  // We'll do this using the SCALE() function when adding in contributions to the jnll

  // Q_struc_time
  // Type age_rw2_tau = exp(age_rw2_log_tau);

  // Q_unstruc_time
  // Type time_sd = 1/sqrt(exp(time_log_tau));

  // Q_struc_spacetime
  // Type spacetime_tau = exp(spacetime_log_tau);

  // Q_struc_space
  // Need to create BYM2 precision matrix
  // Make spatial precision matrix
  SparseMatrix<Type> Q_space = bym2_Q(space_logit_phi, space_log_tau, Q_struc_space);
  int full_N = N + N;      // top half is total effect, bottom half is structured effect

  // Set up objective function -- joint negative log-likelihood
  // Define it in three parts (just for ease of reading code): (1) priors (2) space/time/age fields (3) data contribution
  vector<Type> jnll_comp(3);
  Type jnll = 0; // sum of jnll_comp

  // Initialize vectors for risk (space/time/age and space/time)
  vector<Type> risk_ia( N * A_m );          // expit(total latent field) risk for each region/age group point, arranged in order arrange(region, age)
  // vector<Type> risk_i( N );   // expit(total latent field) risk for each region

  // Initialize vectors for survey, space (also total and structured components), spacetime, rw2s by age group...
  // vector<Type> survey_coefs_vec( K ); // survey fixed effect
  vector<Type> epsilon_s( 2*N );     // joint total (followed by) structured effects vec
  vector<Type> total_latent_i( N );  // struc + unstruc in space
  vector<Type> struc_latent_i( N );  // struc in space
  // vector<Type> epsilon_st( N * S ); // spacetime
  // vector<Type> age_rw2_1_vec( S ); // rw2s by age group
  // vector<Type> age_rw2_2_vec( S );
  // vector<Type> age_rw2_3_vec( S );

  // Fill in vectors for survey, space, spacetime, rw2s by age group...
  // for(int s = 0; s < K; s++) {
  //   survey_coefs_vec[s] = survey_coefs(s);
  // }
  for(int s = 0; s < N; s++){
    epsilon_s[s] = Epsilon_s(s);
    epsilon_s[N + s] = Epsilon_s(N + s);
  }
  // for(int s = 0; s < (N * S); s++) {
  //   epsilon_st[s] = Epsilon_st(s);
  // }
  // for(int s = 0; s < S; s++){
  //   age_rw2_1_vec[s] = age_rw2_1(s);
  //   age_rw2_2_vec[s] = age_rw2_2(s);
  //   age_rw2_3_vec[s] = age_rw2_3(s);
  // }



  //////////////////////////////////////////////////////
  // Constrain things that need to be constrained     //
  //////////////////////////////////////////////////////



  // a constrained vector x_c for an unconstrained vector x is calculated as:
  // x_c = x - Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e)
  // for constraint Ax=e, and for GMRF x with precision Q
  // we'll build up the matrices we need piece by piece

  // Invert precision matrices (possible since we added a small value to the diagonal)
  matrix<Type> Q_inv_space = invertSparseMatrix(Q_space);
  // matrix<Type> Q_inv_st = invertSparseMatrix(Q_struc_spacetime);
  // matrix<Type> Q_inv_rw2 = invertSparseMatrix(Q_struc_time);
  // matrix<Type> Q_inv_survey = invertSparseMatrix(Q_unstruc_survey);

  // Create constraint matrices A
  matrix<Type> A_space(1, full_N); // sum all components
  for(int i = 0; i < N; i++){
    A_space(0, i)     = 0; // no constraints on total effects
    A_space(0, N + i) = 1; // sum-to-0 on structured effects
  }
  // matrix<Type> A_rw2(1, S);
  // for(int i = 0; i < S; i++) {
  //   A_rw2(0, i) = 1; // sum-to-0 constraint
  //   // A_rw2(1, i) = i + 1; // constraint on slopes (1, 2, ..., S)
  // }
  // matrix<Type> A_survey(1, K); // survey fixed effect - sum to zero constraint
  // for(int i = 0; i < K; i++) {
  //   A_survey(0, i) = 1; 
  // }

  // Create A^T
  matrix<Type> A_space_T = A_space.transpose(); // bym2 in space
  // matrix<Type> A_st_T = A_st.transpose(); // A_st input as data from R
  // matrix<Type> A_rw2_T = A_rw2.transpose();
  // matrix<Type> A_survey_T = A_survey.transpose();

  // Create Q^{-1}A^T
  matrix<Type> QinvA_space_mat = Q_inv_space * A_space_T;
  // matrix<Type> QinvA_st = Q_inv_st * A_st_T;
  // matrix<Type> QinvA_rw2= Q_inv_rw2 * A_rw2_T;
  // matrix<Type> QinvA_survey_mat = Q_inv_survey * A_survey_T;

  // For precision matrices with sum-to-zero constraints (survey, space), convert matrix to vector
  vector<Type> QinvA_space(full_N);
  for(int i = 0; i < full_N; i++) {
    QinvA_space(i) = QinvA_space_mat(i,0);
  }
  // vector<Type> QinvA_survey(K);
  // for(int i=0;i<K;i++){
  //   QinvA_survey(i) = QinvA_survey_mat(i,0);
  // }

  // Create AQ^{-1}A^T
  matrix<Type> AQinvA_space = A_space * QinvA_space_mat;
  // matrix<Type> AQinvA_rw2 = A_rw2 * QinvA_rw2;
  // matrix<Type> AQinvA_st = A_st * QinvA_st;
  // matrix<Type> AQinvA_survey = A_survey * QinvA_survey_mat;

  // Create (AQ^{-1}A^T)^{-1}
  matrix<Type> AQinvA_space_inv = AQinvA_space.inverse(); // okay for small matrices
  // matrix<Type> AQinvA_rw2_inv = AQinvA_rw2.inverse(); // okay for small matrices
  // matrix<Type> AQinvA_st_inv = atomic::matinv(AQinvA_st); // better for large matrices
  // matrix<Type> matrix<Type> AQinvA_survey_inv = AQinvA_survey.inverse(); // okay for small matrices

  // Create Ax
  matrix<Type> Ax_space = (A_space * epsilon_s.matrix());
  // matrix<Type> Ax_rw2_1 = (A_rw2 * age_rw2_1_vec.matrix());
  // matrix<Type> Ax_rw2_2 = (A_rw2 * age_rw2_2_vec.matrix());
  // matrix<Type> Ax_rw2_3 = (A_rw2 * age_rw2_3_vec.matrix());
  // matrix<Type> Ax_st = (A_st * epsilon_st.matrix());
  // matrix<Type> Ax_survey = (A_survey * survey_coefs_vec.matrix());

  // Convert Ax from matrix to vector form - needed for dnorm & MVNORM, so only needed for survey fixed effect, rw2s by age group, and space/time interaction
  vector<Type> Ax_space_vec(1);
  Ax_space_vec(0) = Ax_space(0, 0);
  // vector<Type> Ax_rw2_1_vec(1);
  // vector<Type> Ax_rw2_2_vec(1);
  // vector<Type> Ax_rw2_3_vec(1);
  // for(int i = 0; i < 1; i++) {
  //   Ax_rw2_1_vec(i) = Ax_rw2_1(i,0);
  //   Ax_rw2_2_vec(i) = Ax_rw2_2(i,0);
  //   Ax_rw2_3_vec(i) = Ax_rw2_3(i,0);
  // }
  // vector<Type> Ax_st_vec( N + S - 1 );
  // for(int i = 0; i < (N + S - 1); i++) {
  //   Ax_st_vec(i) = Ax_st(i,0);
  // }
  // vector<Type> Ax_survey_vec(1);
  // Ax_survey_vec(0) = Ax_survey(0, 0);
 
  // Convert Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e) to vector form for conditioning by kriging correction
  matrix<Type> krig_correct_space = QinvA_space_mat * AQinvA_space_inv * Ax_space;
  // matrix<Type> krig_correct_1 = QinvA_rw2 * AQinvA_rw2_inv * Ax_rw2_1;
  // matrix<Type> krig_correct_2 = QinvA_rw2 * AQinvA_rw2_inv * Ax_rw2_2;
  // matrix<Type> krig_correct_3 = QinvA_rw2 * AQinvA_rw2_inv * Ax_rw2_3;
  // matrix<Type> krig_correct_st = QinvA_st * AQinvA_st_inv * Ax_st;
  // matrix<Type> krig_correct_survey = QinvA_survey_mat * AQinvA_survey_inv * Ax_survey;
  vector<Type> krig_correct_space_vec(full_N);
  // vector<Type> krig_correct_1_vec(S);
  // vector<Type> krig_correct_2_vec(S);
  // vector<Type> krig_correct_3_vec(S);
  // vector<Type> krig_correct_st_vec(N * S);
  // vector<Type> krig_correct_survey_vec(K);
  for (int i = 0; i < full_N; i++) {
    krig_correct_space_vec(i) = krig_correct_space(i,0);
  }
  // for (int i = 0; i < S; i++) {
  //   krig_correct_1_vec(i) = krig_correct_1(i,0);
  //   krig_correct_2_vec(i) = krig_correct_2(i,0);
  //   krig_correct_3_vec(i) = krig_correct_3(i,0);
  // }
  // for (int i = 0; i < (N * S); i++) {
  //   krig_correct_st_vec(i) = krig_correct_st(i,0);
  // }
  // for (int i = 0; i < K; i++) {
  //   krig_correct_survey_vec(i) = krig_correct_survey(i,0);
  // }

  // Construct constrained vector x_c = x - Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e)
  vector<Type> epsilon_s_c = epsilon_s - krig_correct_space_vec;
  // vector<Type> age_rw2_1_c = age_rw2_1_vec - krig_correct_1_vec;
  // vector<Type> age_rw2_2_c = age_rw2_2_vec - krig_correct_2_vec;
  // vector<Type> age_rw2_3_c = age_rw2_3_vec - krig_correct_3_vec;
  // vector<Type> epsilon_st_c = epsilon_st - krig_correct_st_vec;
  // vector<Type> survey_coefs_c = survey_coefs_vec - krig_correct_survey_vec;

  // For BYM2s (only space at the moment), fill in total and struc vectors for use in risk calculation
  for(int s = 0; s < N; s++){
    total_latent_i[s] = epsilon_s_c(s);     // first half is total
    struc_latent_i[s] = epsilon_s_c(N + s); // second half is structured effect
  }

  // Multiply precision matrices by their precision parameters
  // for (int i = 0; i < S; i++) {
  //   for (int j = 0; j < S; j++) {
  //     Q_struc_time.coeffRef(i, j) = Q_struc_time.coeffRef(i, j) * age_rw2_tau;
  //   }
  // }
  // for (int i = 0; i < (N*S); i++) {
  //   for (int j = 0; j < (N*S); j++) {
  //     Q_struc_spacetime.coeffRef(i, j) = Q_struc_spacetime.coeffRef(i, j) * spacetime_tau;
  //   }
  // }
 


  //////////////////////////////////////////////////////
  // Add to the joint negative log-likelihood         //
  // (1) space/time/age/survey fields                 //
  // (2) priors                                       //
  // (3) data                                         //
  //////////////////////////////////////////////////////



  /////////
  // (1) //
  /////////



  // NOTE: likelihoods from namespace 'density' already return NEGATIVE log-liks so we add
  //       other likelihoods return positive log-liks so we subtract

  // NOTE on constrained densities:
  // With constrained precisions, the density is:
  // pi(x | Ax) = pi(Ax | x) * pi(x) / pi(Ax)
  // where
  // pi(x)  is the unconstrained GMRF with precision Q_s
  // pi(Ax | x) is a gaussian with mean A*0 = 0, and covariance AQinvA^T
  // pi(Ax) is 0 if the constraint is not met, or the constant |AA^T|^{-.5} if it is metabolic
  //        with A = (0, 0, ..., 0, 1, 1, ..., 1), |AA^T| = 1, and log(1) = 0, so this contributes nothing for sum-to-zero constraints
  // HOWEVER. pi(Ax) does contribute something when A contains the two identifiability constraints for the RW2. We should not need
  // to add to the negative log likelihood log(det(AA^T)^(-1/2)) since this shouldn't affect optimization

  // Add in pi(Ax | x)
  // jnll_comp[0] -= dnorm(Ax_space_vec(0), Type(0.0), Type(sqrt(AQinvA_space(0,0))), true); // N(mean, sd)
  // jnll_comp[0] += MVNORM(AQinvA_rw2)(Ax_rw2_1_vec);
  // jnll_comp[0] += MVNORM(AQinvA_rw2)(Ax_rw2_2_vec);
  // jnll_comp[0] += MVNORM(AQinvA_rw2)(Ax_rw2_3_vec);
  // jnll_comp[0] += MVNORM(AQinvA_st)(Ax_st_vec);
  // jnll_comp[0] -= dnorm(Ax_survey_vec(0), Type(0.0), Type(sqrt(AQinvA_survey(0,0))), true); // N(mean, sd)

  // Add in pi(x)
  jnll_comp[0] += GMRF(Q_space)(epsilon_s);
  // jnll_comp[0] += GMRF(Q_struc_time)(age_rw2_1_vec);
  // jnll_comp[0] += GMRF(Q_struc_time)(age_rw2_2_vec);
  // jnll_comp[0] += GMRF(Q_struc_time)(age_rw2_3_vec);
  // jnll_comp[0] += GMRF(Q_struc_spacetime)(epsilon_st);
  // jnll_comp[0] += GMRF(Q_unstruc_survey)(survey_coefs_vec); // survey
  // for (int i = 0; i < S; i++) {
  //   jnll_comp[0] -= dnorm(Epsilon_t(i), Type(0.0), time_sd, true);
  // }
  // Type beta_sd = sqrt(1/exp(beta_log_tau));
  // for (int i = 0; i < N; i++) {
  //   jnll_comp[1] -= dnorm(beta(i), Type(0.0), beta_sd, true);
  // }
  Type cluster_sd = sqrt(1/exp(cluster_log_tau));
  for (int i = 0; i < nclusts; i++) {
    jnll_comp[1] -= dnorm(cluster_effect(i), Type(0.0), cluster_sd, true);
  }



  /////////
  // (2) //
  /////////



  // priors for intercepts
  for (int i = 0; i < A_m; i++) {
    jnll_comp[1] -= dnorm(alpha(i), Type(0.0), Type(31.62278), true); // corresponds to a precision of 0.001, as in INLA default
  }

  // hyperpriors for bym2 in space
  jnll_comp[1] -= dlogitbeta(space_logit_phi, Type(0.5), Type(0.5), true);
  // Type space_tau = exp(space_log_tau);
  // jnll_comp[1] -= dgamma(space_tau, Type(1.0), Type(1/0.00005), true);
  jnll_comp[1] -= dpcprec(space_log_tau, Type(1.0), Type(0.01), true);

  // hyperprior for IID re in time
  // Type time_tau = exp(time_log_tau);
  // jnll_comp[1] -= dlgamma(time_log_tau, Type(1.0), Type(1/0.00005), true);
  // jnll_comp[1] -= dpcprec(time_log_tau, Type(1.0), Type(0.01), true);

  // hyperprior for RW2s by age group
  // jnll_comp[1] -= dgamma(age_rw2_tau, Type(1.0), Type(1/0.00005), true);
  // jnll_comp[1] -= dpcprec(age_rw2_log_tau, Type(1.0), Type(0.01), true);
  
  // hyperprior for space/time interaction
  // jnll_comp[1] -= dgamma(spacetime_tau, Type(1.0), Type(1/0.00005), true);
  // jnll_comp[1] -= dpcprec(spacetime_log_tau, Type(1.0), Type(0.01), true);

  // hyperprior for beta (random slopes in time for each area)
  // Type beta_tau = exp(beta_log_tau);
  // jnll_comp[1] -= dgamma(beta_tau, Type(1.0), Type(1/0.00005), true); // Note that this is not the same pc prior as in SUMMER, which according to Richard does not have a closed form
  // jnll_comp[1] -= dpcprec(beta_log_tau, Type(1.0), Type(0.01), true); 

  // hyperprior for iid cluster random effect
  jnll_comp[1] -= dpcprec(cluster_log_tau, Type(1.0), Type(0.01), true);


  /////////
  // (3) //
  /////////


  int idx_space;
  int idx_age;
  int idx_cluster;
  int idx_sa;
  // Likelihood contribution from each datapoint ic
  for (int i = 0; i < nobs; i++){

    // figure out which region and time this cluster is in
    idx_space = region_id(i) - 1; // subtract one because c++ uses 0-indexing
    // int idx_time = time_id(i) - 1; // subtract one because c++ uses 0-indexing
    idx_age = age_id(i) - 1; // subtract one because c++ uses 0-indexing
    // int idx_sta = spacetimeage_id(i) - 1; // subtract one because c++ uses 0-indexing
    // int idx_st = spacetime_id(i) - 1; // subtract one because c++ uses 0-indexing
    idx_sa = spaceage_id(i) - 1;
    // int idx_age_group = age_group_id(i); // will be 1, 2, or 3.
    idx_cluster = cluster_id(i) - 1; // subtract one because c++ uses 0-indexing
    // int idx_survey = survey_id(i) - 1; // subtract one because c++ uses 0-indexing
    // int rw2_1 = 0; 
    // int rw2_2 = 0; 
    // int rw2_3 = 0;
    // if (idx_age_group == 1) {
    //   rw2_1 = 1;
    // }
    // if (idx_age_group == 2) {
    //   rw2_2 = 1;
    // }
    // if (idx_age_group == 3) {
    //   rw2_3 = 1;
    // }

    // calculate risk
    risk_ia(idx_sa) = exp(alpha(idx_age) + cluster_effect(idx_cluster) + total_latent_i(idx_space) + log(hiv_adj(0))) / (1 + exp(alpha(idx_age) + cluster_effect(idx_cluster) + total_latent_i(idx_space) + log(hiv_adj(0))));

    // and add data contribution to jnll
    jnll_comp[2] -= dbinom( y_iac(i), N_iac(i), risk_ia(idx_sa), true);

  } 

  // combine all parts of contributions to the jnll
  jnll += jnll_comp.sum();

  // ADREPORT section - currently unused
  // REPORT(Q_unstruc_time);

  return jnll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
