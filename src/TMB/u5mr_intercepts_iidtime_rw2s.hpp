// U5MR model: 
//     Age-specific intercepts, 
//     BYM2 for space, 
//     RW2s in time for three age groups, 
//     random area-specific slopes in time, 
//     Knorr-Held Type IV Space-Time Interaction

#ifndef u5mr_intercepts_iidtime_rw2s_hpp
#define u5mr_intercepts_iidtime_rw2s_hpp

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
Type u5mr_intercepts_iidtime_rw2s(objective_function<Type>* obj)
{

  // Data and parameters passed in from R

  // Total number of space/time/age/survey/cluster 
  DATA_INTEGER( S );   // Number of time points
  DATA_INTEGER( A_m ); // Number of age groups
  DATA_INTEGER( nclusts ); // Number of clusters 
  DATA_INTEGER( nobs ); // Number of observations 

  // Data
  DATA_VECTOR( y_tac );   // Number of deaths in time/age/cluster
  DATA_VECTOR( N_tac );   // Ntrials for binomial in time/age/cluster
  DATA_VECTOR( hiv_adj );  // HIV adjustment ratios - log offsets for each time point. If no offset wanted, hiv_adj can be a vector of 1's
  
  // Ids for space/time/age/survey
  DATA_IVECTOR( time_id ); // Numeric indicator for which time point an observation belongs to
  DATA_IVECTOR( age_id ); // Numeric indicator for which age group an observation belongs to
  DATA_IVECTOR( timeage_id ); // Numeric indicator for which time/age point an observation belongs to
  DATA_IVECTOR( age_group_id ); // Numeric indicator for which age group an observation belongs to
  DATA_IVECTOR( cluster_id );

  
  // Precision matrices for time (structured)
  DATA_SPARSE_MATRIX( Q_struc_time ); // structured precision matrices must be scaled

  // Constraint matrices

  // Prior specifications - currently unused. 

  // Fixed effects & hyperpars
  PARAMETER_ARRAY( alpha );           // Intercept. An array because there are different intercepts for each age group, a A_m x 1 array
  PARAMETER( time_log_tau );            // Log of tau param (precision of temporal covariance matrix - the iid bit)
  PARAMETER( age_rw2_log_tau );          // Log of tau param (precision on temporal covariance matrix - the rw2s for age groups)
  PARAMETER( cluster_log_tau );          // Log of tau param for cluster random effect

  // Random effects
  PARAMETER_ARRAY( Epsilon_t );    // effects for each time point.
  PARAMETER_ARRAY( age_rw2_1 );    // rw2s for 3 different age groups
  PARAMETER_ARRAY( age_rw2_2 );
  PARAMETER_ARRAY( age_rw2_3 );
  PARAMETER_ARRAY( cluster_effect );



  /////////////////////////////////
  // Define internal objects     //
  /////////////////////////////////

  // Multiply certain precision matrices by their precision parameters
  // Not needed for anything BYM2 since that is done in the bym2_Q() function
  // Not needed for Q_unstruc_survey since we fix the precision here and have it multiplied in already in R (prior to tmb)
  // We'll do this using the SCALE() function when adding in contributions to the jnll

  // Q_struc_time
  Type age_rw2_tau = exp(age_rw2_log_tau);

  // Q_unstruc_time
  Type time_sd = 1/sqrt(exp(time_log_tau));

  // Multiply precision matrices by their precision parameters
  SparseMatrix<Type> Q_struc_time_internal(S, S);
  for (int i = 0; i < S; i++) {
    for (int j = 0; j < S; j++) {
      Q_struc_time_internal.coeffRef(i, j) = Q_struc_time.coeffRef(i, j) * age_rw2_tau;
    }
  }

  // Set up objective function -- joint negative log-likelihood
  // Define it in three parts (just for ease of reading code): (1) priors (2) space/time/age fields (3) data contribution
  vector<Type> jnll_comp(3);
  Type jnll = 0; // sum of jnll_comp

  // Initialize vectors for risk (space/time/age and space/time)
  vector<Type> risk_ta( S * A_m );          // expit(total latent field) risk for each time/age group point, arranged in order arrange(region, time, age)

  // Initialize vectors for survey, space (also total and structured components), spacetime, rw2s by age group...
  vector<Type> age_rw2_1_vec( S ); // rw2s by age group
  vector<Type> age_rw2_2_vec( S );
  vector<Type> age_rw2_3_vec( S );

  // Fill in vectors for rw2s by age group...
  for(int s = 0; s < S; s++){
    age_rw2_1_vec[s] = age_rw2_1(s);
    age_rw2_2_vec[s] = age_rw2_2(s);
    age_rw2_3_vec[s] = age_rw2_3(s);
  }



  //////////////////////////////////////////////////////
  // Constrain things that need to be constrained     //
  //////////////////////////////////////////////////////



  // a constrained vector x_c for an unconstrained vector x is calculated as:
  // x_c = x - Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e)
  // for constraint Ax=e, and for GMRF x with precision Q
  // we'll build up the matrices we need piece by piece

  // Invert precision matrices (possible since we added a small value to the diagonal)
  matrix<Type> Q_inv_rw2 = invertSparseMatrix(Q_struc_time_internal);

  // Create constraint matrices A
  matrix<Type> A_rw2(1, S);
  for(int i = 0; i < S; i++) {
    A_rw2(0, i) = 1; // sum-to-0 constraint
    // A_rw2(1, i) = i + 1; // constraint on slopes (1, 2, ..., S)
  }

  // Create A^T
  matrix<Type> A_rw2_T = A_rw2.transpose();

  // Create Q^{-1}A^T
  matrix<Type> QinvA_rw2= Q_inv_rw2 * A_rw2_T;

  // Create AQ^{-1}A^T
  matrix<Type> AQinvA_rw2 = A_rw2 * QinvA_rw2;

  // Create (AQ^{-1}A^T)^{-1}
  matrix<Type> AQinvA_rw2_inv = AQinvA_rw2.inverse(); // okay for small matrices

  // Create Ax
  matrix<Type> Ax_rw2_1 = (A_rw2 * age_rw2_1_vec.matrix());
  matrix<Type> Ax_rw2_2 = (A_rw2 * age_rw2_2_vec.matrix());
  matrix<Type> Ax_rw2_3 = (A_rw2 * age_rw2_3_vec.matrix());

  // Convert Ax from matrix to vector form - needed for dnorm & MVNORM, so only needed for survey fixed effect, rw2s by age group, and space/time interaction
  vector<Type> Ax_rw2_1_vec(1);
  vector<Type> Ax_rw2_2_vec(1);
  vector<Type> Ax_rw2_3_vec(1);
  for(int i = 0; i < 1; i++) {
    Ax_rw2_1_vec(i) = Ax_rw2_1(i,0);
    Ax_rw2_2_vec(i) = Ax_rw2_2(i,0);
    Ax_rw2_3_vec(i) = Ax_rw2_3(i,0);
  }
 
  // Convert Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e) to vector form for conditioning by kriging correction
  matrix<Type> krig_correct_1 = QinvA_rw2 * AQinvA_rw2_inv * Ax_rw2_1;
  matrix<Type> krig_correct_2 = QinvA_rw2 * AQinvA_rw2_inv * Ax_rw2_2;
  matrix<Type> krig_correct_3 = QinvA_rw2 * AQinvA_rw2_inv * Ax_rw2_3;
  vector<Type> krig_correct_1_vec(S);
  vector<Type> krig_correct_2_vec(S);
  vector<Type> krig_correct_3_vec(S);
  for (int i = 0; i < S; i++) {
    krig_correct_1_vec(i) = krig_correct_1(i,0);
    krig_correct_2_vec(i) = krig_correct_2(i,0);
    krig_correct_3_vec(i) = krig_correct_3(i,0);
  }

  // Construct constrained vector x_c = x - Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e)
  vector<Type> age_rw2_1_c = age_rw2_1_vec - krig_correct_1_vec;
  vector<Type> age_rw2_2_c = age_rw2_2_vec - krig_correct_2_vec;
  vector<Type> age_rw2_3_c = age_rw2_3_vec - krig_correct_3_vec;



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
  // jnll_comp[0] += MVNORM(AQinvA_rw2)(Ax_rw2_1_vec);
  // jnll_comp[0] += MVNORM(AQinvA_rw2)(Ax_rw2_2_vec);
  // jnll_comp[0] += MVNORM(AQinvA_rw2)(Ax_rw2_3_vec);

  // Add in pi(x)
  jnll_comp[0] += GMRF(Q_struc_time_internal)(age_rw2_1_vec);
  jnll_comp[0] += GMRF(Q_struc_time_internal)(age_rw2_2_vec);
  jnll_comp[0] += GMRF(Q_struc_time_internal)(age_rw2_3_vec);

  for (int i = 0; i < S; i++) {
    jnll_comp[0] -= dnorm(Epsilon_t(i), Type(0.0), time_sd, true);
  }

  Type cluster_sd = sqrt(1 / exp(cluster_log_tau));
  // Type cluster_sd = sqrt(1/ Type(22.81349));
  for (int i = 0; i < nclusts; i++) {
    // jnll_comp[0] -= dnorm(cluster_effect(i), Type(0.0), cluster_sd, true);
    jnll_comp[0] -= dnorm(cluster_effect(i), Type(0.0), cluster_sd, true);
  }



  /////////
  // (2) //
  /////////



  // priors for intercepts
  for (int i = 0; i < A_m; i++) {
    jnll_comp[1] -= dnorm(alpha(i), Type(0.0), Type(31.62278), true); // corresponds to a precision of 0.001, as in INLA default
  }

  // hyperprior for IID re in time
  jnll_comp[1] -= dpcprec(time_log_tau, Type(1.0), Type(0.01), true);
  // Type time_tau = exp(time_log_tau);
  // jnll_comp[1] -= dlgamma(time_log_tau, Type(1.0), Type(1/0.00005), true);

  // hyperprior for RW2s by age group
  jnll_comp[1] -= dlgamma(age_rw2_log_tau, Type(1.0), Type(1/0.00005), true);
  // jnll_comp[1] -= dpcprec(age_rw2_log_tau, Type(1.0), Type(0.01), true);
  
  // hyperprior for iid cluster random effect
  jnll_comp[1] -= dpcprec(cluster_log_tau, Type(1.0), Type(0.01), true);
  // Type cluster_tau = exp(cluster_log_tau);
  // jnll_comp[1] -= dlgamma(cluster_log_tau, Type(1.0), Type(1.0), true);


  /////////
  // (3) //
  /////////



  // Likelihood contribution from each datapoint ic
  int idx_time;
  int idx_age;
  int idx_ta;
  int idx_age_group;
  int idx_cluster;
  int rw2_1;
  int rw2_2;
  int rw2_3;
  Type logit_risk;
  for (int i = 0; i < nobs; i++){

    // figure out which time this cluster is in
    idx_time = time_id(i) - 1; // subtract one because c++ uses 0-indexing
    idx_age = age_id(i) - 1; // subtract one because c++ uses 0-indexing
    idx_ta = timeage_id(i) - 1; // subtract one because c++ uses 0-indexing
    idx_age_group = age_group_id(i); // will be 1, 2, or 3.
    idx_cluster = cluster_id(i) - 1; // subtract one because c++ uses 0-indexing
    rw2_1 = 0; 
    rw2_2 = 0; 
    rw2_3 = 0;
    if (idx_age_group == 1) {
      rw2_1 = 1;
    }
    if (idx_age_group == 2) {
      rw2_2 = 1;
    }
    if (idx_age_group == 3) {
      rw2_3 = 1;
    }

    // calculate risk
    logit_risk = alpha(idx_age) + cluster_effect(idx_cluster) + rw2_1 * age_rw2_1_c(idx_time) + rw2_2 * age_rw2_2_c(idx_time) + rw2_3 * age_rw2_3_c(idx_time) + Epsilon_t(idx_time); //+ log(hiv_adj(idx_time));
    risk_ta(idx_ta) = exp(logit_risk)/(1 + exp(logit_risk));

    // and add data contribution to jnll
    jnll_comp[2] -= dbinom( y_tac(i), N_tac(i), risk_ta(idx_ta), true);

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
