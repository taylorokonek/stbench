// Normal model: 
//     (single time point assumed)
//     Single intercepts, 
//     BYM2 for space, 
//     Benchmarked

#ifndef normal_intercepts_bym2_hpp
#define normal_intercepts_bym2_hpp

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
Type normal_intercepts_bym2(objective_function<Type>* obj)
{

  // Data and parameters passed in from R

  // Total number of space/time/age/survey/cluster 
  DATA_INTEGER( N );   // Number of spatial regions
  DATA_INTEGER( nclusts ); // Number of clusters in the dataset - corresponds to number of rows
  DATA_INTEGER( nobs ); // Number of observations

  // Data
  DATA_VECTOR( y_ic );   // Number of deaths in space/time/age/cluster
  DATA_VECTOR( SE_i );   // Fixed standard deviation for each area
  DATA_VECTOR( hiv_adj );
  
  // Ids for space/time/age/survey
  DATA_IVECTOR( region_id ); // Numeric indicator for which region an observation belongs to
  DATA_IVECTOR( cluster_id );

  // Benchmarking things
  // DATA_VECTOR( natl ); // National-level estimates we want to benchmark to, in order of time
  // DATA_VECTOR( natl_sd ); // SDs for national-level estimates, in order of time
  // DATA_VECTOR( pop_weights ); // population weights for benchmarking constraint. Must sum to 1 at each time point and be in order arrange(region, time)
  
  // Precision matrices for space (structured), time (unstructured and structured), survey
  DATA_SPARSE_MATRIX( Q_struc_space ); // structured precision matrices must be scaled

  // 'Eigenvalues' of Q_struc_space^{-}, for use in pc prior for phi
  DATA_VECTOR( gamma_tilde ); 

  // Ids for space/time/age in risk - used for benchmarking 
  // DATA_IVECTOR( risk_region_id ); // Numeric indicator for which region an observation belongs to

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
  DATA_VECTOR( alpha_pri ); // prior for intercept alpha, with alpha ~ N(alpha_pri[0], alpha_pri[1]) (mean, sd)

  // Fixed effects & hyperpars
  PARAMETER( alpha );           // Intercept
  PARAMETER( space_logit_phi );          // logit of BYM2 mixing parameter
  PARAMETER( space_log_tau );            // Log of INLA tau param (precision of space covariance matrix)
  PARAMETER( cluster_log_tau );          // Log of tau param for cluster random effect

  // Random effects
  PARAMETER_ARRAY( Epsilon_s );    // effects for each spatial region (length 2N). first half is total, second half is structured
  PARAMETER_ARRAY( cluster_effect );



  /////////////////////////////////
  // Define internal objects     //
  /////////////////////////////////



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
  // vector<Type> temp_risk_i( N );
  vector<Type> risk_i( N );          // expit(total latent field) risk for each region

  // Initialize vectors for survey, space (also total and structured components), spacetime, rw2s by age group...
  vector<Type> epsilon_s( 2*N );     // joint total (followed by) structured effects vec
  vector<Type> total_latent_i( N );  // struc + unstruc in space
  vector<Type> struc_latent_i( N );  // struc in space

  // Fill in vectors for survey, space, spacetime, rw2s by age group...
  for(int s = 0; s < N; s++){
    epsilon_s[s] = Epsilon_s(s);
    epsilon_s[N + s] = Epsilon_s(N + s);
  }



  //////////////////////////////////////////////////////
  // Constrain things that need to be constrained     //
  //////////////////////////////////////////////////////



  // a constrained vector x_c for an unconstrained vector x is calculated as:
  // x_c = x - Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e)
  // for constraint Ax=e, and for GMRF x with precision Q
  // we'll build up the matrices we need piece by piece

  // Invert precision matrices (possible since we added a small value to the diagonal)
  matrix<Type> Q_inv_space = invertSparseMatrix(Q_space);

  // Create constraint matrices A
  matrix<Type> A_space(1, full_N); // sum all components
  for(int i = 0; i < N; i++){
    A_space(0, i)     = 0; // no constraints on total effects
    A_space(0, N + i) = 1; // sum-to-0 on structured effects
  }

  // Create A^T
  matrix<Type> A_space_T = A_space.transpose(); // bym2 in space

  // Create Q^{-1}A^T
  matrix<Type> QinvA_space_mat = Q_inv_space * A_space_T;

  // For precision matrices with sum-to-zero constraints (survey, space), convert matrix to vector
  vector<Type> QinvA_space(full_N);
  for(int i = 0; i < full_N; i++) {
    QinvA_space(i) = QinvA_space_mat(i,0);
  }

  // Create AQ^{-1}A^T
  matrix<Type> AQinvA_space = A_space * QinvA_space_mat;

  // Create (AQ^{-1}A^T)^{-1}
  matrix<Type> AQinvA_space_inv = AQinvA_space.inverse(); // okay for small matrices

  // Create Ax
  matrix<Type> Ax_space = (A_space * epsilon_s.matrix());

  // Convert Ax from matrix to vector form - needed for dnorm & MVNORM, so only needed for survey fixed effect, rw2s by age group, and space/time interaction
  vector<Type> Ax_space_vec(1);
  Ax_space_vec(0) = Ax_space(0, 0);
 
  // Convert Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e) to vector form for conditioning by kriging correction
  matrix<Type> krig_correct_space = QinvA_space_mat * AQinvA_space_inv * Ax_space;
  vector<Type> krig_correct_space_vec(full_N);
  for (int i = 0; i < full_N; i++) {
    krig_correct_space_vec(i) = krig_correct_space(i,0);
  }

  // Construct constrained vector x_c = x - Q^{-1}A'(AQ^{-1}A')^{-1}(Ax-e)
  vector<Type> epsilon_s_c = epsilon_s - krig_correct_space_vec;

  // For BYM2s (only space at the moment), fill in total and struc vectors for use in risk calculation
  for(int s = 0; s < N; s++){
    total_latent_i[s] = epsilon_s_c(s);     // first half is total
    struc_latent_i[s] = epsilon_s_c(N + s); // second half is structured effect
  }

  // Multiply precision matrices by their precision parameters
  // Not needed for any random effects in this model
 


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

  // Add in pi(x)
  jnll_comp[0] += GMRF(Q_space)(epsilon_s);
  Type cluster_sd = sqrt(1/exp(cluster_log_tau));
  for (int i = 0; i < nclusts; i++) {
    jnll_comp[1] -= dnorm(cluster_effect(i), Type(0.0), cluster_sd, true);
  }



  /////////
  // (2) //
  /////////



  // prior for intercepts
  jnll_comp[1] -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true); // corresponds to a precision of 0.001, as in INLA default

  // hyperpriors for bym2 in space
  jnll_comp[1] -= dlogitbeta(space_logit_phi, Type(0.5), Type(0.5), true);
  // jnll_comp[1] -= dpcphi(space_logit_phi, gamma_tilde, Type(0.5), Type(2/3), true);

  // Type space_tau = exp(space_log_tau);
  // jnll_comp[1] -= dgamma(space_tau, Type(1.0), Type(1/0.00005), true);
  jnll_comp[1] -= dpcprec(space_log_tau, Type(1.0), Type(0.01), true);

  // hyperprior for iid cluster random effect
  jnll_comp[1] -= dpcprec(cluster_log_tau, Type(1.0), Type(0.01), true);
  


  /////////
  // (3) //
  /////////



  // Likelihood contribution from each datapoint ic
  for (int i = 0; i < nobs; i++){

    // figure out which region and time this cluster is in
    int idx_space = region_id(i) - 1; // subtract one because c++ uses 0-indexing
    int idx_cluster = cluster_id(i) - 1; // subtract one because c++ uses 0-indexing

    // calculate risk
    risk_i(idx_space) = alpha + cluster_effect(idx_cluster) + total_latent_i(idx_space) + log(hiv_adj(0));

    // and add data contribution to jnll
    jnll_comp[2] -= dnorm(y_ic(i), risk_i(idx_space), SE_i(idx_space), true);

  } 

  // combine all parts of contributions to the jnll
  jnll += jnll_comp.sum();

  // ADREPORT section - currently unused
  // REPORT(Q_struc_space);

  return jnll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif