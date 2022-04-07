// Various helper functions for models:
//     d_func: used in pc prior density for phi parameter in BYM2
//     bym2_Q: create joint (sparse) matrix of BYM2 joint precision, as in riebleretal2016

// include libraries

#ifndef helper_funcs_hpp
#define helper_funcs_hpp 1

#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function to create joint (sparse) matrix of the BYM2 joint precision
template<class Type>
//SparseMatrix<Type> bym2_Q(Type logit_phi, Type log_prec,
SparseMatrix<Type> bym2_Q(Type logit_phi, Type log_prec,
                          SparseMatrix<Type> Q_struc) {

  int N = int(Q_struc.rows()); // number of spatial regions
  int full_N = N + N;      // top half is total effect, bottom half is structured effect
  Type tau = exp(log_prec);
  Type phi = exp(logit_phi)/(1 + exp(logit_phi));

  SparseMatrix<Type> Q(full_N, full_N);

  // block 11
  for(int i = 0; i < N; i++){
    Q.coeffRef(i, i) = (tau / (1 - phi));
  }

  // block 12 and block 21
  for(int i = 1; i < N; i++){
    Q.coeffRef(i, N+i) = (-sqrt(phi * tau) / (1 - phi)); // diag of block 12
    Q.coeffRef(N+i, i) = (-sqrt(phi * tau) / (1 - phi)); // diag of block 21
  }

  // block 22

  // first, correct diagonal as needed for the joint specification
  for(int i = 0; i < N; i++){
    Q_struc.coeffRef(i, i) = Q_struc.coeffRef(i, i) + (phi / (1-phi));
  }
  // and then fill it into the joint precision
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      Q.coeffRef(N + i, N + j) = Q_struc.coeffRef(i, j);
    }
  }

  return Q;
}

#endif