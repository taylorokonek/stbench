// Additional densities for hyperpriors and priors:
//     dpcphi: pc prior for phi parameter in BYM2
//     dpcprec: pc prior for precision parameter in mvnormal
//     dlogitbeta: logit beta prior on theta = logit(probability) s.t. probability = p ~ Beta(a, b)
//     dlogtgaussian: prior on log(precision) of multi-var normal where stdev sigma ~ N(u, s) & sigma > 0

// pc prior density for phi parameter in BYM2
#ifndef addtl_densities_hpp
#define addtl_densities_hpp 1

#include <Eigen/Sparse>
#include <vector>

// helper function d_func(phi) for use in pc prior density for phi parameter in BYM2
template<class Type>
Type d_func(Type a, vector<Type> gamma_tilde) {
  Type part1 = a * (gamma_tilde - 1).sum();
  Type part2 = (log(a * (gamma_tilde - 1) + 1)).sum();
  return sqrt(part1 - part2); 
}

template<class Type>
Type dpcphi(Type logit_phi, vector<Type> gamma_tilde, Type U, Type alpha, int give_log = 0) {

  // set up
  Type phi = exp(logit_phi)/(1 + exp(logit_phi));
  Type lambda = -log(1 - alpha) / d_func(U, gamma_tilde);

  // compute log density
  Type part1 = log(lambda) + logit_phi - log(1 + exp(logit_phi)) - log(2);
  Type part2 = -log(d_func(phi, gamma_tilde)) - lambda * d_func(phi, gamma_tilde);
  vector<Type> part3_numerator(gamma_tilde.size());
  for (int i = 0; i < gamma_tilde.size(); i++) {
    part3_numerator(i) = pow(gamma_tilde(i) - 1, 2);
  }
  vector<Type> part3_denominator = 1 + phi * (gamma_tilde - 1);
  Type temp = (part3_numerator/part3_denominator).sum();
  if (temp < 0) {
    temp = temp * (-1);
  }
  Type part3 = log( temp );
  Type part4 = logit_phi - 2 * log(1 + exp(logit_phi));

  Type logres = part1 + part2 + part3 + part4;

  if(give_log)return logres; else return exp(logres);
}

// pc prior for normal precisions
template<class Type>
Type dpcprec(Type log_tau, Type U, Type alpha, int give_log = 0) {
  Type lambda = -log(alpha) / U;

  Type logres = log(lambda) - log(2) - 3 * log_tau/2 - lambda * pow(exp(log_tau), -1/2);

  if(give_log)return logres; else return exp(logres);
}

// helper function for logit beta prior on theta = logit(probability)
//   s.t. probability = p ~ Beta(a, b)
template<class Type>
Type dlogitbeta(Type logit_p, Type a, Type b, int give_log=0)
{
  Type part1 = lgamma(a + b) - lgamma(a)  - lgamma(b);
  Type part2 = (a - 1) * (logit_p - log(1 + exp(logit_p)));
  Type part3 = (b - 1) * log( 1 - exp(logit_p)/(1 + exp(logit_p)));
  Type part4 =  logit_p - 2 * log( 1 + exp(logit_p));

  Type logres = part1 + part2 + part3 + part4;

  if(give_log)return logres; else return exp(logres);
}

// helper function for pc prior on log(precision) of multi-var normal
// where stdev sigma ~ N(u, s) & sigma > 0
// internally, INLA places the prior on log(precision)
template<class Type>
Type dlogtgaussian(Type log_prec, Type u, Type s, int give_log=0)
{
  Type part1 = -0.5 * log(8 * M_PI) - log(s);
  Type part2 = -0.5 * 1 / pow(s, 2) * pow( (exp(-log_prec / 2) - u ), 2) - log_prec / 2;

  Type logres = part1 + part2;

  if(give_log)return logres; else return exp(logres);
}

// helper function for betabinomial distribution. takes parameters (Ntrials, risk, overdispersion), 
// transforms these to (Ntrials, alpha, beta) and goes from there
//   s.t. probability = x ~ Betabinom(Ntrials, risk, overdispersion)
template<class Type>
Type dbetabinom(Type n, Type k, Type p, Type d, int give_log=0)
{
  Type a = p/d - p; // alpha
  Type b = (1-p)/d - 1 + p; // beta
  // basically just a ton of lgamma terms all together
  Type part1 = lgamma(n + 1) - lgamma(k + 1)  - lgamma(n - k + 1);
  Type part2 = lgamma(k + a) + lgamma(n - k + b)  - lgamma(n + a + b);
  Type part3 = lgamma(a + b) - lgamma(a)  - lgamma(b);
  
  Type logres = part1 + part2 + part3;
  
  if(give_log)return logres; else return exp(logres);
}

#endif
