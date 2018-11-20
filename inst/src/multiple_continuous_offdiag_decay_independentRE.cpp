#include <TMB.hpp>
#include <cppad/cppad.hpp>
#include <Eigen/Dense>

using namespace Eigen;
using namespace density;

// Templates for fitting a two state Markov model to multiple individuals with decay formulation and re
// Form of each offdiagonal transition matirix is beta0 + beta_1 * exp(beta_2*cov) + re

// state_list, time_list, and covariate_list (of length one) are all templates that read a list in from R of
// states, times, and matricies of covariates respectively for
// use by the objective function defined below.
// Each element of these lists refers to one individual.

template<class Type>
struct state_list : vector<vector <Type> >  {
  state_list(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asVector<Type>(sm);
    }
  }
};

template<class Type>
struct time_list : vector<vector <Type> > {
  time_list(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asVector<Type>(sm);
    }
  }
};
template<class Type>
struct covariate_list : vector<vector <Type> > {
  covariate_list(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asVector<Type>(sm);
    }
  }
};

template<class Type>
Type objective_function<Type>::operator() (){
  DATA_FACTOR(ID);
  DATA_STRUCT(states, state_list); //an array of numeric states (i.e., in whale example 1s and 2s) for each whale
  DATA_STRUCT(times, time_list); //an array of times for each whale
  DATA_STRUCT(covariates, covariate_list); //an array covriate vectors for each whale
  vector<Type> q(2); // declare q
  PARAMETER_VECTOR(log_baseline);//a vector of the log off diagonal transition baselines
  vector<Type> baseline = exp(log_baseline); // declare intercept
  PARAMETER_MATRIX(betas_matrix); // coefficients for the intensity jump and exp decay
  Type logb1_12 = betas_matrix(0,0); Type logb1_21 = betas_matrix(1,0);
  Type b1_12 = exp(logb1_12); Type b1_21 = exp(logb1_21);
  Type b2_12 = -exp(betas_matrix(0,1)); Type b2_21 = -exp(betas_matrix(1,1));
  // Declaring random effects
  PARAMETER_MATRIX(u);
  // random effect log sigma
  PARAMETER(log_sigma);
  Type sigma = exp(log_sigma);
  int wh =  NLEVELS(ID); // number of whales
  // Initialize log-likelihood variable for parallel summation:
  parallel_accumulator<Type> ll(this);
  matrix<Type> Q(2,2); // declare transition matrix
  //contribution from observed data
  for (int j = 0; j < wh; j++){
    vector<Type> tem = times(j);
    vector<Type> sem = states(j); // times, states,
    vector<Type> covs = covariates(j); //covariates
    int t = tem.size();
    for (int i = 0; i < (t-1); i++){
      // MVN latent variables u for each individual j
      q(0) = exp(log_baseline(0) + logb1_12*exp(b2_12*covs(i)) + u(j,0));
      q(1) = exp(log_baseline(1) + logb1_21*exp(b2_21*covs(i)) + u(j,1));
      Q(0,0) = - q(0); Q(0,1) = q(0); Q(1,0) = q(1); Q(1,1) = -q(1); 
      Type temp = tem(i+1) - tem(i);
      int x = CppAD::Integer(sem(i));
      int y = CppAD::Integer(sem(i+1));
      matrix<Type> Qt = Q*temp;
      matrix<Type> P = atomic::expm(Qt); // Prob transition matrix
      Type p = P(x-1,y-1);
      ll -= log(p);
    }
    ll -= dnorm(u(j,0), Type(0), sigma,true); // contribution from 1--2 transition for individual j
    ll -= dnorm(u(j,1), Type(0), sigma,true); // contribution from 2--1 transition for individual j
  }
  ADREPORT(Q); ADREPORT(b1_12); ADREPORT(b1_21); ADREPORT(b2_12); ADREPORT(b2_21); ADREPORT(sigma);
  return ll;
}
