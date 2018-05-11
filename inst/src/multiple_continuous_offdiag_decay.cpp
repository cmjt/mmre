#include <TMB.hpp>
#include <cppad/cppad.hpp>
#include <Eigen/Dense>

using namespace Eigen;
using namespace density;

// Templates for fitting a two state Markov model to multiple individuals with a single covariate
//  Form of each offdiagonal transition matirix is beta0 + beta_1 * exp(beta_2*cov) +

// state_list, time_list, and covariate_list (of length one) are all templates that read a list in from R of
// states, times, and covariates respectively for
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
  PARAMETER_MATRIX(betas_matrix); // the coefficients for the covariates for transition 1->2 and 2->1
  int wh =  NLEVELS(ID); // number of whales
  Type ll = 0; //declare log-likelihood
  matrix<Type> Q(2,2); // declare transition matrix
  //contribution from observed data
  for (int j = 0; j < wh; j++){
    vector<Type> tem = times(j);
    vector<Type> sem = states(j); // times, states,
    vector<Type> covs = covariates(j); //covariates
    vector<Type> beta0 = betas_matrix.row(0); // transiions 1--2
    vector<Type> beta1 = betas_matrix.row(1); // transitions 2--1
    int t = tem.size();
      for (int i = 0; i < (t-1); i++){
	// MVN latent variables u for each individual j
	q(0) = exp(beta0(0) + beta0(1)*exp(-exp(beta0(2))*covs(i)));
	q(1) = exp(beta1(0) + beta1(1)*exp(-exp(beta1(2))*covs(i)));
	Q(0,0) = - q(0); Q(0,1) = q(0); Q(1,0) = q(1); Q(1,1) = -q(1); 
      	Type temp = tem(i+1) - tem(i);
	int x = CppAD::Integer(sem(i));
	int y = CppAD::Integer(sem(i+1));
      	matrix<Type> Qt = Q*temp;
      	matrix<Type> P = atomic::expm(Qt); // Prob transition matrix
	Type p = P(x-1,y-1);
	ll += log(p);
      }
      vector<Type> q0(2); // declare q0
      q0(0) = exp(beta0(0)); q0(1) = exp(beta1(0));
      vector<Type> beta_1(2);
      beta_1(0) = beta0(1); beta_1(1) = beta1(1);
      vector<Type> beta_2(2);
      beta_2(0) = -exp(beta0(2)); beta_2(1) = -exp(beta1(2));
      ADREPORT(q0);
      ADREPORT(beta_1);
      ADREPORT(beta_2);
  }
  return -ll;
}
