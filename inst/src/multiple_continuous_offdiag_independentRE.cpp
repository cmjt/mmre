#include <TMB.hpp>
#include <cppad/cppad.hpp>
#include <Eigen/Dense>

using namespace Eigen;
using namespace density;

// Templates for fitting a two state Markov model to multiple individuals with a multiple covariates
// on the transition matrix elements state 1 -> 2 and state 2 -> 1 for
// each individual. 

// state_list, time_list, and covariate_list are all templates that read a list in from R of
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
struct covariate_list : vector<matrix <Type> > {
  covariate_list(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

template<class Type>
Type objective_function<Type>::operator() (){
  DATA_FACTOR(ID);
  DATA_STRUCT(states, state_list); //an array of numeric states (i.e., in whale example 1s and 2s) for each whale
  DATA_STRUCT(times, time_list); //an array of times for each whale
  DATA_STRUCT(covariates, covariate_list); //an array covriate matricies for each whale refering to each spline
  vector<Type> q(2); // declare q
  PARAMETER_VECTOR(log_baseline);//a vector of the log off diagonal transition baselines
  vector<Type> baseline = exp(log_baseline); // declare intercept
  PARAMETER_MATRIX(betas_matrix); // the coefficients for the covariates for transition 1->2 and 2->1
  vector<Type> log_coef1_2 = betas_matrix.row(0);
  vector<Type> log_coef2_1 = betas_matrix.row(1);
  vector<Type> coef1_2 = exp(log_coef1_2);
  vector<Type> coef2_1 = exp(log_coef2_1);
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
    vector<Type> covariates1_2 = covariates(j)*log_coef1_2;
    vector<Type> covariates2_1 = covariates(j)*log_coef2_1;
    int t = tem.size();
      for (int i = 0; i < (t-1); i++){
	q(0) = baseline(0)*exp(covariates1_2(i) + u(j,0));
        q(1) = baseline(1)*exp(covariates2_1(i) + u(j,1));
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
  ADREPORT(Q);ADREPORT(coef1_2);ADREPORT(coef2_1);// HR and baseline coefs
  ADREPORT(sigma);
  return ll;
}
