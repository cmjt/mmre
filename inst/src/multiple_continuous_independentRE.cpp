#include <TMB.hpp>
#include <cppad/cppad.hpp>
#include <Eigen/Dense>

using namespace Eigen;
using namespace density;

// Templates for fitting a two state Markov model to multiple individuals with a transition matrix Q
// and an independent MVN random effect on the transition matrix elements state 1 -> 2 and state 2 -> 1 for
// each individual.

// state_list and time_list are both templates that read a list in from R of
// states and times respectively for use by the objective function defined below.
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
Type objective_function<Type>::operator() (){
  DATA_FACTOR(ID);
  DATA_STRUCT(states, state_list); //an array of numeric states (i.e., in whale example 1s and 2s) for each whale
  DATA_STRUCT(times, time_list); //an array of times for each whale
  PARAMETER_VECTOR(log_baseline);//a vector of the log off diagonal transition intensities to be estimated using ML
  vector<Type> q(2); // declare q
  // Declaring random effects
  PARAMETER_MATRIX(u);
  // random effect log sigma
  PARAMETER(log_sigma);
  Type sigma = exp(log_sigma);
  int wh =  NLEVELS(ID); // number of whales
  Type ll(0);
  matrix<Type> Q(2,2); // declare transition matrix
  //contribution from observed data
  for (int j = 0; j < wh; j++){
    vector<Type> tem = times(j);
    vector<Type> sem = states(j);
    q(0) = exp(log_baseline(0) + u(j,0));
    q(1) = exp(log_baseline(1) + u(j,1));
    Q(0,0) = - q(0); Q(0,1) = q(0); Q(1,0) = q(1); Q(1,1) = -q(1); 
    int t = tem.size();
      for (int i = 0; i < (t-1); i++){
      	Type temp = tem(i+1) - tem(i);
	int x = CppAD::Integer(sem(i));
	int y = CppAD::Integer(sem(i+1));
      	matrix<Type> Qt = Q*temp;
      	matrix<Type> P = atomic::expm(Qt); // Prob transition matrix
	Type p = P(x-1,y-1);
	ll -= log(p);
      }
      // MVN latent variables u for each individual j
      ll -= dnorm(u(j,0), Type(0), sigma,true); // contribution from 1--2 transition for individual j
      ll -= dnorm(u(j,1), Type(0), sigma,true); // contribution from 2--1 transition for individual j
  }
  ADREPORT(Q); ADREPORT(sigma);
  return ll;
}
