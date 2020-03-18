#include <TMB.hpp>
#include <cppad/cppad.hpp>
#include <Eigen/Dense>

using namespace Eigen;
// Templates for fitting an N state Markov model to multiple individuals with transition matrix Q
// NO random effects

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
  DATA_STRUCT(times, time_list); //an array of times for each individual
  DATA_INTEGER(nstates); // number of states
  PARAMETER_ARRAY(qarray); // Q matrix (diagonal entries will be ignored)
  matrix<Type> Q(nstates,nstates); // make Q
  for (int k = 0; k < nstates; k++){
    // upper diag
    for(int l = 1; l < nstates; l++){
      if(k < l){
	Type val = qarray(l-1,0,k);
    	Q(k,l) = val;
      }
    }
    // lower diag
    for(int l = 0; l < (nstates - 1); l++){
      if(k > l){
    	Type val = qarray(l,0,k);
    	Q(k,l) = val;
      }
    }
    // diagonal
    Type sum = Q.row(k).sum();
    Q(k,k) = -sum;
  }
  int wh =  NLEVELS(ID); // number of individuals
  vector<Type> ll(wh);
  //contribution from observed data
  for (int j = 0; j < wh; j++){
    vector<Type> tem = times(j);
    vector<Type> sem = states(j); 
    int t = tem.size();
    for (int i = 0; i < (t-1); i++){
      Type temp = tem(i+1) - tem(i);
      	matrix<Type> Qt = Q*temp;
      	matrix<Type> P = atomic::expm(Qt);
	int x = CppAD::Integer(sem(i));
	int y = CppAD::Integer(sem(i+1));
	Type p = P(x-1,y-1);
	ll(j) += log(p);
      }
  }
  //ADREPORT(Q);
  return -sum(ll);
}
