/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "SpecialFunctions.h"
#include "../math_meigen/mEigen.h"
#include "../math_random/librandom.h"


//================== Functions ==========================

/// liblibra namespace
namespace liblibra{

/// libspecialfunctions namespace
namespace libspecialfunctions{


MATRIX mean(MATRIX& X){
/** 
  Compute the mean of the data in each row of X:

  mean_i = <X_i> = (1/nsampl) sum_k { X_ik }

  X (ndof x nsampl)
  where ndof - the number of degrees of freedom 
        nsampl - the number of samples on ndof-dimensional data

*/

  int ndof = X.n_rows;
  int sz = X.n_cols;

  MATRIX res(ndof, 1);
  
  for(int i=0;i<ndof;i++){
    for(int t=0;t<sz;t++){
      res.M[i] += X.get(i,t);
    }
    res.M[i] /= double(sz);
  }

  return res;
}

CMATRIX mean(CMATRIX& X){
/** 
  Compute the mean of the data in each row of X:

  mean_i = <X_i> = (1/nsampl) sum_k { X_ik }

  X (ndof x nsampl)
  where ndof - the number of degrees of freedom 
        nsampl - the number of samples on ndof-dimensional data

*/

  int ndof = X.n_rows;
  int sz = X.n_cols;

  CMATRIX res(ndof, 1);
  
  for(int i=0;i<ndof;i++){
    for(int t=0;t<sz;t++){
      res.M[i] += X.get(i,t);
    }
    res.M[i] /= double(sz);
  }

  return res;
}


MATRIX deviation(MATRIX& X){
/**
  Returns the deviation of each component from its mean value
*/

  int ndof = X.n_rows;  
  int sz = X.n_cols;

  MATRIX res(ndof, sz);
  MATRIX E(ndof, sz); // Expectation value

  E = mean(X);
  
  for(int i=0;i<ndof;i++){
    for(int t=0;t<sz;t++){
      res.set(i,t,  X.get(i,t) - E.get(i));
    }
  }

  return res;
}

CMATRIX deviation(CMATRIX& X){
/**
  Returns the deviation of each component from its mean value
*/

  int ndof = X.n_rows;  
  int sz = X.n_cols;

  CMATRIX res(ndof, sz);
  CMATRIX E(ndof, sz); // Expectation value

  E = mean(X);
  
  for(int i=0;i<ndof;i++){
    for(int t=0;t<sz;t++){
      res.set(i,t,  X.get(i,t) - E.get(i));
    }
  }

  return res;
}

MATRIX variance(MATRIX& X, int opt){
/** 
  opt: 0 - population; 1 - sample;
*/

  int ndof = X.n_rows;
  int sz = X.n_cols;
  double denom = 1;
  if(opt==0){ denom = sz; }
  else if(opt==1){ denom = sz-1; }  

  MATRIX dx(ndof, sz);
  dx = deviation(X);

  MATRIX res(ndof, 1);
  for(int i=0; i<ndof; i++){
    double tmp = 0.0;

    for(int t=0; t<sz; t++){
      tmp += dx.get(i, t) * dx.get(i, t); 
    }
    res.set(i, 0,  tmp/denom );
  }// for i

  return res;

}

MATRIX variance(CMATRIX& X, int opt){
/**
  opt: 0 - population; 1 - sample;
*/

  int ndof = X.n_rows;
  int sz = X.n_cols;
  double denom = 1;
  if(opt==0){ denom = sz; }
  else if(opt==1){ denom = sz-1; }

  CMATRIX dx(ndof, sz);
  dx = deviation(X);

  MATRIX res(ndof, 1);
  for(int i=0; i<ndof; i++){
    double tmp = 0.0;

    for(int t=0; t<sz; t++){
      tmp += (std::conj(dx.get(i, t)) * dx.get(i, t)).real();
    }
    res.set(i, 0,  tmp/denom );
  }// for i

  return res;

}

MATRIX std_dev(MATRIX& X, int opt){

  int ndof = X.n_rows;
  MATRIX var(ndof, 1);  
  var = variance(X, opt);

  for(int i=0; i<ndof; i++){  var.set(i, 0,  sqrt(var.get(i, 0)) ); }
  return var;

}

MATRIX std_dev(CMATRIX& X, int opt){

  int ndof = X.n_rows;
  MATRIX var(ndof, 1);
  var = variance(X, opt);

  for(int i=0; i<ndof; i++){  var.set(i, 0,  sqrt(var.get(i, 0)) ); }
  return var;

}

MATRIX covariance(MATRIX& X){
/** 
  Compute a covariance of the data in X:

  cov_ij = <X_i * X_j> = (1/nsampl) sum_k { X_ik * X_jk }

  X (ndof x nsampl)
  where ndof - the number of degrees of freedom 
        nsampl - the number of samples on ndof-dimensional data

*/


  int ndof = X.n_rows;
  int sz = X.n_cols;

  MATRIX res(ndof, ndof);

  for(int i=0;i<ndof;i++){
    for(int j=0;j<ndof;j++){

      double tmp = 0.0;
      for(int t=0;t<sz;t++){  tmp += X.get(i,t) * X.get(j,t);  }

      res.set(i,j, tmp/double(sz) );

    }
  }

  return res;
}


MATRIX covariance(MATRIX& X, MATRIX& Y){
/** 
  Compute a covariance of the data in X and Y:

  cov_ij = <X_i * Y_j> = (1/nsampl) sum_k { X_ik * Y_jk }

  X (ndof x nsampl) and Y (ndof x nsampl)
  where ndof - the number of degrees of freedom 
        nsampl - the number of samples on ndof-dimensional data

*/

  int nx = X.n_rows;
  int ny = Y.n_rows;
  int sz = X.n_cols;

  if(Y.n_cols!=sz){
    cout<<"Error in covariance(MATRIX&, MATRIX&): the sizes of the matrices are incompatible\n";
    cout<<"X is a "<<nx<<" x "<<sz<<" matrix\n";
    cout<<"Y is a "<<ny<<" x "<<Y.n_cols<<" matrix\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  MATRIX res(nx, ny);

  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){

      double tmp = 0.0;
      for(int t=0;t<sz;t++){  tmp += X.get(i,t) * Y.get(j,t);  }

      res.set(i,j, tmp/double(sz) );

    }// for j
  }// for i

  return res;
}


CMATRIX covariance(CMATRIX& X){
/** 
  Compute a covariance of the data in X:

  cov_ij = <X_i * X_j> = (1/nsampl) sum_k { X_ik * conj(X_jk) }

  According to: https://en.wikipedia.org/wiki/Complex_random_variable

  X (ndof x nsampl)
  where ndof - the number of degrees of freedom 
        nsampl - the number of samples on ndof-dimensional data

*/


  int ndof = X.n_rows;
  int sz = X.n_cols;

  CMATRIX res(ndof, ndof);


  for(int i=0;i<ndof;i++){
    for(int j=0;j<ndof;j++){

      complex<double> tmp(0.0, 0.0);
      for(int t=0;t<sz;t++){  tmp += X.get(i,t) * std::conj(X.get(j,t));  }

      res.set(i,j, tmp/double(sz) );

    }
  }

  return res;
}


CMATRIX covariance(CMATRIX& X, CMATRIX& Y){
/** 
  Compute a covariance of the data in X and Y:

  cov_ij = <X_i * Y_j> = (1/nsampl) sum_k { X_ik * conj(Y_jk) }

  X (ndof x nsampl) and Y (ndof x nsampl)
  where ndof - the number of degrees of freedom 
        nsampl - the number of samples on ndof-dimensional data

*/

  int nx = X.n_rows;
  int ny = Y.n_rows;
  int sz = X.n_cols;

  if(Y.n_cols!=sz){
    cout<<"Error in covariance(CMATRIX&, CMATRIX&): the sizes of the matrices are incompatible\n";
    cout<<"X is a "<<nx<<" x "<<sz<<" matrix\n";
    cout<<"Y is a "<<ny<<" x "<<Y.n_cols<<" matrix\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  CMATRIX res(nx, ny);

  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){

      complex<double> tmp(0.0,0.0);
      for(int t=0;t<sz;t++){  tmp += X.get(i,t) * std::conj(Y.get(j,t));  }

      res.set(i,j, tmp/double(sz) );

    }// for j
  }// for i

  return res;
}




void sample(MATRIX& x, MATRIX& mean_x, MATRIX& sigma_x, Random& rnd){
/**
    """
    This function generates ntraj ndof-dimensional vectors sampled from a 
    normal distribution with a given mean and variance

    Args: 
        x ( MATRIX(ndof, ntraj) ): Each column of the matrix corresponds to 
            a vector of certain properties (e.g. coordinates, momenta, of all DOFs) for 
            a given trajectory (element of ensemble)
        mean_x ( MATRIX(ndof, 1) ):  The mean of the ndof-dimensional vector (component-wise)
        sigma_x ( MATRIX(ndof, 1) ): The variance width for each component
        rnd ( Random ): The random number generator object

    Returns:
        None: but changes the matrix ```x```

    """
*/
  int nr = x.n_rows;
  int nc = x.n_cols;

  for(int i=0;i<nr;i++){
    for(int j=0; j<nc; j++){
      x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() );
      }
  }

}


void sample(MATRIX* x, MATRIX& mean_x, MATRIX& sigma_x, Random& rnd){
/**
    """
    This function generates ntraj ndof-dimensional vectors sampled from a 
    normal distribution with a given mean and variance

    Args: 
        x ( MATRIX(ndof, ntraj) ): Each column of the matrix corresponds to 
            a vector of certain properties (e.g. coordinates, momenta, of all DOFs) for 
            a given trajectory (element of ensemble)
        mean_x ( MATRIX(ndof, 1) ):  The mean of the ndof-dimensional vector (component-wise)
        sigma_x ( MATRIX(ndof, 1) ): The variance width for each component
        rnd ( Random ): The random number generator object

    Returns:
        None: but changes the matrix ```x```

    """
*/
  int nr = x->n_rows;
  int nc = x->n_cols;

  for(int i=0;i<nr;i++){
    for(int j=0; j<nc; j++){
      x->set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() );
      }
  }

}



int set_random_state(vector<double>& prob, double ksi){
/**
    """
    This function implements a simple random state selection procedure. 
    Each state is selected with a given probability

    Args:
        prob ( list of N doubles ): The probabilities of all N states 
        ksi ( double ): A random number uniformly distributed in the range of (0.0, 1.0).
            It determines the outcome of this function.

    Returns:
        integer: finstate: The index of the selected state

    """
*/

  int nstates = prob.size();
  int finstate = 0;

  double left = 0.0;
  double right = 0.0;

  for(int i=0; i<nstates; i++){
    if(i==0){
      left = 0.0;
      right = prob[i];
    }
    else{
      left = right;
      right = right + prob[i];
    }
 
    if( (left<ksi) && (ksi<=right) ){  finstate = i; }

  }

  return finstate;

}



}// namespace libspecialfunctions
}// namespace liblibra



