/*********************************************************************************
* Copyright (C) 2024 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file matrix_gradients.cpp
  \brief The file implements various kinds of contributions to Lagrangian and
  their gradients needed for optimization

*/

#include "opt.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;

/// libopt namespace
namespace libopt{


void normalize_transition_matrix(MATRIX& flux, MATRIX& num_coeff, MATRIX& denom_coeff){
/**
   J_ij = num_coeff_ij * flux_ij / sum_k {  denom_coeff_kj * flux_kj }

*/

  int i,j,n,m;
  n = flux.n_rows; 
  m = flux.n_cols;

  for(j=0;j<m;j++){  

    double sum = 0.0;

    for(i=0;i<n;i++){  sum += denom_coeff.get(i,j)*flux.get(i,j);   }
    if( fabs(sum)>0.0){ 
      for(i=0;i<n;i++){ flux.scale(i, j, num_coeff.get(i,j)/sum ); }
    }
  }// for j

}
/*
void normalize_transition_matrix(MATRIX& flux, int opt){
  normalize_transition_matrix(flux, opt, 0);
}

void normalize_transition_matrix(MATRIX& flux){
  normalize_transition_matrix(flux, 1, 0);
}
*/
void hyperbolic(double x, double alpha, double& val, double& deriv){
/**
*/

  double arg = alpha*x;
  if(arg<-50.0){ val = 0.0; deriv = 1.0; }
  else if(arg>50.0){ val = x; deriv = 1.0; }
  else{
    double ex = exp(-arg);
    double den = (1.0 + ex);
    val = x/den;
    deriv = (1.0 + (1.0 + arg) * ex )/(den*den);
  }

}

MATRIX transform_to_hyperbolic(MATRIX& A, double alp){
  int i,j;
  int n = A.n_cols;
  MATRIX res(n,n);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      double aij = A.get(i,j);
      double val, deriv;
      hyperbolic(aij, alp, val, deriv);
      res.set(i,j, val);
    }
  }
  return res;
}


double derivs0(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA){
/**
  L0 = |y - Ax|^2 

  dL0/dA = 2(Ax-y)x^T 

*/
  int n = x.n_rows;
  MATRIX dx(n,1);
  dx = A*x - y;
  double Lagrangian = (dx.T() * dx).get(0,0);
  gradA = 2.0 * dx * x.T();

  return Lagrangian;
}

double derivs1(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, double& L, double& gradL, int opt, int opt2){
/**

   opt2 = 0

   L1 = lambda * sum_{i,j, i \neq j} A_ij^2  <-- linear          (opt = 0)
 
   L1 = lambda^2 * sum_{i,j, i \neq j} A_ij^2  <-- quadratic     (opt = 1)

   dL1/dA = 2[A - diag(A)]^T * lambda ( or lambda^2)

   dL1/dlambda =  sum_{i,j, i \neq j} A_ij^2 or                  (opt = 0)

                 2*lambda   sum_{i,j, i \neq j} A_ij^2           (opt = 1)


   opt2 = 1

   L1 = lambda * sum_{i,j} A_ij^2  <-- linear          (opt = 0)

   L1 = lambda^2 * sum_{i,j} A_ij^2  <-- quadratic     (opt = 1)

   dL1/dA = 2A^T  * lambda ( or lambda^2)

   dL1/dlambda =  sum_{i,j} A_ij^2 or                  (opt = 0)

                 2*lambda   sum_{i,j} A_ij^2           (opt = 1)
   

*/

  int i, j, k;
  int n = A.n_cols;
  double Lagrangian;

  MATRIX diagA(n,n);
  double sum = 0.0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){ if (i!=j){  sum += A.get(i,j)* A.get(i,j); }   }

    diagA.set(i,i, A.get(i,i)); 
    if(opt2==1){   for(i=0;i<n;i++){ sum += A.get(i,i)* A.get(i,i); }   }
  }

  if(opt2==0){   gradA = 2.0*(A - diagA).T(); }
  else if(opt2==1){    gradA = 2.0*A.T();  }

  if(opt==0){   Lagrangian = L * sum;  gradL = sum;  gradA *= L; }
  else if(opt==1){   Lagrangian = L*L * sum; gradL = 2.0*L*sum; gradA *= L*L; }

  return Lagrangian;

}

double derivs1_new(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, double& L, double& gradL, MATRIX& num_coeff, MATRIX& denom_coeff, int opt){
/**

   L1 = lambda * sum_ij { J_ij^2}         (opt = 0)
   L1 = lambda^2 * sum_ij { J_ij^2 }      (opt = 1)


   J_ij = a_ij * f(A_ij) / sum_k { b_kj * f(A_kj) }


   where a_ij and b_ij are the num_coeff and denom_coeff respectively



   (dL3/dA)_{b,a} = 2 * lambda * f'(A_ab)/N_b^2 * Delta_ab    (opt = 0)

   (dL3/dA)_{b,a} = 2 * lambda^2 f'(A_ab)/N_b^2 * Delta_ab    (opt = 1)


   dL3/dlambda =  sum_ij { J_ij^2  }                         (opt = 0)

               = 2*lambda * sum_ij { J_ij^2 }                (opt = 1)


  Delta_ab = sum_i { ( J_ab * a_ab * b_ib - J_ib * a_ib * b_ab )  * f(A_ib) }

*/

  int i, j, k, a, b;
  int n = A.n_cols;
  double Lagrangian = 0.0;

  double alp = 10.0;

  MATRIX f(n,n);
  MATRIX df(n,n);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      double aij = A.get(i,j);
      double val, deriv;
      hyperbolic(aij, alp, val, deriv);
      f.set(i,j, val);
      df.set(i,j, deriv);
    }
  }

  MATRIX J(f);
  normalize_transition_matrix(J, num_coeff, denom_coeff);

  double sum = 0.0;
  for(i=0;i<n;i++){ for(j=0;j<n;j++){  sum += J.get(i,j)*J.get(i,j);  }   }

  vector<double> N(n, 0.0);
  MATRIX Delta(n,n);

  for(b=0;b<n;b++){
    for(a=0;a<n;a++){ N[b] +=  denom_coeff.get(a, b) * f.get(a, b);  }
  }// for b

  for(a=0; a<n; a++){ 
    for(b=0;b<n;b++){

      double val = 0.0;
      for(i=0;i<n;i++){
        double val1 = J.get(a,b)*num_coeff.get(a,b)*denom_coeff.get(i,b);
        double val2 = J.get(i,b)*num_coeff.get(i,b)*denom_coeff.get(a,b);

        val += ( val1 - val2 )*f.get(i,b);
      }// for k

      Delta.set(a,b, val);
    }// for j
  }// for i


  gradA = 0.0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){

      double val = Delta.get(i,j) * df.get(i,j) / (N[j] * N[j]);

      if(opt==0){ gradA.set(j,i, val*L); }
      else if(opt==1){ gradA.set(j,i, val*L*L); }

    }
  }


  if(opt==0){
    Lagrangian = L * sum;
    gradL = sum;
  }
  else if(opt==1){
    Lagrangian = L * L * sum;
    gradL = 2.0*L*sum;
  }

  return Lagrangian;

}

double derivs2(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, MATRIX& mu, MATRIX& gradmu, int opt, int opt2){
/**

   opt2 = 0 - constraint on columns
   opt2 = 1 - constraint on rows (just use the transposed A)


   L2 = sum_j { mu_i * (\sum_i{A_ij} - 1)^2 }   <--- linear      (opt = 0)

   L2 = sum_j { mu_i^2 * (\sum_i{A_ij} - 1)^2 }   <--- quadratic  (opt = 1)

   (dL2/dA)_{b,a} =  2 *mu_i^2 B_{b}; (for quadratic)
   dL2/dmu_i = 2 mu_i * B_i^2 (for quadratic)

*/
    int i, j, k;
    int n = A.n_cols;
    double Lagrangian = 0.0;

    MATRIX a(n,n);

    if(opt2==0){  a = A; }
    else if(opt2==1){  a = A.T(); }

    vector<double> b(n, 0.0);
    for(i=0;i<n;i++){  b[i] =  a.sum_col(i) - 1.0; }

    gradmu = 0.0; 
    gradA = 0.0;

    for(i=0;i<n;i++){
      double b2 =  b[i] * b[i];
      double mu2 = mu.get(i, 0) *  mu.get(i, 0);

      if(opt==0){  
        Lagrangian += mu.get(i, 0) * b2; 
        gradmu.set(i,0, b2); 

        for(j=0;j<n;j++){ gradA.set(i,j, 2.0* mu.get(i, 0) * b[i]); }

      }
      else if(opt==1){  
        Lagrangian +=  mu2 * b2; 
        gradmu.set(i,0, 2.0*mu.get(i, 0)*b2);

        for(j=0;j<n;j++){ gradA.set(i,j, 2.0* mu2 * b[i]); }
      }
    }// for i

  return Lagrangian;

}


double derivs3(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, double& L, double& gradL, int opt, int opt2, MATRIX& s){
/**
   opt = 0 -linear in lambda,
   opt = 1 - quadratic in lambda


   opt2 = 0 - no off-diagonal elements


   L3 = lambda * sum_{i,j, i \neq j} s_ij^2 * J_ij^2  <-- linear          (opt = 0)

   L3 = lambda^2 * sum_{i,j, i \neq j} s_ij^2 * J_ij^2  <-- quadratic     (opt = 1)


   J_ij = s_ij * f(A_ij) / N_j  where N_j = sum_{k} { f(A_kj) }

   f(x) = x/(1 + exp(-alpha*x)) 



   (dL3/dA)_{b,a} = 2 * lambda * s_ab^2 * f'(A_ab)/N_b * Delta_ab    (opt = 0)

   (dL3/dA)_{b,a} = 2 * lambda^2 * s_ab^2 * f'(A_ab)/N_b * Delta_ab  (opt = 1)

   where I_b = sum_{i\neq b} { J_ib^2 }

   dL3/dlambda =  sum_{i,j, i \neq j} s_ij^2 * J_ij^2 or                  (opt = 0)

               = 2*lambda * sum_{i,j, i \neq j} s_ij^2 * J_ij^2           (opt = 1)


  Delta_ab = sum_i { ( s_ab^2 * (1 - delta_ab) * J_ab - s_ib^2 * (1 - delta_ib) * J_ib  )  * J_ib} 


   opt2 = 1 - with diagonal elements (not yet developed)

*/

  int i, j, k;
  int n = A.n_cols;
  double Lagrangian = 0.0;

  double alp = 10.0;

  MATRIX f(n,n);
  MATRIX df(n,n);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      double aij = A.get(i,j);
      double val, deriv;
      hyperbolic(aij, alp, val, deriv);
      f.set(i,j, val);
      df.set(i,j, deriv);
    }
  }

  MATRIX J(f);
//  normalize_transition_matrix(J,1); SHOULD NORMALIZE IT!

  double sum = 0.0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      double sij = s.get(i,j);
      if(opt2==0){ 
        if(i!=j){ sum += J.get(i,j)*J.get(i,j) * sij * sij; }
      }
      else{ sum += J.get(i,j)*J.get(i,j) * sij * sij; }
    }
  }

  vector<double> I(n, 0.0);
  vector<double> N(n, 0.0);
  MATRIX Delta(n,n);

  for(i=0;i<n;i++){
      N[i] = f.sum_col(i);
//      I[i] = J.sum_col(i, 2) - J.get(i,i)*J.get(i,i);

    for(j=0;j<n;j++){ 
      double val = 0.0;
      for(k=0;k<n;k++){
        double sij = s.get(i,j);
        double skj = s.get(k,j);
        double val1, val2;
        if(opt2==0){
          val1 = 0.0;  if(i==j){ val1 = 0.0; }else{ val1 = J.get(i,j) * sij * sij; }
          val2 = 0.0;  if(k==j){ val2 = 0.0; }else{ val2 = J.get(k,j) * skj * skj; }
        }
        else{  
          val1 = J.get(i,j) * sij * sij; 
          val2 = J.get(k,j) * skj * skj;
        }
        val += (val1 - val2)*J.get(k,j);
 
      }// for k 

      Delta.set(i,j, val);
    }// for j
  }// for i

  gradA = 0.0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){

      double val = Delta.get(i,j) * df.get(i,j) / N[j];

      if(opt==0){ gradA.set(j,i, val*L); }
      else if(opt==1){ gradA.set(j,i, val*L*L); }

    }
  }

  if(opt==0){  
    Lagrangian += L * sum;
    gradL = sum;
  }
  else if(opt==1){ 
    Lagrangian += L*L * sum; 
    gradL = 2.0*L*sum;
  }

  return Lagrangian;
}


double derivs4(MATRIX& x, MATRIX& y, MATRIX& A, MATRIX& gradA, MATRIX& mu, MATRIX& gradmu, int opt){
/**


   L4 = sum_i { mu_i * (\sum_j{J_ij} - 1)^2 }   <--- linear      (opt = 0)

   L4 = sum_i { mu_i^2 * (\sum_j{J_ij} - 1)^2 }   <--- quadratic  (opt = 1)
   

   b_i = sum_j {J_ij} - 1

   (dL4/dA)_{b,a} =  2 f'(A_ab)/N_b * sum_i { ( mu_a * b_a - mu_i * b_i ) * J_ib } ( opt = 0)

   dL4/dmu_i = b_i^2 (opt = 0)


   (dL4/dA)_{b,a} =  2 f'(A_ab)/N_b * sum_i { ( mu_a**2 * b_a - mu_i**2 * b_i ) * J_ib } ( opt = 1)

   dL4/dmu_i = 2.0 * mu_i * b_i^2 (opt = 1)

*/
  int i, j, k;
  int n = A.n_cols;

  MATRIX f(n,n);
  MATRIX df(n,n);

  double alp = 10.0; 

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      double aij = A.get(i,j);
      double val, deriv;
      hyperbolic(aij, alp, val, deriv);
      f.set(i,j, val);
      df.set(i,j, deriv);
    }
  }

  MATRIX J(f);
//  normalize_transition_matrix(J,1); SHOULD NORMALIZE

  vector<double> b(n, 0.0);
  for(i=0;i<n;i++){  b[i] =  J.sum_row(i) - 1.0; }

  gradmu = 0.0;
  gradA = 0.0;
  double Lagrangian = 0.0;

  for(i=0;i<n;i++){
    double b2 =  b[i] * b[i];
    double mu2 = mu.get(i, 0) *  mu.get(i, 0);

    if(opt==0){
      Lagrangian += mu.get(i, 0) * b2;
      gradmu.set(i,0, b2);

      for(j=0;j<n;j++){ 
        double val = 0.0;

        for(k=0;k<n;k++){ val += (mu.get(i,0) * b[i] - mu.get(k,0) * b[k]) * J.get(k,j);  }

        gradA.set(j, i, val);
      }// for j
    }// opt = 0

    else if(opt==1){
      Lagrangian +=  mu2 * b2;
      gradmu.set(i,0, 2.0*mu.get(i, 0)*b2);

      for(j=0;j<n;j++){ 
        double val = 0.0;

        for(k=0;k<n;k++){ val += (mu.get(i,0) * mu.get(i,0) * b[i] - mu.get(k,0) * mu.get(k,0) * b[k]) * J.get(k,j);  }

        gradA.set(j, i, val); 
      }// for j
    }// opt = 1
  }// for i

  return Lagrangian;
  
}




}// namespace libopt
}// liblibra

