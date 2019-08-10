/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef SPECIALFUNCTIONS_H
#define SPECIALFUNCTIONS_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

/// libspecialfunctions namespace
namespace libspecialfunctions{



using liblinalg::MATRIX;
using liblinalg::CMATRIX;
using liblinalg::MATRIX3x3;
using liblinalg::QUATERNION;

#define POW(x,y) ((x>=0)?pow(x,y):-pow(-x,y))
#define DELTA(i,j) ((i==j)?1.0:0.0)
#define SIGN(x) ((x==0.0)?(0.0):(x/fabs(x)))
#define MAX(x,y) ((x>=y)?x:y)
#define MIN(x,y) ((x<=y)?x:y)

double FAST_POW(double x,int n);

double sinh_(double);
double sin_(double);
double ERF(double);
double ERFC(double);
double gamma_lower(double s,double x); 
double Fn(int n,double t);

// Integrals of Gaussian functions
double gaussian_int(int n, double alp);
double gaussian_norm2(int n,double alp);
double gaussian_norm1(int n,double alp);
double gaussian_normalization_factor(int n,double alp);

  
double FACTORIAL(int);
double DFACTORIAL(int);
double BINOM(int,int);
void zero_array(double* X,int n);
void binomial_expansion(int n1,int n2,double x1,double x2,double* f,double* dfdx1, double* dfdx2,int is_derivs);
boost::python::list binomial_expansion(int n1,int n2,double x1,double x2,int is_derivs);

  
void LEGENDRE(int n,double x,double a,double b,double& p,double& q);
void CHEBYSHEV(int n,double x,double a,double b,double& p,double& q);
void LAGUERRE(int n,double x,double& p,double& q);
void HERMITE(int n,double x,double& p,double& q);

double Ellipe(double,double,int);
void   Ellipe2(double,double,double&,double&,double&);
void Jacobi_Elliptic(double u,double m,double tolerance,double& am, double& sn, double& cn, double& dn);
double Km(double m,double tol);
void Ellint(double m,double sinphi,double tol,double& Km_,double& value);


int randperm(int size,int of_size,vector<int>& result);

double RANDOM(double a,double b);

// === exp(x*dt) via eigendecomposition ====
MATRIX3x3 exp_(MATRIX3x3&, double);
MATRIX exp_(MATRIX&, double);
CMATRIX exp_(CMATRIX&, complex<double>);

// === exp(x*dt) via Taylor sum ====
MATRIX exp_2(MATRIX& x, double dt, int nterms, double max_tol);
MATRIX exp_2(MATRIX& x, double dt, int nterms);
MATRIX exp_2(MATRIX& x, double dt);

CMATRIX exp_2(CMATRIX& x, complex<double> dt, int nterms, double max_tol);
CMATRIX exp_2(CMATRIX& x, complex<double> dt, int nterms);
CMATRIX exp_2(CMATRIX& x, complex<double> dt);


// ==== exp(x*dt)*sinh(x*t)/(x*t) via eigendecomposition =====
MATRIX3x3 exp1_(MATRIX3x3&,double);
MATRIX exp1_(MATRIX&,double);






int merge_sort(vector< pair<int,double> >&, vector< pair<int,double> >&);
boost::python::list merge_sort(boost::python::list inp);


MATRIX mean(MATRIX& X);
CMATRIX mean(CMATRIX& X);

MATRIX deviation(MATRIX& X);
CMATRIX deviation(CMATRIX& X);

MATRIX covariance(MATRIX& X);
CMATRIX covariance(CMATRIX& X);


}// namespace libspecialfunctions
}// namespace liblibra

#endif // SPECIALFUNCTIONS_H


