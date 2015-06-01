#ifndef SPECIALFUNCTIONS_H
#define SPECIALFUNCTIONS_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>



//#include "MATRIX.h"
//#include "MATRIX3x3.h"
//#include "QUATERNION.h"
#include "../linalg/liblinalg.h"


//using libmmath::liblinalg;


namespace libmmath{
namespace libspecialfunctions{

//using namespace libmmath::liblinalg;


using libmmath::liblinalg::MATRIX;
using libmmath::liblinalg::MATRIX3x3;
using libmmath::liblinalg::QUATERNION;

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

MATRIX exp_(MATRIX& ,double);
MATRIX exp1_(MATRIX&,double);
MATRIX3x3 exp_(MATRIX3x3& ,double);
MATRIX3x3 exp1_(MATRIX3x3&,double);


void MATRIX_TO_QUATERNION(MATRIX&,QUATERNION&);
void QUATERNION_TO_MATRIX(QUATERNION&,MATRIX&);
void MATRIX_TO_QUATERNION(MATRIX3x3&,QUATERNION&);
void QUATERNION_TO_MATRIX(QUATERNION&,MATRIX3x3&);
void solve_linsys(MATRIX&,MATRIX&, MATRIX&,double,int);

int merge_sort(vector< pair<int,double> >&, vector< pair<int,double> >&);


}// namespace libspecialfunctions
}// namespace libmmath

#endif // SPECIALFUNCTIONS_H


