/*********************************************************************************
* Copyright (C) 2017-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file FT.h
  \brief The file describes a set of functions for complex Fourier transforms and convolution
    
*/


#ifndef FT_H
#define FT_H

#include "CMATRIX.h"



/// liblibra 
namespace liblibra{

using namespace std;

/// liblinalg namespace
namespace liblinalg{


void solve_linsys(CMATRIX& C,CMATRIX& D, CMATRIX& X,double eps,int maxiter,double omega);
void solve_linsys1(CMATRIX& C,CMATRIX& X,double eps,int maxiter,double omega);


//---------- Fourier transforms ----------------
void dft(CMATRIX& in,CMATRIX& out);
void inv_dft(CMATRIX& in,CMATRIX& out);

void cft(CMATRIX& in,CMATRIX& out,double xmin,double dx);
void inv_cft(CMATRIX& in,CMATRIX& out,double xmin,double dx);

void cft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);
void inv_cft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);

void cft2(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk);
void inv_cft2(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk);

void cft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);
void inv_cft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);

void convolve(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx);
void convolve_2D(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx,double dy);

//-------- Fast Fourier Transforms -------------
void cfft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);  
void inv_cfft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx);

void cfft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);
void inv_cfft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy);



}//namespace liblinalg
}// namespace liblibra

#endif // FT_H
