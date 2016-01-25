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
/**
  \file Grid_functions.h
  \brief The file describes some auxiliary functions for grid operations
    
*/

#ifndef GRID_FUNCTIONS_H
#define GRID_FUNCTIONS_H

#include "../../mmath/libmmath.h"
using namespace libmmath;

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid{


//--------------------- General ---------------------

int find_grid_size(double xmin,double xmax, double dx);
CMATRIX init_grid(double xmin,double xmax, double dx);


//------------------- 1D specific -------------------
void init_gauss_1D(vector<CMATRIX>& wfc, CMATRIX& X, double x_, double p_, double dx, int nstates, int occ_state);

void print_1D(CMATRIX& X, vector<CMATRIX>& PSI,string filename);
void print_1D(CMATRIX& X, vector<CMATRIX>& PSI,string prefix, int frame);

void ft_1D(vector<CMATRIX>& psi,vector<CMATRIX>& reci_psi,int opt,double xmin,double kxmin,double dx);


//------------------ 2D specific --------------------
void init_gauss_2D(vector<CMATRIX>& wfc, 
                   CMATRIX& X,double x_,double px_,double dx, 
                   CMATRIX& Y,double y_,double py_,double dy, 
                   int nstates, int occ_state);

void print_2D(CMATRIX& X, CMATRIX& Y, vector<CMATRIX>& PSI, string filename);
void print_2D(CMATRIX& X, CMATRIX& Y, vector<CMATRIX>& PSI, string prefix, int frame);

void ft_2D(vector<CMATRIX>& psi, vector<CMATRIX>& reci_psi,int opt,
           double xmin,double kxmin,double dx,double ymin,double kymin,double dy);




}// namespace libwfcgrid
}// namespace libdyn

#endif // GRID_FUNCTIONS_H

