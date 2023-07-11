/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Potentials_mb_vdw.h
  This file implements the many-body potentials (with few-body potentials as sepcial case) involving
  vdW interactions. Something like the lattice sums, or just summing all the 2-body pairs without creating
  large number of auxiliary data.
*/


#ifndef POTENTIALS_MB_VDW_H
#define POTENTIALS_MB_VDW_H

#include "../math_linalg/liblinalg.h"
#include "../cell/libcell.h"
#include "../math_specialfunctions/libspecialfunctions.h"
#include "Switching_functions.h"
#include "Potentials_vdw.h"
#include "Potentials_elec.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libcell;
using namespace libspecialfunctions;
using libcell::triple;
using libcell::quartet;
using libcell::excl_scale;


namespace libpot{

//--------------------- Many-body potentials ----------------------------------

double VdW_Ewald3D(vector<VECTOR>& r, vector<int>& types, int max_type, vector<double>& Bij, MATRIX3x3& box, /* Inputs */ 
                   vector<VECTOR>& f, MATRIX3x3& at_stress,  /* Outputs*/
                   int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */                   
                   );

double VdW_Ewald3D(vector<VECTOR>& r, vector<double>& q, MATRIX3x3& box, /* Inputs */ 
                   vector<VECTOR>& f, MATRIX3x3& at_stress,  /* Outputs*/
                   int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */
                   );


double Elec_Ewald3D(VECTOR* r,         /* Inputs */
                    VECTOR* g,
                    VECTOR* m,
                    VECTOR* f,
                    MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs */
                    int sz,double* q,
                    int nexcl, int* excl1, int* excl2, double* scale,
                    MATRIX3x3* box,int rec_deg,int pbc_deg,
                    double etha,int is_cutoff, double R_on, double R_off,
                    int& time,vector< vector<triple> >& images, vector<triple>& central_translation,
                    double* dr2,double dT, int& is_update);  /* Parameters */                  

double Vdw_LJ(VECTOR* r,                                               /* Inputs */
              VECTOR* g,
              VECTOR* m,
              VECTOR* f,
              MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
              int sz,double* epsilon, double* sigma,
              int nexcl, int* excl1, int* excl2, double* scale,
              MATRIX3x3* box,int rec_deg,int pbc_deg,
              double etha,int is_cutoff, double R_on, double R_off,
              int& time, vector< vector<triple> >& images, vector<triple>& central_translation,
              double* dr2, double dT, int& is_update     /* Parameters */
             );

double Vdw_LJ1(VECTOR* r,                                               /* Inputs */
               VECTOR* g,
               VECTOR* m,
               VECTOR* f,
               MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
               int sz,double* epsilon, double* sigma,
               int nexcl, int* excl1, int* excl2, double* scale,
               MATRIX3x3* box,int rec_deg,int pbc_deg,
               double etha,int is_cutoff, double R_on, double R_off,
               int& time, vector< vector<quartet> >& images, vector<triple>& central_translation,
               double* dr2, double dT, int& is_update     /* Parameters */
             );

double Vdw_LJ2_no_excl(VECTOR* r,                                               /* Inputs */
                       VECTOR* g,
                       VECTOR* m,
                       VECTOR* f,
                       MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                       int sz,double* epsilon, double* sigma,
                       int nexcl, int* excl1, int* excl2, double* scale,
                       MATRIX3x3* box,int rec_deg,int pbc_deg,
                       double etha,int is_cutoff, double R_on, double R_off,
                       int& time,vector< vector<excl_scale> >& excl_scales
                      );

double Vdw_LJ2_excl(VECTOR* r,                                               /* Inputs */
                    VECTOR* g,
                    VECTOR* m,
                    VECTOR* f,
                    MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                    int sz,double* epsilon, double* sigma,
                    int nexcl, int* excl1, int* excl2, double* scale,
                    MATRIX3x3* box,int rec_deg,int pbc_deg,
                    double etha,int is_cutoff, double R_on, double R_off,
                    int& time,vector< vector<excl_scale> >& excl_scales
                   );


double LJ_Coulomb(VECTOR* r, VECTOR* g, VECTOR* m, VECTOR* f,
                  MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress,
                  int sz,double* epsilon, double* sigma,double* q,int is_cutoff, double R_on, double R_off,
                  int nexcl, int* excl1, int* excl2, double* scale);



} // namespace libpot
} //liblibra

#endif // POTENTIALS_MB_VDW_H

