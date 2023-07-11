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

#ifndef POTENTIALS_MB_ELEC_H
#define POTENTIALS_MB_ELEC_H


#include "../math_linalg/liblinalg.h"
#include "../math_specialfunctions/libspecialfunctions.h"
#include "../cell/libcell.h"
#include "Switching_functions.h"
#include "Potentials_elec.h"
#include "../Units.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libcell;
using namespace libspecialfunctions;

namespace libpot{


double Elec_Ewald3D(vector<VECTOR*>& r, vector<double>& q, MATRIX3x3& box,    /* Inputs */ 
                    vector<VECTOR*>& f, MATRIX3x3& at_stress,  /* Outputs*/
                    vector< std::pair<int, int> >& exclusions, 
                    vector< VECTOR* >& excl_vector1,
                    vector< VECTOR* >& excl_vector2,
                    vector<double>& excl_scales,
                    int rec_deg,int pbc_deg, double etha, double R_on, double R_off    /* Parameters */
                   );

double Elec_Ewald3D(vector<VECTOR>& r, vector<double>& q, MATRIX3x3& box, double epsilon,   /* Inputs */ 
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


}//namespace libpot
}// liblibra

#endif //POTENTIALS_MB_ELEC_H