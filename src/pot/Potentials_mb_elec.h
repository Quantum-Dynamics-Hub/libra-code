#ifndef POTENTIALS_MB_ELEC_H
#define POTENTIALS_MB_ELEC_H

#include "../mmath/libmmath.h"
#include "../cell/libcell.h"
#include "Switching_functions.h"
#include "Potentials_elec.h"

using namespace libmmath;
using namespace libcell;

namespace libpot{

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

#endif //POTENTIALS_MB_ELEC_H