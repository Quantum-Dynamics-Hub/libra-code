/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef IVR_H
#define IVR_H

#include "../math_linalg/liblinalg.h"
#include "../math_random/librandom.h"
#include "../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace librandom;
using namespace libmeigen;


namespace libivr{

///============ Sampling Initial Distributions ===========================
///  In ivr_sampling.cpp

MATRIX ivr_Husimi(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd);
vector<MATRIX> ivr_Husimi(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size);

MATRIX ivr_LSC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd);
vector<MATRIX> ivr_LSC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size);

MATRIX ivr_DHK(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd);
vector<MATRIX> ivr_DHK(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size);

MATRIX ivr_FB_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd);
vector<MATRIX> ivr_FB_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd, int sample_size);

MATRIX ivr_FF_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd);
vector<MATRIX> ivr_FF_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd, int sample_size);


///============ Matrix Elements and Overlaps ===========================
///  In ivr_matrix_elements.cpp

complex<double> CS_overlap(MATRIX& q, MATRIX& p, MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& invWidth0);
complex<double> mat_elt_FB_B(MATRIX& q, MATRIX& p, int opt, int lab);
complex<double> mat_elt_FF_B(MATRIX& q, MATRIX& p, MATRIX& qp, MATRIX& pp, MATRIX& WidthT, MATRIX& invWidthT, int opt, int lab);
complex<double> mat_elt_HUS_B(MATRIX& q, MATRIX& p, int opt, int lab);
complex<double> mat_elt_LSC_B(MATRIX& q, MATRIX& p, int opt, int lab);


///============ SC Prefactors  ===========================
///  In ivr_prefactors.cpp

complex<double> MQC_prefactor_FB_G
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
 CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
 CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
);

complex<double> MQC_prefactor_FF_G
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
 CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
 CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
);

vector<complex<double> > DHK_prefactor
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
 CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
 CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
);


///============ Propagators  ===========================
///  In ivr_propagators.cpp

void Integrator(MATRIX& q, MATRIX& p, vector<CMATRIX>& M, double& action, MATRIX& mass, double dt);


}// namespace libivr
}// liblibra

#endif // IVR_H
