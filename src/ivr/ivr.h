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


class ivr_params{

public:

  int Ndof;             ///< The number of DOFs

  MATRIX* qIn;          ///<  Ndof x 1 matrix of center of the distributions
  MATRIX* pIn;          ///<  Ndof x 1 matrix of center of the distributions

  MATRIX* Width0;       ///<  Ndof x Ndof matrix of width parameters at t = 0
  MATRIX* WidthT;       ///<  Ndof x Ndof matrix of width parameters at t = T
  MATRIX* invWidth0;    ///<  Ndof x Ndof matrix of inverse width parameters at t = 0
  MATRIX* invWidthT;    ///<  Ndof x Ndof matrix of inverse width parameters at t = T
                       

  MATRIX* TuningQ;      ///<  Ndof x Ndof matrix of coordinate tuning parameters
  MATRIX* TuningP;      ///<  Ndof x Ndof matrix of momentum tuning parameters
  MATRIX* invTuningQ;   ///<  Ndof x Ndof inverse matrix of coordinate tuning parameters
  MATRIX* invTuningP;   ///<  Ndof x Ndof inverse matrix of momentum tuning parameters

  ivr_params(int ndof_);               ///< Constructor
  ivr_params(const ivr_params& par_);  ///< copy Constructor
  ~ivr_params();                       ///< destructor

  void set_qIn(MATRIX& qIn_);
  void set_pIn(MATRIX& pIn_);
  void set_Width0(double w);
  void set_WidthT(double w);
  void set_TuningQ(double cq);
  void set_TuningP(double cp);

  MATRIX get_qIn(){ return *qIn; }
  MATRIX get_pIn(){ return *pIn; }
  MATRIX get_Width0(){ return *Width0; }
  MATRIX get_WidthT(){ return *WidthT; }
  MATRIX get_invWidth0(){ return *invWidth0; }
  MATRIX get_invWidthT(){ return *invWidthT; }
  MATRIX get_TuningQ(){ return *TuningQ; }
  MATRIX get_TuningP(){ return *TuningP; }
  MATRIX get_invTuningQ(){ return *invTuningQ; }
  MATRIX get_invTuningP(){ return *invTuningP; }




};



class ivr_observable{

public:
  int observable_type;    // 0 - for q, 1 - for p
  int observable_label;   // index of the DOF

  ivr_observable(){ observable_type = 0; observable_label = 0; }
  ivr_observable(int obs_typ_, int obs_lab_){ observable_type = obs_typ_; observable_label = obs_lab_; }

};


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
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, ivr_params& prms);

complex<double> MQC_prefactor_FF_G
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, ivr_params& prms);


vector<complex<double> > DHK_prefactor
(vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, ivr_params& prms);


///============ Propagators  ===========================
///  In ivr_propagators.cpp
void Integrator(MATRIX& q, MATRIX& p, vector<MATRIX>& M, double& action, MATRIX& mass, double dt);



///============ TCF calculators  ===========================
///  In ivr_timecorr.cpp

void compute_tcf(vector< complex<double> >& TCF, vector<int>& MCnum,
                 vector<MATRIX>& q, vector<MATRIX>& p, vector<int>& status,
                 int ivr_opt, int observable_type, int observable_label);


void compute_tcf(vector< complex<double> >& TCF, vector<int>& MCnum,
                 ivr_params& prms,
                 vector<MATRIX>& q,  vector<MATRIX>& p,  vector<int>& status, vector<double>& action,vector< vector<MATRIX> >& Mono,
                 vector<MATRIX>& qp, vector<MATRIX>& pp, vector<int>& statusp,vector<double>& actionp,vector< vector<MATRIX> >& Monop,
                 int ivr_opt, int observable_type, int observable_label);


void normalize_tcf(vector< complex<double> >& TCF, vector<int>& MCnum);

}// namespace libivr
}// liblibra

#endif // IVR_H
