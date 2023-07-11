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
/**
  \file ivr_timecorr.cpp
  \brief These are the C++ re-implementation of the Fortran codes from the Ananth group.
         The original codes can be found here:
         https://github.com/AnanthGroup/SC-IVR-Code-Package    

  According to original documentation:

! This program computes the multidim. Time Correlation Function (TCF)
  
! The output file  'TCF.out'  contains the real and imaginary
! parts of the correlation function as a function of time.

  Notes on the reimplementation: this is no longer a standalone program, but
  rather a module. The results are returned in a digital form (and can pe easily
  printed out). The propagation must have been done outside the functions, 
  so we assume all the trajectories and the corresponding dependent quantities
  have been pre-computed and are now passed into these routines.

*/

#include "ivr.h"


/// liblibra namespace
namespace liblibra{

/// libivr namespace
namespace libivr{


void compute_tcf(vector< complex<double> >& TCF, vector<int>& MCnum,
                 vector<MATRIX>& q, vector<MATRIX>& p, vector<int>& status,
                 int ivr_opt, int observable_type, int observable_label){
/**

  \brief A generic routine to compute unnormalized TCF

  This function adds a contribution to the TCF due a given Monte Carlo configuration
  Call the function with the same TCF variable many times to get a convergence w.r.t.
  the Monte Carlo sampling

  \param[in/out] TCF - the storage for the unnormalized TCF being computed
  \param[in/out] MCnum - a counter of the successful trajectories at every time step
  \param[in] q - coordinates for all the propagated timesteps (even if the run was stopped earlier)
  \param[in] p - momenta for all the propagated timesteps (even if the run was stopped earlier)
  \param[in] status - the vector describing the status of all the time-points in the supplied
             trajectory. This is a list of 1 (success) and 0 (bad trajectory at that point). If we have
             status[t] == 0, then the data points q[t] and p[t] are disregarded and not used in the 
             calculations of TCF. 
  \param[in] ivr_opt - option for slection of the type of the IVR: 0 (Husimi), 1 (LS-IVR), etc.
  \param[in] observable_type - a selector of the observable type: 0 - coordinates, 1 - momenta
  \param[in] observable_label - a selector of a specific DOF 


*/

  int Ntime = TCF.size();

  for(int i = 0; i < Ntime; i++){    MCnum[i] += status[i];   }


  if(ivr_opt==0){  /// Husimi

    for(int i = 0; i < Ntime; i++){  
      if(status[i]==1){

        TCF[i] = TCF[i] + mat_elt_HUS_B(q[i], p[i], observable_type, observable_label);

      }
    } // for i

  }/// Husimi


  else if(ivr_opt==1){  /// LS-IVR

    for(int i = 0; i < Ntime; i++){  
      if(status[i]==1){

        TCF[i] = TCF[i] + mat_elt_LSC_B(q[i], p[i], observable_type, observable_label);

      }
    } // for i

  }/// LS-IVR




}




void compute_tcf(vector< complex<double> >& TCF, vector<int>& MCnum,  ivr_params& prms,
                 vector<MATRIX>& q,  vector<MATRIX>& p,  vector<int>& status, vector<double>& action,vector< vector<MATRIX> >& Mono,
                 vector<MATRIX>& qp, vector<MATRIX>& pp, vector<int>& statusp,vector<double>& actionp,vector< vector<MATRIX> >& Monop,
                 int ivr_opt, int observable_type, int observable_label){
/**

  \brief A generic routine to compute unnormalized TCF

  This function adds a contribution to the TCF due a given Monte Carlo configuration
  Call the function with the same TCF variable many times to get a convergence w.r.t.
  the Monte Carlo sampling

  \param[in/out] TCF - the storage for the unnormalized TCF being computed
  \param[in/out] MCnum - a counter of the successful trajectories at every time step

  \param[in] q - coordinates for all the propagated timesteps (even if the run was stopped earlier)
  \param[in] p - momenta for all the propagated timesteps (even if the run was stopped earlier)
  \param[in] status - the vector describing the status of all the time-points in the supplied
             trajectory. This is a list of 1 (success) and 0 (bad trajectory at that point). If we have
             status[t] == 0, then the data points q[t] and p[t] are disregarded and not used in the 
             calculations of TCF. 
  \param[in] action - action computed for all the time points of the trajectory
  \param[in] Mono - monodromy matrices sampled along the trajectory
  \param[in] qp - same meaning as q, just another point in the pair of the initial starting phase space points
  \param[in] pp - same meaning as p, just another point in the pair of the initial starting phase space points
  \param[in] statusp - same meaning as status, but for the "primed" trajectories (qp, pp)
  \param[in] actionp - action computed for all the time points of the trajectory, just another point ---
  \param[in] Monop - monodromy matrices sampled along the trajectory, just another point ---
  \param[in] ivr_opt - option for slection of the type of the IVR: 0 (FF_MQC),  1 (DHK)
  \param[in] observable_type - a selector of the observable type: 0 - coordinates, 1 - momenta
  \param[in] observable_label - a selector of a specific DOF 


*/

  int i;
  int Ndof = q[0].n_rows;
  int Ntime = TCF.size();

  for(i = 0; i < Ntime; i++){   MCnum[i] += status[i] * statusp[i];     }


  MATRIX qIn(prms.get_qIn());
  MATRIX pIn(prms.get_pIn());
  MATRIX Width0(prms.get_Width0()); 
  MATRIX invWidth0(prms.get_invWidth0());
  MATRIX WidthT(prms.get_WidthT()); 
  MATRIX invWidthT(prms.get_invWidthT());
  MATRIX TuningQ(prms.get_TuningQ());
  MATRIX invTuningQ(prms.get_invTuningQ());
  MATRIX TuningP(prms.get_TuningP());
  MATRIX invTuningP(prms.get_invTuningP());



  if(ivr_opt==0){  /// FF_MQC


    MATRIX norm(Ndof, Ndof); norm  = invTuningQ * invTuningP;
    double normC = 1.0;
    for(i = 0; i<Ndof; i++){   normC = normC * norm.get(i,i);   }
    normC = sqrt(normC);


    MATRIX qav(Ndof, 1); qav = 0.5*(q[0] + qp[0]);
    MATRIX pav(Ndof, 1); pav = 0.5*(p[0] + pp[0]);

    complex<double> ovlp = CS_overlap(qav, pav, qIn, pIn, Width0, invWidth0);
    double sampling = (std::conj(ovlp) * ovlp).real();


    // INITIAL COHERENT STATE OVERLAPS
    ovlp  = CS_overlap(q[0],  p[0],  qIn, pIn, Width0, invWidth0);
    complex<double> ovlpp = CS_overlap(qp[0], pp[0], qIn, pIn, Width0, invWidth0);

    // OVERLAP RATIO
    complex<double> OverlapR = (std::conj(ovlpp) * ovlp)/sampling;


    int Maslov  = 0;
    complex<double> prev(1.0, 0.0);

    for(i = 0; i < Ntime; i++){  

      if(status[i]==1 && statusp[i]==1){


        // POSITION MATRIX ELEMENT AT TIME t
        complex<double> posn = mat_elt_FF_B(q[i], p[i], qp[i], pp[i], WidthT, invWidthT, observable_type, observable_label);

        // MONODROMY MATRICES FOR PREFACTOR
        vector<CMATRIX> Mfwd(4);  
        Mfwd[0] = CMATRIX(Mono[i][0]);
        Mfwd[1] = CMATRIX(Mono[i][1]);
        Mfwd[2] = CMATRIX(Mono[i][2]);
        Mfwd[3] = CMATRIX(Mono[i][3]);

        vector<CMATRIX> Mbck(4);  
        Mbck[0] = CMATRIX(Mono[i][3]).T();
        Mbck[1] =-1.0 * CMATRIX(Mono[i][1]).T();
        Mbck[2] =-1.0 * CMATRIX(Mono[i][2]).T();
        Mbck[3] = CMATRIX(Mono[i][0]).T();

        // CALCULATE PREFACTOR
        complex<double> pref;
        pref = MQC_prefactor_FF_G( Mfwd, Mbck, prms);


        // TRACK MASLOV INDEX
        if ((std::abs(pref)<0.0) && (pref.imag()*prev.imag()<0.0) ) { Maslov = Maslov + 1; }
        prev  = pref;

        // CALCULATE TCF
        double argg = action[i] - actionp[i];
        TCF[i] = TCF[i] + normC * OverlapR * pow(-1.0, Maslov) * posn * sqrt(pref) * complex<double>(cos(argg), sin(argg));

      }
    } // for i

  }/// FF_QMC


  if(ivr_opt==1){  /// DHK

    vector<int> Maslov(2, 0);
    vector<complex<double> > pref(2, complex<double>(1.0, 0.0)); 
    vector<complex<double> > prev(2, complex<double>(1.0, 0.0)); 


    // INITIAL COHERENT STATE OVERLAPS
    complex<double> ovlp  = CS_overlap(q[0],  p[0],  qIn, pIn, Width0, invWidth0);
    complex<double> ovlpp = CS_overlap(qp[0], pp[0], qIn, pIn, Width0, invWidth0);

    // OVERLAP RATIO
    complex<double> OverlapR = 1.0/(std::conj(ovlpp) * ovlp);



    for(i = 0; i < Ntime; i++){  

      if(status[i]==1 && statusp[i]==1){

        // POSITION MATRIX ELEMENT AT TIME t
        complex<double> posn = mat_elt_FF_B(q[i], p[i], qp[i], pp[i], WidthT, invWidthT, observable_type, observable_label);

        // MONODROMY MATRICES FOR PREFACTOR
        vector<CMATRIX> Mfwd(4);  
        Mfwd[0] = CMATRIX(Mono[i][0]);
        Mfwd[1] = CMATRIX(Mono[i][1]);
        Mfwd[2] = CMATRIX(Mono[i][2]);
        Mfwd[3] = CMATRIX(Mono[i][3]);

        vector<CMATRIX> Mbck(4);  
        Mbck[0] = CMATRIX(Mono[i][3]).T();
        Mbck[1] =-1.0 * CMATRIX(Mono[i][1]).T();
        Mbck[2] =-1.0 * CMATRIX(Mono[i][2]).T();
        Mbck[3] = CMATRIX(Mono[i][0]).T();

        // CALCULATE PREFACTOR
        pref = DHK_prefactor(Mfwd, Mbck, prms);

        // TRACK MASLOV INDEX
        if ((std::abs(pref[0])<0.0) && (pref[0].imag()*prev[0].imag()<0.0) ) { Maslov[0] = Maslov[0] + 1; }
        if ((std::abs(pref[1])<0.0) && (pref[1].imag()*prev[1].imag()<0.0) ) { Maslov[1] = Maslov[1] + 1; }
        prev  = pref;


        // CALCULATE TCF
        double argg = action[i] - actionp[i];
        TCF[i] = TCF[i] + OverlapR * pow(-1.0, Maslov[0]+Maslov[1]) * posn * sqrt(pref[0]*pref[1]) * complex<double>(cos(argg), sin(argg));


      }// if trajectories are included

    }// for i

  }// DHK



}


void normalize_tcf(vector< complex<double> >& TCF, vector<int>& MCnum){
/**
  \brief A generic routine to normalized TCF
  Call this function once you have had enough Monte Carlo sampling events executed
     
  \param[in/out] TCF - the initially raw TCF, vector of the length Ntime, where Ntime is the 
  length of the trajectory. Once the functions is executed, the values in this  variable will be 
  replaced by the normalized values
  \param[in] MCnum - the number of successful (meaning the real-time propagation was within the
  specified tolerance level of energy and symplecticity conservation) MC events (trajectories)
  survived at each given time. So, MCnum is a vector of length Ntime, each element MCnum[t] gives
  the number of acceptable trajectories at the time t

*/

  if(TCF.size()!=MCnum.size()){
    cout<<"Error in normalize_tcf: The number of elements in TCF ("<<TCF.size()<<") is not consistent\
    with the number of elements in MCnum ("<<MCnum.size()<<")\n";
    exit(0);
  }

  int Ntime = TCF.size();

   for(int j = 0; j < Ntime; j++){

     if(MCnum[j]==0){  TCF[j] = complex<double>(0.0, 0.0); }
     else{ TCF[j] = TCF[j] * (1.0 / float(MCnum[j])); }

   }

}


}/// namespace libivr
}/// liblibra

