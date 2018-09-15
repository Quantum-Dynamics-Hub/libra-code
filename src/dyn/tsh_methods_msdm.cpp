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
  \file tsh_methods_msdm.cpp
  \brief The file implements the modified simplified decay of mixing (MSDM) algorithm. 
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{



CMATRIX msdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates){
    /**
    \brief Modified Simplified Decay of Mixing (MSDM) method.
    This function implements the experimental modification (by Alexey Akimov) of the simplified decay of 
    mixing algorithm for decoherence correction ( Granucci, G.; Persico, M. J. Chem. Phys. 2007, 126, 134114)
    
    \param[in]       Coeff [ CMATRIX ] An object containig electronic DOFs. 
    \param[in]          dt [ float ] The integration timestep. Units = a.u. of time
    \param[in]      act_st [ integer ] The active state index
    \param[in]      decoh_rates [ MATRIX ] The matrix of decoherence (pure dephasing) rates between all pairs of states

    The function returns:
    # C [ CMATRIX ] - the updated state of the electronic DOF, in the same data type as the input

    */

    double kb = 3.166811429e-6;  // Hartree/K
    double sclf;

    CMATRIX C(Coeff);

    // Population of the active state
    double p_aa_old = (std::conj(C.get(act_st)) * C.get(act_st)).real();

    if(p_aa_old>1.0){
      cout<<"=== Place 1 =====\n";
      cout<<"Error in CMATRIX msdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
      cout<<"active state is larger than 1: p_aa_old = "<< p_aa_old << endl;
      cout<<"C = \n"; C.show_matrix();
      cout<<"act_st = "<<act_st<<endl;
      cout<<"Coeff = \n"; Coeff.show_matrix();
      cout<<"decoh_rates = \n"; decoh_rates.show_matrix();
      cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
      exit(0);
    }

    if(p_aa_old>0.0){

      // First - update all the coefficients for the non-active states        
      int N = Coeff.n_elts;

      double inact_st_pop = 0.0; // population of the inactive states after rescaling

      for(int i=0; i<N; i++){
        if(i != act_st){
          double itau = decoh_rates.get(i, act_st); 
          sclf = exp(-dt*itau);
          C.scale(i, 0, sclf);

          inact_st_pop += (std::conj(C.get(i)) * C.get(i)).real();
        }
      }

      if(inact_st_pop>1.0){
        cout<<"=== Place 2 =====\n";
        cout<<"Error in CMATRIX msdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
        cout<<"Total population of inactive states after rescaling is larger than 1: inact_st_pop = "<<inact_st_pop<<endl;
        cout<<"C = \n"; C.show_matrix();
        cout<<"act_st = "<<act_st<<endl;
        cout<<"Coeff = \n"; Coeff.show_matrix();
        cout<<"decoh_rates = \n"; decoh_rates.show_matrix();
        cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
        exit(0);
      }

      double p_aa_new = 1.0 - inact_st_pop;

     
      if(p_aa_new<0.0){
        cout<<"=== Place 3 =====\n";
        cout<<"Error in CMATRIX msdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
        cout<<"new population of the active state "<< p_aa_new <<" is negative...\n";
        cout<<"inact_st_pop = "<<inact_st_pop<<endl;
        cout<<"p_aa_old = "<<p_aa_old<<endl;
        cout<<"C = \n"; C.show_matrix();
        cout<<"act_st = "<<act_st<<endl;
        cout<<"Coeff = \n"; Coeff.show_matrix();
        cout<<"decoh_rates = \n"; decoh_rates.show_matrix();     
        cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
        exit(0);
      }

      sclf = sqrt( p_aa_new / p_aa_old );  // scaling factor for the active state
        
      // Rescale the active state
      C.scale(act_st, 0, sclf);

    }// if p_aa_old > 0.0

    double new_norm = (C.H() * C).get(0,0).real();
    //cout<<"new_norm = "<<new_norm<<endl;


    if(fabs(new_norm-1.0)>0.1){
      cout<<"=== Place 4 =====\n";
      cout<<"Error in CMATRIX msdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
    //  cout<<"new population of the active state "<< p_aa_new <<" is negative...\n";
    //  cout<<"inact_st_pop = "<<inact_st_pop<<endl;
      cout<<"p_aa_old = "<<p_aa_old<<endl;
      cout<<"C = \n"; C.show_matrix();
      cout<<"act_st = "<<act_st<<endl;
      cout<<"Coeff = \n"; Coeff.show_matrix();
      cout<<"decoh_rates = \n"; decoh_rates.show_matrix();     
      cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
      exit(0);
    }

        
    return C;

}


Electronic msdm(Electronic& Coeff, double dt, int act_st, MATRIX& decoh_rates){
    /**
    \brief Modified Simplified Decay of Mixing (MSDM) method.
    This function implements the experimental modification (by Alexey Akimov) of the simplified decay of 
    mixing algorithm for decoherence correction ( Granucci, G.; Persico, M. J. Chem. Phys. 2007, 126, 134114)
    
    \param[in]       Coeff [ Electronic ] An object containig electronic DOFs. 
    \param[in]          dt [ float ] The integration timestep. Units = a.u. of time
    \param[in]      act_st [ integer ] The active state index
    \param[in]      decoh_rates [ MATRIX ] The matrix of decoherence (pure dephasing) rates between all pairs of states

    The function returns:
    # C [ Electronic ] - the updated state of the electronic DOF, in the same data type as the input

    */

    double kb = 3.166811429e-6;  // Hartree/K
    double sclf;

    Electronic C(Coeff);

    // First - update all the coefficients for the non-active states        
    int N = C.nstates;
    double new_norm = 0.0;
    for(int i = 0; i<N; i++){
      if(i != act_st){
        double itau = decoh_rates.get(i, act_st); 
        sclf = exp(-dt*itau);

        C.q[i] = C.q[i] * sclf;
        C.p[i] = C.p[i] * sclf;

        new_norm += C.rho(i, i).real();

      }
    }

    // new_norm now contains the total population of all inactive states after rescaling
    // How much of population is left for the new active state
    double p_aa_new = 1.0 - new_norm;
    double p_aa_old = C.rho(act_st,act_st).real();

    sclf = 1.0;
    if(p_aa_old > 0.0){
      sclf = sqrt( p_aa_new / p_aa_old );  // scaling factor for the active state
    }
        
    // Rescale the active state
    C.q[act_st] = C.q[act_st] * sclf;
    C.p[act_st] = C.p[act_st] * sclf;

    return C;

}




}// namespace libdyn
}// liblibra

