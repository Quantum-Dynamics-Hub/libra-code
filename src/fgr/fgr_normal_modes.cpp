/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This code is partially based on the code of Xiang Sun:
* Copyright (C) Xiang Sun 2015
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file fgr_normal_modes.cpp
  \brief This file implements the normal modes analysis needed in the FGR calculations
  dynamics rates. 
*/

#include "fgr.h"


/// liblibra namespace
namespace liblibra{

/// libivr namespace
namespace libfgr{



///=========spectral densities ==========

double S_omega_ohmic(double omega, double etha, double omega_c) {
/** Ohmic spectral density
 omega - frequency
 etha - friction coefficient that dictates the coupling strength between 
        primary and secondary modes
 omega_c - frequency exponential cutoff 

 Ref:  JCP 2016, 144, 244105  Eq. 21

*/
  return etha * omega * exp(- omega  / omega_c);
}

double S_omega_drude(double omega, double etha) {
/** Drude spectral density
 omega - frequency
 etha - friction coefficient that dictates the coupling strength between 
        primary and secondary modes
 omega_c - frequency exponential cutoff 

 Ref: similar to JPCA 2016, 120, 2976  Eq. 40

*/

  return etha * omega /(1.0 + omega*omega);
}

double S_omega_gaussian(double omega, double etha, double sigma, double omega_op) {
/** Optical Gaussian spectral density
 omega - frequency
 etha - friction coefficient that dictates the coupling strength between 
        primary and secondary modes
 omega_c - frequency exponential cutoff 

 Ref:  JPCA 2016, 120, 2976  Eq. 39

*/
  double argg = (omega - omega_op)/sigma;
  return   0.5 * (etha * omega/(sqrt(2.0*M_PI)*sigma)) * exp(-0.5*argg*argg);
}

/**
double J_omega_ohmic(double omega, double etaa) {
    //notice definition J(omega) is different from S(omega)
    //J_omega = pi/2 * sum_a c_a^2 / omega_a delta(omega - omega_a)
    return etaa * omega * exp(-1 * omega / omega_c);
}

double J_omega_ohmic_eff(double omega, double etaa) {
    //(normal mode) effective SD for Ohmic bath DOF
    //J_omega = pi/2 * sum_a c_a^2 / omega_a delta(omega - omega_a)
    return etaa * omega * pow(Omega,4) / ( pow(Omega*Omega - omega*omega, 2) + etaa*etaa*omega*omega);
}



*/




double eq_shift(double Er, double Omega){
/**
  Compute the (y0) - half of the shift in 
  equilibrium geometry along the primary mode 
  between the donor and acceptor states from 
  the reorganization energy and primary mode frequency

  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, Eq. (44)

  \param[in] Er        - reorganization energy
  \param[in] Omega     - the primary mode mass-weighted frequency 

*/

  double y0 = sqrt(0.5*Er/(Omega*Omega));
  return y0;

}

double reorganization_energy(double y0, double Omega){
/**
  Compute the reorganization energy from the
  y0 - half of the shift in 
  equilibrium geometry along the primary mode 
  between the donor and acceptor states from 
  and primary mode frequency

  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, Eq. (44)

  \param[in] y0        - half of the shift in eq. geom. along the primary mode 
  \param[in] Omega     - the primary mode mass-weighted frequency 

*/

  y0 = Omega*y0;
  double Er = 2.0*y0*y0;
  return Er;

}


double reorganization_energy(vector<double>& omega_nm, vector<double>& req_nm){
/**
  Compute the reorganization energy 
 
  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, ??

  \param[in] omega_nm  - normal-mode frequencies [a.u.]
  \param[in] req_nm    - normal-mode displacements upon charge transfer [a.u.]

*/
  int sz = omega_nm.size();

  if(req_nm.size()!=sz){
    cout<<"ERROR in reorganization_energy: the dimensions of the input arguments do not match\n";
    cout<<"len(omega_nm) = "<<omega_nm.size()<<endl;
    cout<<"len(req_nm) = "<<req_nm.size()<<endl;
    cout<<"Exiting...\n";
    exit(0);
  }

  double Er = 0.0;
  for(int i=0;i<sz;i++){ 
    double d = omega_nm[i]*req_nm[i];
    Er += d*d;
  }
  return 0.5*Er;

}




double diabat_crossing(double dE, double Er, double y0){
/**
  Compute the point of the diabat intersection assuming (y_A = -y_D)

  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, Eq. (47)

  \param[in] dE        - the energy difference between minima of the donor-acceptor diabats
                          dE = E_D - E_A = E_D(R_D) - E_A(R_A)
  \param[in] Er        - the reorganization energy
  \param[in] y0        - the half of the shift in equilibrium geometry along 
                         the primary mode between the donor and acceptor states

*/

  double y_c = y0 * dE/Er;
  return y_c;
}




double coupling_Condon(double gamma, double dE, double Er, double y0){
/**
  Compute the donor-acceptor coupling in a non-Condon approximation

  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, Eq. (47)

  \param[in] gamma     - the linear coupling coefficient for the primary mode
  \param[in] dE        - the energy difference between minima of the donor-acceptor diabats
                          dE = E_D - E_A = E_D(R_D) - E_A(R_A)
  \param[in] Er        - the reorganization energy
  \param[in] y0        - the half of the shift in equilibrium geometry along 
                         the primary mode between the donor and acceptor states

*/

  double Gamma = gamma*y0*(Er + dE)/Er;
  return Gamma;
}


double coupling_non_Condon(double y, double gamma, double dE, double Er, double y0){
/**
  Compute the donor-acceptor coupling in a non-Condon approximation

  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, Eq. (46)

  \param[in] y         - the primary mode mass-weighted coordinate
  \param[in] gamma     - the linear coupling coefficient for the primary mode
  \param[in] dE        - the energy difference between minima of the donor-acceptor diabats
                          dE = E_D - E_A = E_D(R_D) - E_A(R_A)
  \param[in] Er        - the reorganization energy
  \param[in] y0        - the half of the shift in equilibrium geometry along 
                         the primary mode between the donor and acceptor states

*/

  double Gamma_cond = coupling_Condon(gamma, dE, Er, y0);
  double y_c = diabat_crossing(dE, Er, y0);

  double Gamma = Gamma_cond + gamma*(y-y_c);
  return Gamma;
}





void error_msg(std::string name1, int sz1, std::string name2, int sz2){

    if(sz1!=sz2){
      cout<<"ERROR: The size of the "<<name1<<" ("<<sz1
          <<") is not consistent with the size of "<<name2<<" ("<<sz2<<")\n";
      exit(0);
    }

}


vector<double> normal_modes(vector<double>& omega, vector<double>& coeff, MATRIX& T){
/**
  This function computes the normal modes using the primary more, bath modes and bath coupling coefficients
  
  \param[in] omega     - the first element is the primary mode (Omega), other elements are the bath modes
  \param[in] coeff     - coeff[0] = 0 - corresponds to the primary mode, other are the couplings of the
                         primary mode to the bath modes
  \param[out] T        - the transformation coordinates (normal modes)

  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109

*/

    int i;

    /// Check the dimensions of the input and output
    int ndof = omega.size();
    error_msg("coeff", coeff.size(), "omega", ndof);
    error_msg("T.num_of_cols", T.n_cols, "omega", ndof);
    error_msg("T.num_of_rows", T.n_rows, "omega", ndof);


    /// Compute the D-matrix, Eq. (49)
    MATRIX D(ndof, ndof);

    D.set(0, 0, omega[0]*omega[0]);
    for(i=1; i<ndof; i++){  
      D.add(0, 0, coeff[i]*coeff[i]/(omega[i]*omega[i]) ); 
      D.set(i, i, omega[i]*omega[i] );
      D.set(0, i, coeff[i] );   
      D.set(i, 0, coeff[i] ); 
    }
   
    /// diagonalize matrix, the eigenvectors transpose is in result matrix => T:  T^+ * D * T = W 
    /// D * T = T * W
    MATRIX W(ndof, ndof);    
    solve_eigen(D, W, T, 0);  // W contains omega_nm^2 at this point

    /// Normal mode frequencies
    vector<double> omega_nm(ndof, 0.0);
    for(i=0; i<ndof; i++){  omega_nm[i] = sqrt(W.get(i,i)); }


    return omega_nm;

}


vector<double> compute_req(vector<double>& omega, vector<double>& coeff, double y0, MATRIX& T){
/**
  \param[in] omega     - the first element is the primary mode (Omega), other elements are the bath modes
  \param[in] coeff     - coeff[0] = 0 - corresponds to the primary mode, other are the couplings of the
                         primary mode to the bath modes
  \param[in] y0        - the half of the shift in equilibrium geometry along 
                         the primary mode between the donor and acceptor states
  \param[in] T         - the transformation coordinates (normal modes)

  Ref: 
  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109

*/

    int i;

    /// Check the dimensions of the input and output
    int ndof = omega.size();
    error_msg("coeff", coeff.size(), "omega", ndof);
    error_msg("T.num_of_cols", T.n_cols, "omega", ndof);
    error_msg("T.num_of_rows", T.n_rows, "omega", ndof);


    /// Equilibrium donor-to-acceptor displacement in equilibrium geometry along the normal
    /// mode coordinates - (acceptor's potential energy min shift) Eq. (51)
    MATRIX tmp(ndof, 1);
    tmp.set(0, 0, 1.0);
    for(i=1; i<ndof; i++){ tmp.set(i, 0, -coeff[i]/(omega[i]*omega[i]));    }
    tmp = 2.0*y0*T.T()*tmp;

    vector<double> req(ndof, 0.0);
    for(i=0; i<ndof; i++){ req[i] = tmp.get(i, 0);  }
    
    return req;
}




vector<double> compute_TT_scaled(MATRIX& T, double scl){
/**

  \param[in] T         - the transformation coordinates (normal modes)
  \param[in] scl       - scaling coefficient 

  This function is auxiliary to compute: 
  - the coefficients of linear electronic coupling in normal modes - Eq. 53
  - the non-equilibrium shifts along the normal modes - Eq. 55
  
  Ref: 
  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109

*/

    /// Check the dimensions of the input and output
    int ndof = T.n_cols;

    vector<double> res(ndof, 0.0);
    for(int i=0; i<ndof; i++){ res[i] = T.get(0, i) * scl; }

    return res;
}


double LVC2GOA_dE(double E0, double E1, vector<double>& omega_nm, vector<double>& d1, vector<double>& d2){
/**
  Computes the dE in the GOA model from the LVC parameters:
  
  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, Eq. (64)

  \param[in] E0        - donor energy level (electronic)
  \param[in] E1        - acceptor energy level (electronic)
  \param[in] omega_nm  - normal modes frequencies
  \param[in] d1        - displacements of equilibrium geometry of donor
  \param[in] d2        - displacements of equilibrium geometry of acceptor

*/

    double dE = E0 - E1;
   
    /// Check the dimensions of the input
    int ndof = omega_nm.size();
    error_msg("d1", d1.size(), "omega_nm", ndof);
    error_msg("d2", d2.size(), "omega_nm", ndof);

    for(int i=0;i<ndof;i++){
      dE -= 0.5*(d1[i]*d1[i] - d2[i]*d2[i])/(omega_nm[i]*omega_nm[i]);
    }

    return dE;

}

vector<double> LVC2GOA_req(vector<double>& omega_nm, vector<double>& d1, vector<double>& d2){
/**
  Computes the dE in the GOA model from the LVC parameters:
  
  Ref:  (1) Sun, X.; Geva, E. JCP, 2016, 145, 064109, Eq. (65)

  \param[in] E0        - donor energy level (electronic)
  \param[in] E1        - acceptor energy level (electronic)
  \param[in] omega_nm  - normal modes frequencies
  \param[in] d1        - displacements of equilibrium geometry of donor
  \param[in] d2        - displacements of equilibrium geometry of acceptor

*/

   
    /// Check the dimensions of the input
    int ndof = omega_nm.size();
    error_msg("d1", d1.size(), "omega_nm", ndof);
    error_msg("d2", d2.size(), "omega_nm", ndof);

    vector<double> req(ndof, 0.0);
    for(int i=0;i<ndof;i++){  req[i] = (d1[i] - d2[i])/(omega_nm[i]*omega_nm[i]);  }

    return req;

}






}/// namespace libfgr
}/// namespace liblibra


