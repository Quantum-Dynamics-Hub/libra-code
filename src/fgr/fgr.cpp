/*********************************************************************************
* Copyright (C) 2019 Xiang Sun, Alexey V. Akimov
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
  \file fgr.cpp
  \brief This file implements the Fermi Golden Rule (FGR)-based formula for quantum dynamics rates. 
*/

#include "fgr.h"


/// liblibra namespace
namespace liblibra{

/// libivr namespace
namespace libfgr{




complex<double> Integrand_NE_exact(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double wtp = omega*tp;
  double wb = omega*beta;
  double pref = omega*req*req*0.5;

  re = -pref*(1.0-cos(wt))/tanh(0.5*wb);
  im = -pref*sin(wt) - omega*req*shift*(sin(wtp) - sin(wtp - wt));

  return complex<double>(re, im);
}


complex<double> Linear_NE_exact(double tp, double tau, double gamma, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wtp = omega*tp;
  double wt  = omega*tau;
  double wb  = omega*beta;

  double cs = cos(wt);
  double si = sin(wt);
  double Coth = 1.0 / tanh(wb*0.5);

  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);

/*
  Xiang's stuff
  double u1_re = 0.5/omega * Coth * cos(wt);
  double u1_im =-0.5/omega * sin(wt);
  double u2_re = req * (1.0-cos(wt));
  double u2_im = req * Coth * sin(wt);

  re = u1_re + 0.25 * (u2_re - 2.0*shift*cos(wtp)) * (u2_re - 2.0*shift*cos(wtp - wt)) - 0.25 * u2_im * u2_im;
  im = u1_im + 0.25 * (u2_re - 2.0*shift*cos(wtp)) * u2_im + 0.25 * (u2_re - 2.0*shift*cos(wtp - wt)) * u2_im;
*/
  // Explicit way
  complex<double> res;

  res = ( one*(req - 2.0 * shift * cos(wtp) - req * cs) + eye*(req * Coth * si)  );
  res *= 0.25*( one * (req - 2.0* shift * cos(wtp - wt)  - req * cs ) + eye * (Coth * si) );
  res += (0.5/omega) * ( Coth * cs * one - si * eye);

  res *= gamma * gamma;

  

  return res; //gamma*gamma*complex<double>(re, im);

}


complex<double> ACF_NE_exact(double tp, double tau, double omega_DA, double V,
                             vector<double>& omega_nm, vector<double>& gamma_nm,
                             vector<double>& req_nm, vector<double>& shift_NE, 
                             double beta, int type){
    /**
    Computes the ACF, according to:
    Sun, X.; Geva, E. J. Chem. Phys. 145, 064109, 2016 -  Eq. 56

    tp   [a.u. of time]        - time at which the rate constant is to be determined
    tau  [a.u. of time]        - lag time - the main argument of the ACF
    type                       - flag that detemines if the ACF if for Condon (type=0), or non-Condon (type=1) case
    omega_DA [Ha]              - the energy gap of between donor and acceptor levels (hbar = 1)
    omega_nm [Ha]              - energies of the normal modes (hbar = 1)
    gamma_nm [Ha]              - the electron-phonon coupling of each normal mode
    req_nm   [Bohr]            - shift of the equilibrium position of each noraml mode coordinate upon the charge transfer
    shift_nm [Bohr]            - non-equilibrium displacement of each normal mode
    beta [Ha^-1]  1/(kB*T)     - the inverse of thermal energy

    */

    int n_omega = omega_nm.size();
    complex<double> argg(0.0, 0.0);
    complex<double> ampl(0.0, 0.0);
    argg = complex<double>(0.0, omega_DA * tau);

    for(int w=0; w<n_omega; w++){
      argg += Integrand_NE_exact(tp, tau, omega_DA, omega_nm[w], req_nm[w], shift_NE[w], beta); 

      if(type==1){
          ampl += Linear_NE_exact(tp, tau, gamma_nm[w], omega_nm[w], req_nm[w], shift_NE[w], beta);
      }
    }

    if(type==0){  ampl = complex<double>(V*V, 0.0); }

    // C = exp(argg) * ampl = exp[Re(argg)+i*Im(argg)] * [Re(ampl) + i*Im(ampl)] = 
    // = exp(Re(argg)) * [cos(Im(argg)) + i*sin(Im(argg))] * [Re(ampl) + i*Im(ampl)]
    double C_re = exp(argg.real()) * ( cos(argg.imag()) * ampl.real() - sin(argg.imag()) *  ampl.imag() ); // this is Re(C)
    double C_im = exp(argg.real()) * ( cos(argg.imag()) * ampl.imag() + sin(argg.imag()) *  ampl.real() ); // this is Im(C)

    return complex<double>(C_re, C_im);

}




complex<double> Integrand_NE_LSC(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta){
/**
*/

  return Integrand_NE_exact(tp, tau, omega_DA, omega, req, shift, beta);

}

complex<double> Linear_NE_LSC(double tp, double tau, double gamma, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wtp = omega*tp;
  double wt  = omega*tau;
  double wb  = omega*beta;

  double Coth = 1.0 / tanh(wb*0.5);
  re  = shift*shift*cos(wtp)*cos(wtp-wt);
  re += (Coth/omega)*0.5*cos(wt);
  re -= req*req*0.5*Coth*Coth * pow(sin(0.5*wt),2) * (cos(4.0*wtp - 2.0*wt) + cos(wt));

  im =  (1.0-2.0*cos(wt))*sin(wtp) + sin(wtp-2.0*wt) - 4.0*cos(3.0*wtp - 1.5*wt)*sin(0.5*wt);
  im *= 0.25*req*shift*Coth;

  return gamma*gamma*complex<double>(re, im);

}

complex<double> ACF_NE_LSC(double tp, double tau, double omega_DA, double V,
                           vector<double>& omega_nm, vector<double>& gamma_nm,
                           vector<double>& req_nm, vector<double>& shift_NE, 
                           double beta, int type){
    /**
    Computes the ACF, according to:
    Sun, X.; Geva, E. J. Chem. Phys. 145, 064109, 2016 -  Eq. 57

    tp   [a.u. of time]        - time at which the rate constant is to be determined
    tau  [a.u. of time]        - lag time - the main argument of the ACF
    type                       - flag that detemines if the ACF if for Condon (type=0), or non-Condon (type=1) case
    omega_DA [Ha]              - the energy gap of between donor and acceptor levels (hbar = 1)
    omega_nm [Ha]              - energies of the normal modes (hbar = 1)
    gamma_nm [Ha]              - the electron-phonon coupling of each normal mode
    req_nm   [Bohr]            - shift of the equilibrium position of each noraml mode coordinate upon the charge transfer
    shift_nm [Bohr]            - non-equilibrium displacement of each normal mode
    beta [Ha^-1]  1/(kB*T)     - the inverse of thermal energy

    */

    int n_omega = omega_nm.size();
    complex<double> argg(0.0, 0.0);
    complex<double> ampl(0.0, 0.0);
    argg = complex<double>(0.0, omega_DA*tau);

    for(int w=0; w<n_omega; w++){
      argg += Integrand_NE_LSC(tp, tau, omega_DA, omega_nm[w], req_nm[w], shift_NE[w], beta); 

      if(type==1){
          ampl += Linear_NE_LSC(tp, tau, gamma_nm[w], omega_nm[w], req_nm[w], shift_NE[w], beta);
      }
    }

    if(type==0){  ampl = complex<double>(V*V, 0.0); }

    // C = exp(argg) * ampl = exp[Re(argg)+i*Im(argg)] * [Re(ampl) + i*Im(ampl)] = 
    // = exp(Re(argg)) * [cos(Im(argg)) + i*sin(Im(argg))] * [Re(ampl) + i*Im(ampl)]
    double C_re = exp(argg.real()) * ( cos(argg.imag()) * ampl.real() - sin(argg.imag()) *  ampl.imag() ); // this is Re(C)
    double C_im = exp(argg.real()) * ( cos(argg.imag()) * ampl.imag() + sin(argg.imag()) *  ampl.real() ); // this is Im(C)

    return complex<double>(C_re, C_im);

}




complex<double> Integrand_NE_CAV(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double wtp = omega*tp;

  re = -(req*req/beta)*(1.0-cos(wt));
  im = - omega*req*req*0.5*sin(wt) - omega*req*shift * (sin(wtp) + sin(wt - wtp));

  return complex<double>(re, im);
}


complex<double> Linear_NE_CAV(double tp, double tau, double gamma, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wtp = omega*tp;
  double wt  = omega*tau;
  double wb  = omega*beta;

  re = shift*shift*cos(wtp)*cos(wtp-wt) + 1.0/(wb*omega) * cos(wt) 
       - req*req* 0.5 / pow(wb*0.5,2) * pow(sin(0.5*wt),2) * (cos(4.0*wtp - 2.0*wt) + cos(wt));
  im = 0.5*req*shift/wb * ( (1.0-2.0*cos(wt))*sin(wtp) + sin(wtp-2.0*wt) - 4.0*cos(3.0*wtp - 1.5*wt)*sin(0.5*wt) );

  return gamma*gamma*complex<double>(re, im);

}

complex<double> ACF_NE_CAV(double tp, double tau, double omega_DA, double V,
                           vector<double>& omega_nm, vector<double>& gamma_nm,
                           vector<double>& req_nm, vector<double>& shift_NE,
                           double beta, int type){
    /**
    Computes the ACF, according to:
    Sun, X.; Geva, E. J. Chem. Phys. 145, 064109, 2016 -  Eq. 58

    tp   [a.u. of time]        - time at which the rate constant is to be determined
    tau  [a.u. of time]        - lag time - the main argument of the ACF
    type                       - flag that detemines if the ACF if for Condon (type=0), or non-Condon (type=1) case
    omega_DA [Ha]              - the energy gap of between donor and acceptor levels (hbar = 1)
    omega_nm [Ha]              - energies of the normal modes (hbar = 1)
    gamma_nm [Ha]              - the electron-phonon coupling of each normal mode
    req_nm   [Bohr]            - shift of the equilibrium position of each noraml mode coordinate upon the charge transfer
    shift_nm [Bohr]            - non-equilibrium displacement of each normal mode
    beta [Ha^-1]  1/(kB*T)     - the inverse of thermal energy

    */

    int n_omega = omega_nm.size();
    complex<double> argg(0.0, 0.0);
    complex<double> ampl(0.0, 0.0);
    argg = complex<double>(0.0, omega_DA*tau);

    for(int w=0; w<n_omega; w++){
      argg += Integrand_NE_CAV(tp, tau, omega_DA, omega_nm[w], req_nm[w], shift_NE[w], beta); 

      if(type==1){
          ampl += Linear_NE_CAV(tp, tau, gamma_nm[w], omega_nm[w], req_nm[w], shift_NE[w], beta);
      }
    }

    if(type==0){  ampl = complex<double>(V*V, 0.0); }

    // C = exp(argg) * ampl = exp[Re(argg)+i*Im(argg)] * [Re(ampl) + i*Im(ampl)] = 
    // = exp(Re(argg)) * [cos(Im(argg)) + i*sin(Im(argg))] * [Re(ampl) + i*Im(ampl)]
    double C_re = exp(argg.real()) * ( cos(argg.imag()) * ampl.real() - sin(argg.imag()) *  ampl.imag() ); // this is Re(C)
    double C_im = exp(argg.real()) * ( cos(argg.imag()) * ampl.imag() + sin(argg.imag()) *  ampl.real() ); // this is Im(C)

    return complex<double>(C_re, C_im);

}






complex<double> Integrand_NE_CD(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double wtp = omega*tp;
  double prefactor = 0.5*omega*req*req*wt;

  re = -(req*req/beta)*(1.0-cos(wt));
  im = - prefactor - omega*req*shift * (sin(wtp) + sin(wt - wtp));

  return complex<double>(re, im);
}

complex<double> Linear_NE_CD(double tp, double tau, double gamma, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wtp = omega*tp;
  double wt  = omega*tau;
  double wb  = omega*beta;

  re = shift*shift*cos(wtp)*cos(wtp-wt) + cos(wt)/(wb*omega) -  pow((req*sin(wt)/wb),2);
  im = -4.0*shift*req/wb * cos(wtp-0.5*wt) * pow(cos(0.5*wt),2) * sin(0.5*wt);

  return gamma*gamma*complex<double>(re, im);

}

complex<double> ACF_NE_CD(double tp, double tau, double omega_DA, double V,
                          vector<double>& omega_nm, vector<double>& gamma_nm,
                          vector<double>& req_nm, vector<double>& shift_NE, 
                          double beta, int type){
    /**
    Computes the ACF, according to:
    Sun, X.; Geva, E. J. Chem. Phys. 145, 064109, 2016 -  Eq. 59

    tp   [a.u. of time]        - time at which the rate constant is to be determined
    tau  [a.u. of time]        - lag time - the main argument of the ACF
    type                       - flag that detemines if the ACF if for Condon (type=0), or non-Condon (type=1) case
    omega_DA [Ha]              - the energy gap of between donor and acceptor levels (hbar = 1)
    omega_nm [Ha]              - energies of the normal modes (hbar = 1)
    gamma_nm [Ha]              - the electron-phonon coupling of each normal mode
    req_nm   [Bohr]            - shift of the equilibrium position of each noraml mode coordinate upon the charge transfer
    shift_nm [Bohr]            - non-equilibrium displacement of each normal mode
    beta [Ha^-1]  1/(kB*T)     - the inverse of thermal energy

    */

    int n_omega = omega_nm.size();
    complex<double> argg(0.0, 0.0);
    complex<double> ampl(0.0, 0.0);
    argg = complex<double>(0.0, omega_DA*tau);
            
    for(int w=0; w<n_omega; w++){
      argg += Integrand_NE_CD(tp, tau, omega_DA, omega_nm[w], req_nm[w], shift_NE[w], beta); 

      if(type==1){
          ampl += Linear_NE_CD(tp, tau, gamma_nm[w], omega_nm[w], req_nm[w], shift_NE[w], beta);
      }
    }

    if(type==0){  ampl = complex<double>(V*V, 0.0); }

    // C = exp(argg) * ampl = exp[Re(argg)+i*Im(argg)] * [Re(ampl) + i*Im(ampl)] = 
    // = exp(Re(argg)) * [cos(Im(argg)) + i*sin(Im(argg))] * [Re(ampl) + i*Im(ampl)]
    double C_re = exp(argg.real()) * ( cos(argg.imag()) * ampl.real() - sin(argg.imag()) *  ampl.imag() ); // this is Re(C)
    double C_im = exp(argg.real()) * ( cos(argg.imag()) * ampl.imag() + sin(argg.imag()) *  ampl.real() ); // this is Im(C)

    return complex<double>(C_re, C_im);

}




complex<double> Integrand_NE_W0(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double prefactor = 0.5*omega*req*req*wt;

  re = -( prefactor / tanh(0.5*beta*omega) ) * wt*wt*0.5;
  im = - prefactor - omega*req*shift * cos(omega*tp) * wt;

  return complex<double>(re, im);
}

complex<double> Linear_NE_W0(double tp, double tau, double gamma, double omega, double req, double shift, double beta){
/**
*/
  double re, im;
  double wtp = omega*tp;
  double argg = req*tau*omega;

  double Coth = 1.0 / tanh(beta*omega*0.5);
  re = Coth*0.5/omega + 0.5*shift*shift*(1.0+cos(2.0*wtp)) - 0.25*argg*argg*Coth*Coth;
  im = - Coth * argg * shift * cos(wtp);

  return gamma*gamma*complex<double>(re, im);

}

complex<double> ACF_NE_W0(double tp, double tau, double omega_DA, double V,
                          vector<double>& omega_nm, vector<double>& gamma_nm,
                          vector<double>& req_nm, vector<double>& shift_NE, 
                          double beta, int type){
    /**
    Computes the ACF, according to:
    Sun, X.; Geva, E. J. Chem. Phys. 145, 064109, 2016 -  Eq. 60

    tp   [a.u. of time]        - time at which the rate constant is to be determined
    tau  [a.u. of time]        - lag time - the main argument of the ACF
    type                       - flag that detemines if the ACF if for Condon (type=0), or non-Condon (type=1) case
    omega_DA [Ha]              - the energy gap of between donor and acceptor levels (hbar = 1)
    omega_nm [Ha]              - energies of the normal modes (hbar = 1)
    gamma_nm [Ha]              - the electron-phonon coupling of each normal mode
    req_nm   [Bohr]            - shift of the equilibrium position of each noraml mode coordinate upon the charge transfer
    shift_nm [Bohr]            - non-equilibrium displacement of each normal mode
    beta [Ha^-1]  1/(kB*T)     - the inverse of thermal energy

    */

    int n_omega = omega_nm.size();
    complex<double> argg(0.0, 0.0);
    complex<double> ampl(0.0, 0.0);
    argg = complex<double>(0.0, omega_DA*tau);
            
    for(int w=0; w<n_omega; w++){
      argg += Integrand_NE_W0(tp, tau, omega_DA, omega_nm[w], req_nm[w], shift_NE[w], beta); 

      if(type==1){
          ampl += Linear_NE_W0(tp, tau, gamma_nm[w], omega_nm[w], req_nm[w], shift_NE[w], beta);
      }
    }

    if(type==0){  ampl = complex<double>(V*V, 0.0); }

    // C = exp(argg) * ampl = exp[Re(argg)+i*Im(argg)] * [Re(ampl) + i*Im(ampl)] = 
    // = exp(Re(argg)) * [cos(Im(argg)) + i*sin(Im(argg))] * [Re(ampl) + i*Im(ampl)]
    double C_re = exp(argg.real()) * ( cos(argg.imag()) * ampl.real() - sin(argg.imag()) *  ampl.imag() ); // this is Re(C)
    double C_im = exp(argg.real()) * ( cos(argg.imag()) * ampl.imag() + sin(argg.imag()) *  ampl.real() ); // this is Im(C)

    return complex<double>(C_re, C_im);

}




complex<double> Integrand_NE_C0(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta){
/**
  omega -  frequency of the normal mode [a.u.] 
  tp    -  [a.u.]
  tau   -  []
  shift -  []
  req   -  []  
*/
  double re, im;
  double wt = omega*tau;
  double prefactor = 0.5*omega*req*req*wt;

  re = -prefactor * tau/beta;
  im = -prefactor - omega*req*shift * cos(omega*tp)*wt;

  return complex<double>(re, im);
}

                                 
complex<double> Linear_NE_C0(double tp, double tau, double gamma, double omega, double req, double shift, double beta){
/**
  omega -  frequency of the normal mode [a.u.] 
  tp    -  [a.u.]
  tau   -  []
  shift -  []
  req   -  []  
*/
  double re, im;
  double wtp = omega*tp;
  double argg = req*tau/beta;

  re = 1.0/(beta*omega*omega) + 0.5*shift*shift*(1.0 + cos(2.0*wtp)) - argg*argg;
  im = -2.0*argg*shift * cos(wtp);

  return gamma*gamma*complex<double>(re, im);
}

complex<double> ACF_NE_C0(double tp, double tau, double omega_DA, double V,
                              vector<double>& omega_nm, vector<double>& gamma_nm,
                              vector<double>& req_nm, vector<double>& shift_NE, 
                              double beta, int type){
    /**
    Computes the ACF, according to:
    Sun, X.; Geva, E. J. Chem. Phys. 145, 064109, 2016 -  Eq. 61

    tp   [a.u. of time]        - time at which the rate constant is to be determined
    tau  [a.u. of time]        - lag time - the main argument of the ACF
    type                       - flag that detemines if the ACF if for Condon (type=0), or non-Condon (type=1) case
    omega_DA [Ha]              - the energy gap of between donor and acceptor levels (hbar = 1)
    omega_nm [Ha]              - energies of the normal modes (hbar = 1)
    gamma_nm [Ha]              - the electron-phonon coupling of each normal mode
    req_nm   [Bohr]            - shift of the equilibrium position of each noraml mode coordinate upon the charge transfer
    shift_nm [Bohr]            - non-equilibrium displacement of each normal mode
    beta [Ha^-1]  1/(kB*T)     - the inverse of thermal energy

    */

    int n_omega = omega_nm.size();
    complex<double> argg(0.0, 0.0);
    complex<double> ampl(0.0, 0.0);
    argg = complex<double>(0.0, omega_DA*tau);
            
    for(int w=0; w<n_omega; w++){
      argg += Integrand_NE_C0(tp, tau, omega_DA, omega_nm[w], req_nm[w], shift_NE[w], beta); 

      if(type==1){
          ampl += Linear_NE_C0(tp, tau, gamma_nm[w], omega_nm[w], req_nm[w], shift_NE[w], beta);
      }
    }

    if(type==0){  ampl = complex<double>(V*V, 0.0); }

    // C = exp(argg) * ampl = exp[Re(argg)+i*Im(argg)] * [Re(ampl) + i*Im(ampl)] = 
    // = exp(Re(argg)) * [cos(Im(argg)) + i*sin(Im(argg))] * [Re(ampl) + i*Im(ampl)]
    double C_re = exp(argg.real()) * ( cos(argg.imag()) * ampl.real() - sin(argg.imag()) *  ampl.imag() ); // this is Re(C)
    double C_im = exp(argg.real()) * ( cos(argg.imag()) * ampl.imag() + sin(argg.imag()) *  ampl.real() ); // this is Im(C)

    return complex<double>(C_re, C_im);

}





double NEFGRL_rate(double tp, double omega_DA, double V, 
                   vector<double>& omega_nm, vector<double>& gamma_nm,
                   vector<double>& req_nm, vector<double>& shift_NE,
                   int method, double beta, int type, double dtau
                  ){
/** 
  Computes the nonequilibrium FGR rates in Condon (type=0) and non-Condon(type=1) cases using normal modes

  Sun, X.; Geva, E. J. Chem. Phys. 145, 064109, 2016 -  Eq. 6


  k(t') = 2/(hbar^2) * Re { int_0^t' { dtau * C(tau) }  }

  \param[in]     V      - D-A coupling [a.u. of energy] 
  \param[in]  omega_DA  - D-A electronic transition frequency [a.u.]
  \param[in]      t     - time since the initial photoexcitation [a.u. of time]
  \param[in]   dtau     - integration timestep for backward propagation, used to take the integral above
                         (tau = dtau*n < t, where n is integer)
  \param[in] omega_nm   - frequencies of the bath normal modes
  \param[in] gamma_nm   - couplings of the bath normal modes to the primary mode
  \param[in] req_nm     - equilibrum position displacement in the excited state vs. the ground state
  \param[in] shift_NE   - non-equilibrium shift of the normal modes
  \param[in] method     - flag that specifies which method to use:
                        0 - Exact quantum result
                        1 - LSC
                        2 - CAV
                        3 - CD
                        4 - W0
                        5 - C0
  \param[in] beta       - inverse thermal energy [unitless]

  Returns: the instantaneous rate at a given time

*/

  double k = 0.0;  // rate constant

  int N = int(tp/dtau) + 1; //static_cast<int>(tp/dtau); 

  complex<double> C(0.0, 0.0);

  for(int n=0; n<N; n++){ //tau index

    double tau = n * dtau;

    if(tau<=tp){

      switch(method){
        case 0: {  C = ACF_NE_exact(tp, tau, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, beta, type);    } break;
        case 1: {  C = ACF_NE_LSC(tp, tau, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, beta, type);    } break;
        case 2: {  C = ACF_NE_CAV(tp, tau, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, beta, type);    } break;
        case 3: {  C = ACF_NE_CD(tp, tau, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, beta, type);    } break;
        case 4: {  C = ACF_NE_W0(tp, tau, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, beta, type);    } break;
        case 5: {  C = ACF_NE_C0(tp, tau, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, beta, type);    } break;
        default: {
          cout<<"Method "<<method<<" is not available. Please choose from the following options:\n";
          cout<<"0 - Exact\n";
          cout<<"1 - LSC\n";
          cout<<"2 - CAV\n";
          cout<<"3 - CD\n";
          cout<<"4 - W0\n";
          cout<<"5 - C0\n";
          cout<<"Exiting now...\n";
          exit(0);
        }// default
      }// switch
    }// tau<=tp

    k += C.real(); // we only need Re(k)  
  }

  k *= (2.0*dtau);  // yields k = k(t')
           
  return k;

}



MATRIX NEFGRL_population(double omega_DA, double V,
                         vector<double>& omega_nm, vector<double>& gamma_nm,
                         vector<double>& req_nm, vector<double>& shift_NE,
                         int method, double beta, int type, double dtau, double tmax, double dt 
                        ){
/**
  Noneq FGR in non-Condon case (linear coupling) using normal modes

  k(t') = 2/(hbar^2) * Re { int_0^t' { dtau * C(tau) }  }

  \param[in]     V      - D-A coupling [a.u. of energy] 
  \param[in]  omega_DA  - D-A electronic transition frequency [a.u.]
  \param[in]   dtau     - integration timestep for backward propagation, used to take the integral above
                         (tau = dtau*n < t, where n is integer)
  \param[in] omega_nm   - frequencies of the bath normal modes
  \param[in] gamma_nm   - couplings of the bath normal modes to the primary mode
  \param[in] req_nm     - equilibrum position displacement in the excited state vs. the ground state
  \param[in] shift_NE   - non-equilibrium shift of the normal modes
  \param[in] method     - flag that specifies which method to use:
                        0 - Exact quantum result
                        1 - LSC
                        2 - CAV
                        3 - CD
                        4 - W0
                        5 - C0
  \param[in]   tmax     - time since the initial photoexcitation [a.u. of time]
  \param[in]     dt     - integration timestep for the trajectory [a.u. of time]
  \param[in]   beta     - hbar*omega_cutoff / kB*T inverse "temperature" [unitless]

  Returns: the matrix with: current time, instantaneous rate, population on the donor state

*/

  int nsteps = (int)(tmax/dt)+1;
  MATRIX res(nsteps, 3);

  double sum = 0.0; //probability of donor state

  for(int step=0; step<nsteps; step++) {  // t' 

    // k = k(t')
    double k = NEFGRL_rate(step*dt, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, method, beta, type, dtau);
        
    sum += k * dt; 
    double P = exp(-sum);  //exp(- int dt' k(t'))

    res.set(step, 0, step*dt);
    res.set(step, 1, k);
    res.set(step, 2, P);

  }

  return res;

}


}/// namespace libfgr
}/// namespace liblibra


