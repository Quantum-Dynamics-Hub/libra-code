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
  \file fgr.cpp
  \brief This file implements the Fermi Golden Rule (FGR)-based formula for quantum
  dynamics rates. 
*/

#include "fgr.h"


/// liblibra namespace
namespace liblibra{

/// libivr namespace
namespace libfgr{




complex<double> Integrand_NE_exact(double omega_DA, double omega, double tp, double tau, double shift, double req, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double wtp = omega*tp;
  double wb = omega*beta;
  double pref = omega*req*req*0.5;

  re = -pref*(1.0-cos(wt))/tanh(0.5*wb);
  im = omega_DA*tau - pref*sin(wt) - omega*req*shift*(sin(wtp) + sin(wt - wtp));

  return complex<double>(re, im);
}



complex<double> Linear_NE_exact(double gamma, double omega, double tp, double tau, double shift, double req, double beta){
/**
*/
  double re, im;
  double wtp = omega*tp;
  double wt  = omega*tau;
  double wb  = omega*beta;

  double Coth = 1.0 / tanh(wb*0.5);
  double u1_re = 0.5/omega * Coth * cos(wt);
  double u1_im =-0.5/omega * sin(wt);
  double u2_re = req * (1.0-cos(wt));
  double u2_im = req * Coth * sin(wt);

  re = u1_re + 0.25 * (u2_re - 2.0*shift*cos(wtp)) * (u2_re - 2.0*shift*cos(wtp - wt)) - 0.25 * u2_im * u2_im;
  im = u1_im + 0.25 * (u2_re - 2.0*shift*cos(wtp)) * u2_im + 0.25 * (u2_re - 2.0*shift*cos(wtp - wt)) * u2_im;

  return gamma*gamma*complex<double>(re, im);

}





complex<double> Integrand_NE_LSC(double omega_DA, double omega, double tp, double tau, double shift, double req, double beta){
/**
*/

  return Integrand_NE_exact(omega_DA, omega, tp, tau, shift, req, beta);

}

complex<double> Linear_NE_LSC(double gamma, double omega, double tp, double tau, double shift, double req, double beta){
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




complex<double> Integrand_NE_CAV(double omega_DA, double omega, double tp, double tau, double shift, double req, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double wtp = omega*tp;

  re = -(req*req/beta)*(1.0-cos(wt));
  im = omega_DA*tau - omega*req*req*0.5*sin(wt) - omega*req*shift * (sin(wtp) + sin(wt - wtp));

  return complex<double>(re, im);
}


complex<double> Linear_NE_CAV(double gamma, double omega, double tp, double tau, double shift, double req, double beta){
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




complex<double> Integrand_NE_CD(double omega_DA, double omega, double tp, double tau, double shift, double req, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double wtp = omega*tp;
  double prefactor = 0.5*omega*req*req*wt;

  re = -(req*req/beta)*(1.0-cos(wt));
  im = omega_DA*tau - prefactor - omega*req*shift * (sin(wtp) + sin(wt - wtp));

  return complex<double>(re, im);
}

complex<double> Linear_NE_CD(double gamma, double omega, double tp, double tau, double shift, double req, double beta){
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



complex<double> Integrand_NE_W0(double omega_DA, double omega, double tp, double tau, double shift, double req, double beta){
/**
*/
  double re, im;
  double wt = omega*tau;
  double prefactor = 0.5*omega*req*req*wt;

  re = -( prefactor / tanh(0.5*beta*omega) ) * wt*wt*0.5;
  im = omega_DA*tau - prefactor - omega*req*shift * cos(omega*tp) * wt;

  return complex<double>(re, im);
}

complex<double> Linear_NE_W0(double gamma, double omega, double tp, double tau, double shift, double req, double beta){
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



complex<double> Integrand_NE_Marcus(double omega_DA, double omega, double tp, double tau, double shift, double req, double beta){
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
  im = omega_DA*tau  -prefactor - omega*req*shift * cos(omega*tp)*wt;

  return complex<double>(re, im);
}

                                 
complex<double> Linear_NE_Marcus(double gamma, double omega, double tp, double tau, double shift, double req, double beta){
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



double NEFGRL_rate(double V, double omega_DA, double t, double dtau,  
                   vector<double>& omega_nm, vector<double>& gamma_nm,
                   vector<double>& req_nm, vector<double>& shift_NE,
                   int method, double beta
                  ){
/** 
  Noneq FGR in non-Condon case (linear coupling) using normal modes

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
                        5 - Marcus
  \param[in] beta       - inverse thermal energy [unitless]

  Returns: the instantaneous rate at a given time

*/

  double k = 0.0;  // rate constant

  int N = static_cast<int>(t/dtau); 
  int n_omega = omega_nm.size();

  for(int n=0; n<N; n++){ //tau index

    double tau = n * dtau;

    complex<double> argg(0.0, 0.0);
    complex<double> ampl(0.0, 0.0);

    for(int w=0; w<n_omega; w++){

      switch(method){

        case 0: {
          argg += Integrand_NE_exact(omega_DA, omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta); 
          ampl += Linear_NE_exact(gamma_nm[w], omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta);
        } break;

        case 1: {
          argg += Integrand_NE_LSC(omega_DA, omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta);
          ampl += Linear_NE_LSC(gamma_nm[w], omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta);
        } break;

        case 2: {
          argg += Integrand_NE_CAV(omega_DA, omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta); 
          ampl += Linear_NE_CAV(gamma_nm[w], omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta);
        } break;

        case 3: {
          argg += Integrand_NE_CD(omega_DA, omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta); 
          ampl += Linear_NE_CD(gamma_nm[w], omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta);
        } break;

        case 4: {
          argg += Integrand_NE_W0(omega_DA, omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta); 
          ampl += Linear_NE_W0(gamma_nm[w], omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta);
        } break;
    
        case 5: {
          argg += Integrand_NE_Marcus(omega_DA, omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta); 
          ampl += Linear_NE_Marcus(gamma_nm[w], omega_nm[w], t, tau, shift_NE[w], req_nm[w], beta);
        } break;

        default: {
          cout<<"Method "<<method<<" is not available. Please choose from the following options:\n";
          cout<<"0 - Exact\n";
          cout<<"1 - LSC\n";
          cout<<"2 - CAV\n";
          cout<<"3 - CD\n";
          cout<<"4 - W0\n";
          cout<<"5 - Marcus\n";
          cout<<"Exiting now...\n";
          exit(0);
        }

      }// switch

    }
 
    // k = exp(argg) * ampl = exp[Re(argg)+i*Im(argg)] * [Re(ampl) + i*Im(ampl)] = 
    // = exp(Re(argg)) * [cos(Im(argg)) + i*sin(Im(argg))] * [Re(ampl) + i*Im(ampl)]
    // but we only need Re(k)
    k += exp(argg.real()) * ( cos(argg.imag()) * ampl.real() - sin(argg.imag()) *  ampl.imag() ); // this is Re(k)
  
  }

  k *= (2.0*dtau);
           
  return k;

}


MATRIX NEFGRL_population(double V, double omega_DA, double dtau,  
                         vector<double>& omega_nm, vector<double>& gamma_nm,
                         vector<double>& req_nm, vector<double>& shift_NE,
                         int method, double T, double dt, double beta
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
                        5 - Marcus
  \param[in]      T     - time since the initial photoexcitation [a.u. of time]
  \param[in]     dt     - integration timestep for the trajectory [a.u. of time]
  \param[in]   beta     - hbar*omega_cutoff / kB*T inverse "temperature" [unitless]

  Returns: the matrix with: current time, instantaneous rate, population on the donor state

*/

  int nsteps = (int)(T/dt)+1;
  MATRIX res(nsteps, 3);

  double sum = 0.0; //probability of donor state

  for(int step=0; step<nsteps; step++) {  // t' 

    double k = NEFGRL_rate(V, omega_DA, step*dt, dtau, omega_nm, gamma_nm, req_nm, shift_NE, method, beta);
        
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


