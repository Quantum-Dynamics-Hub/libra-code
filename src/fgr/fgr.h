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

#ifndef FGR_H
#define FGR_H

#include "../math_linalg/liblinalg.h"
#include "../math_random/librandom.h"
#include "../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace librandom;
using namespace libmeigen;


namespace libfgr{


///=============== In fgr_normal_modes.cpp ============================
double S_omega_ohmic(double omega, double eta); //ohmic with decay spectral density
double S_omega_drude(double omega, double eta);//another spectral density
double S_omega_gaussian(double omega, double eta, double sigma, double omega_op);//gaussian spectral density
/*
double J_omega_ohmic(double omega, double eta);//bath Ohmic SD
double J_omega_ohmic_eff(double omega, double eta); //effective SD for Ohmic bath
void Integrand_LSC(double omega, double t, double &re, double &im);
void Integrand_LSC_inh(double omega, double t, double &re, double &im);
void Integrand_CL_avg(double omega, double t, double &re, double &im);
void Integrand_CL_donor(double omega, double t, double &re, double &im);
void Integrand_2cumu(double omega, double t, double &re, double &im);
void Integrand_2cumu_inh(double omega, double t, double &re, double &im);

*/

double eq_shift(double Er, double Omega);
double reorganization_energy(double y0, double Omega);
double reorganization_energy(vector<double>& omega_nm, vector<double>& req_nm);
double diabat_crossing(double dE, double Er, double y0);
double coupling_Condon(double gamma, double dE, double Er, double y0);
double coupling_non_Condon(double y, double gamma, double dE, double Er, double y0);

vector<double> normal_modes(vector<double>& omega, vector<double>& coeff, MATRIX& T);
vector<double> compute_req(vector<double>& omega, vector<double>& coeff, double y0, MATRIX& T);
vector<double> compute_TT_scaled(MATRIX& T, double scl);

double LVC2GOA_dE(double E0, double E1, vector<double>& omega_nm, vector<double>& d1, vector<double>& d2);
vector<double> LVC2GOA_req(vector<double>& omega_nm, vector<double>& d1, vector<double>& d2);


///=============== In fgr.cpp ============================

complex<double> Integrand_NE_exact(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta);
complex<double> Linear_NE_exact(double tp, double tau, double gamma, double omega, double req, double shift, double beta);
complex<double> ACF_NE_exact(double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
                             vector<double>& req_nm, vector<double>& shift_NE, double beta, int type);

complex<double> Integrand_NE_LSC(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta);
complex<double> Linear_NE_LSC(double tp, double tau, double gamma, double omega, double req, double shift, double beta);
complex<double> ACF_NE_LSC(double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
                           vector<double>& req_nm, vector<double>& shift_NE, double beta, int type);

complex<double> Integrand_NE_CAV(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta);
complex<double> Linear_NE_CAV(double tp, double tau, double gamma, double omega, double req, double shift, double beta);
complex<double> ACF_NE_CAV(double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
                           vector<double>& req_nm, vector<double>& shift_NE, double beta, int type);

complex<double> Integrand_NE_CD(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta);
complex<double> Linear_NE_CD(double tp, double tau, double gamma, double omega, double req, double shift, double beta);
complex<double> ACF_NE_CD(double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
                          vector<double>& req_nm, vector<double>& shift_NE, double beta, int type);

complex<double> Integrand_NE_W0(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta);
complex<double> Linear_NE_W0(double tp, double tau, double gamma, double omega, double req, double shift, double beta);
complex<double> ACF_NE_W0(double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
                          vector<double>& req_nm, vector<double>& shift_NE, double beta, int type);

complex<double> Integrand_NE_C0(double tp, double tau, double omega_DA, double omega, double req, double shift, double beta);
complex<double> Linear_NE_C0(double tp, double tau, double gamma, double omega, double req, double shift, double beta);
complex<double> ACF_NE_C0(double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
                              vector<double>& req_nm, vector<double>& shift_NE, double beta, int type);


double NEFGRL_rate(double t, double omega_DA, double V,
                   vector<double>& omega_nm, vector<double>& gamma_nm,
                   vector<double>& req_nm, vector<double>& shift_NE,
                   int method, double beta, int type, double dtau);

MATRIX NEFGRL_population(double omega_DA, double V,
                         vector<double>& omega_nm, vector<double>& gamma_nm,
                         vector<double>& req_nm, vector<double>& shift_NE,
                         int method, double beta, int type, double dtau, double tmax, double dt);


}// namespace libfgr
}// liblibra

#endif // FGR_H
