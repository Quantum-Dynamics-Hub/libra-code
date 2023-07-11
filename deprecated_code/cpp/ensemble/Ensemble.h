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
/**
  \file Ensemble.h
  \brief The file implements Python export function
    
*/

#ifndef ENSEMBLE_H
#define ENSEMBLE_H


#include "../../math_linalg/liblinalg.h"
#include "../../math_random/librandom.h"
#include "../../math_specialfunctions/libspecialfunctions.h"
#include "../../Units.h"

#include "../../hamiltonian/libhamiltonian.h"
#include "../nuclear/libnuclear.h"
#include "../electronic/libelectronic.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace librandom;
using namespace libspecialfunctions;

using namespace libhamiltonian;
using namespace libhamiltonian::libhamiltonian_generic;
using namespace libhamiltonian::libhamiltonian_model;


/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;


/// libensemble namespace
namespace libensemble{


class Ensemble{
/** 
  \brief Ensemble class that is designed to propagate many trajectories (entangled or independent) at once

  Object of Ensemble type essentially represents a time-dependent electron-nuclear wavefunction
  el - represents electronic component
  mol - represent nuclear component
  This is a semiclassical representation of a wavefunction via a swarm (ensemble) of trajectories
  
  PSI(t) = sum_k {  delta((R,P)-mol[k]) * sum_i { el[k][i] * |i;k> } }
  
  mol[k] - k-th copy of the nuclear system (k-th trajectory)
  el[k][i] - coefficients of i-th electronic basis state associated with k-th copy of the system
  |i;k> - is i-th electronic basis state parameterized by nuclear degrees of freedom mol[k]  

*/


  public:

  // Propagated variables
  int ntraj;                ///< number of trajectories in ensemble 
  int nnucl;                ///< number of nuclear degrees of freedom
  int nelec;                ///< number of electronic DOFs

  vector<int> is_active;    ///< flag stating if the i-th trajectory is active, if not - it is assumed to be fixed: no integration is applied
  vector<Nuclear>    mol;   ///< nuclear subsystems
  vector<Electronic>  el;   ///< electronic subsystems
  vector<Hamiltonian*> ham; ///< Hamiltonian "handlers" - unique for each copy

  // For Python access:
//  Nuclear get_mol(int i){ return *mol[i]; }
//  Electronic get_el(int i){ return *el[i]; }
//  Hamiltonian& get_ham(int i){ return *ham[i]; }

//  void set_mol(int i, Nuclear& _mol){ return mol[i] = &_mol; }
//  void set_el(int i, Electronic& _el){ return el[i] = &_el; }
//  void set_ham(int i, Hamiltonian& _ham){ ham[i] = &_ham; }
//  void set_ham(int i, Hamiltonian_Model& _ham){ ham[i] = &_ham; }

  void ham_set_ham(int i, std::string opt, int mopt);
  void ham_set_ham(std::string opt, int mopt);
  void ham_set_ham(int i, Hamiltonian& _ham);

  void ham_set_rep(int i, int _rep);
  void ham_set_rep(int _rep);

  void ham_set_params(int i, vector<double>& params_);
  void ham_set_params(int i, boost::python::list params_);
  void ham_set_params(vector<double>& params_);
  void ham_set_params(boost::python::list params_);

  void ham_set_v(int i, vector<double>& v);
  void ham_set_v(int i, boost::python::list v);
  void ham_set_v();

  void ham_set_q(int i, vector<double>& v);
  void ham_set_q(int i, boost::python::list v);

  void ham_compute(int i);
  void ham_compute_diabatic(int i);
  void ham_compute_adiabatic(int i);

  std::complex<double> ham_H(int traj, int i,int j);     
  std::complex<double> ham_dHdq(int traj, int i,int j,int n);
  std::complex<double> ham_D(int traj, int i,int j,int n);    
  std::complex<double> ham_nac(int traj,int i,int j);        
  std::complex<double> ham_Hvib(int traj, int i,int j);       



  void el_propagate_electronic(int i,double dt);
  void el_propagate_electronic(double dt);

  void mol_propagate_p(int i,double dt);
  void mol_propagate_p(double dt);

  void mol_propagate_q(int i,double dt);
  void mol_propagate_q(double dt);

  

  // Properties, countable, statistics
  // So far it is specific for 2D systems (and 1D too)
/*
  vector<double> ave_q;    // averages of coordinates q
  vector<double> ave_p;    // averages of momenta p
  vector<double> sigma_q;  // spread in q
  vector<double> sigma_p;  // spread in p
*/


  //----------- Methods -----------
  void _init(int, int, int);

  Ensemble(){ ntraj = 0; nnucl = 0; nelec = 0; }
  Ensemble(int,int,int);
 ~Ensemble();


  void se_pop(vector<double>&,double,double);
  void se_pop(vector<double>&);
  boost::python::list se_pop(double,double);
  boost::python::list se_pop();

  void sh_pop(vector<double>&,double,double);
  void sh_pop(vector<double>&);
  boost::python::list sh_pop(double,double);
  boost::python::list sh_pop();

  void sh_pop1(vector<double>&,double,double);
  void sh_pop1(vector<double>&);
  boost::python::list sh_pop1(double,double);
  boost::python::list sh_pop1();



  void print_map(std::string prefix, double Xmin, double Xmax, double dx, double Ymin, double Ymax, double dy, int snap); // for 2D projections on XY plane
  void integral_flux(vector< vector<double> >& Int_flx, double Xmin, double Xmax, double dx, double Ymin, double Ymax, double dy);

  void compute_averages();


};


}// namespace libensemble
}// namespace libdyn
}// liblibra

#endif // ENSEMBLE
