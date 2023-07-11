/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef STATE_H
#define STATE_H

#include "../../dyn/libdyn.h"
#include "../../chemobjects/libchemobjects.h"
#include "../../atomistic/libatomistic.h"
#include "../../math_linalg/liblinalg.h"
#include "../../math_random/librandom.h"
#include "../../io/libio.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libdyn;
using namespace libdyn::libbarostat;
using namespace libdyn::libelectronic;
using namespace libdyn::libnuclear;


using namespace libchemobjects;
using namespace libatomistic;
using namespace liblinalg;
using namespace librandom;

namespace libscripts{
namespace libstate{


double compute_kinetic_energy(Nuclear* mol);
double compute_kinetic_energy(Nuclear& mol);

double compute_potential_energy(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt);
double compute_potential_energy(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt);

double compute_forces(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt);
double compute_forces(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt);


class MD{

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const MD&); // Copies the content which is defined
  void extract_dictionary(boost::python::dict);

public:

  std::string ensemble;         int is_ensemble;        // ensemble used (NVE, NVT, etc)    
  std::string integrator;       int is_integrator;      // name of MD integrator used

  double dt;                    int is_dt;              // time step in tau (0.02 tau = 1 fs)
  int n_medium;                 int is_n_medium;        // number of medium steps per 1 slow (dt) step - for RESPA
  int n_fast;                   int is_n_fast;          // number of fast steps per 1 medium (dt/n_medium) step - for RESPA
  int n_outer;                  int is_n_outer;         // number of outer steps(barostat,thermostat,etc.) per 1 md step (slow) - for RESPA
  double max_time;              int is_max_time;        // time of the simulation
  int    max_step;              int is_max_step;        // # of MD step

  double curr_time;             int is_curr_time;       // current time since start of MD simulation
  double curr_step;             int is_curr_step;       // current # of MD step

  int terec_exp_size;           int is_terec_exp_size;  // expansion size for Terec method
  int use_vlist;                int is_use_vlist;       // Flag turning on/off the Verlet list option 
  int vlist_upd_freq;           int is_vlist_upd_freq;  // Number of steps between updating the Verlet list
  int vlist_time;               int is_vlist_time;      // Number of steps from the last update of Verlet list

  //----------- Basic class operations ---------------------------
  // Defined in rbmd.cpp
  MD();           // constructor
  MD(boost::python::dict);
  MD(const MD&);  // copy-constructor
 ~MD();           // destructor

  MD& operator=(const MD&); // assignment operator
  void show_info();
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<MD>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<MD>& vt,int& status);



class State{

  //--------- Auxiliary internal functions -------------
  void init_variables();           // Initializes variables
  void copy_content(const State&); // Copies the content which is defined

  //
  int is_md_initialized; 

public:

  // Objects
  System*     syst;      
  Thermostat* thermostat; 
  Barostat*   barostat;     

  // Simulation parameters
  MD*  md; 

  // State veriables
  double E_kin;            int is_E_kin;           // Kinetic energy;
  double E_kin_tr;         int is_E_kin_tr;        // Translational kinetic energy
  double E_kin_rot;        int is_E_kin_rot;       // Rotational kinetic energy
  double E_pot;            int is_E_pot;           // Potential energy;
  double E_tot;            int is_E_tot;           // Total   energy;
  double H;                int is_H;               // Enthalpy             
  double H_NHC;            int is_H_NHC;           // Nose-Hoover Chains Hamiltonian.
  double H_NP;             int is_H_NP;            // Nose-Poincare Hamiltonian
  double H0;               int is_H0;              // Value of the Hamiltonian for time dependent (nonconservative) simulations
  double curr_T;           int is_curr_T;          // Current Temperature;
  double curr_V;           int is_curr_V;          // Current system volume
  double curr_P;           int is_curr_P;          // Current system pressure
  MATRIX3x3 curr_P_tens;   int is_curr_P_tens;     // Current pressure tensor
  int Nf_t;                int is_Nf_t;            // Number of translational degrees of freedom
  int Nf_r;                int is_Nf_r;            // Number of rotational degrees of freedom
  VECTOR L_tot;            int is_L_tot;           // Total angular momentum
  VECTOR P_tot;            int is_P_tot;           // Total linear momentum

  //----------- Basic class operations ---------------------------
  // Defined in State.cpp
  State();              // constructor
  void set_system(System& syst_){ syst = &syst_; }
  void set_thermostat(Thermostat& therm){ thermostat = &therm; }
  void set_barostat(Barostat& baro){ barostat = &baro; }
  void set_md(MD& md_){ md = &md_; }

  State(const State&);  // copy-constructor
 ~State();              // destructor
  State& operator=(const State&); // assignment operator
  void show_info();
  void set(object);

  //----------------------------------------------------------------
  // Defined in State_methods.cpp
  void update();
  void cool();

  // Defined in State_methods1.cpp
  void init_md(Nuclear& mol, Electronic& el, Hamiltonian& ham, Random& rnd);
  void run_md(Nuclear& mol, Electronic& el, Hamiltonian& ham);

  // Defined in State_methods2.cpp
  void init_md(Electronic& el, Hamiltonian& ham, Random& rnd);
  void run_md(Electronic& el, Hamiltonian& ham);


  // Defined in State_methods2.cpp
//  void run_md_respa();

};


}// namespace libstate
}// namespace libscripts

}// liblbira

#endif // STATE_H
