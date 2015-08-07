#ifndef STATE_H
#define STATE_H

#include "System.h"
#include "Thermostat.h"
#include "Barostat.h"

#include "MD.h"

class State{

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const State&); // Copies the content which is defined

  //
  int is_md_initialized; 

public:

  // Objects
  System*     system;      
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
  State(System&);       
  State(System&,Thermostat&);
  State(System&,Barostat&);
  State(System&,Thermostat&,Barostat&);

  State(const State&);  // copy-constructor
 ~State();              // destructor

  State& operator=(const State&); // assignment operator
  void show_info();
  void set(object);

  //----------------------------------------------------------------
  // Defined in State_methods.cpp
  void update();
  void cool();
  void init_velocities(double);
  void init_velocities(double,VECTOR,VECTOR);

  // Defined in State_methods1.cpp
  void set_md(MD&);
  void init_md();
  void run_md();

  // Defined in State_methods2.cpp
  void run_md_respa();

};

#endif // STATE_H
