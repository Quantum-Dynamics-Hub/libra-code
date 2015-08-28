/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef NUCLEAR_H
#define NUCLEAR_H


#include "Mathematics.h"
using namespace std;


class Nuclear{
  // Nuclear system

  public:

  //------------- Data members ---------------------
  int Nnucl;                      // number of nuclei

  // Atomic properties
  vector<std::string> elt_name;   // Element symbols of atoms
  vector<std::string> at_type;    // types of atoms (whatever this means)
  vector<int> Z;                  // nuclear core charges = number in periodic table
  vector<double> Zeff;            // effective charges
  vector<double> mass;            // mass of particles

  // Dynamic variables
  vector<VECTOR> R;               // nuclear positions
  vector<VECTOR> P;               // nuclear momenta
  vector<VECTOR> grad;            // energy gradients w.r.t. nuclear coordinates
  vector<VECTOR> frcs;            // ground-state forces
  vector<double> Mull_charges_gross;    // Mulliken partial charges on all atoms
  vector<double> Mull_charges_net;      // Mulliken partial charges on all atoms


  //------------- Constructors ----------------------
  Nuclear(){ Nnucl = 0; }
  Nuclear(int _Nnucl){  

    Nnucl = _Nnucl;

    elt_name = vector<std::string>(Nnucl,"");
    at_type = vector<std::string>(Nnucl,"");

    Z = vector<int>(Nnucl,0);
    Zeff = vector<double>(Nnucl,0.0);
    mass = vector<double>(Nnucl,1.0); 
    R = vector<VECTOR>(Nnucl,VECTOR(0.0,0.0,0.0));
    P = vector<VECTOR>(Nnucl,VECTOR(0.0,0.0,0.0));
    grad = vector<VECTOR>(Nnucl,VECTOR(0.0,0.0,0.0));
    frcs = vector<VECTOR>(Nnucl,VECTOR(0.0,0.0,0.0));

    Mull_charges_gross = vector<double>(Nnucl,0.0);
    Mull_charges_net = vector<double>(Nnucl,0.0);
    
  }  

  //------------- Function members ---------------------
  int add_atom(std::string, std::string, double, double, double, VECTOR);

  VECTOR get_tv(int,int, int, int, int, const VECTOR&, const VECTOR&, const VECTOR&);

  void dipole_moment(VECTOR&, double&);
  
};

double Energy_nucl(Nuclear&,vector<int>&);


#endif // NUCLEAR_H
