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

#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H

#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <complex>
#include <string.h>
#include "Mathematics.h"
#include "units.h"
#include "Control_Parameters.h"
#include "MOAO.h"
#include "Nuclear.h"

using namespace std;

class HF_integrals{

  // HF - Coulomb and exchange integrals in AO basis
  struct data_element{
    int a, b, c, d;
    double J_abcd;
    double K_adcb;
  };

  vector<data_element> data;

  int find_data(int,int,int,int);  


  public:

  HF_integrals(){ ; ; }

  void set_JK_values(int,int,int,int,double, double);  
  void get_JK_values(int,int,int,int,double&,double&);  

};


class EHT_K{
// This data type realizes a container for a bunch of K constants for EHT
// with some access methods

  struct data_element{
    // Example:   H  1s   Si  3s   1.34   K1_value
    std::string elt1;
    std::string orb_type1;
    std::string elt2;
    std::string orb_type2;

    // H_ij = 0.5*K*(H_ii+H_jj)*S_ij + (K1 + K2*S_ij/sqrt(1+R^2))
    double K_value;
    double K1_value; // units of energy
    double K2_value; // units of energy
    double K3_value; // units of energy
    double K4_value; // units of lengths

    double C0_value; // units of energy
    double C1_value; // units of energy/length
    double C2_value; // units of energy/length^2
    double C3_value; // units of energy/length^3
    double C4_value; // units of lengths

  };

  struct pp_data_element{

    std::string elt1;
    std::string orb_type1;

    double PPa_value; // units of 1/length
    double PP0_value; // units of energy
    double PP1_value; // units of energy/length
    double PP2_value; // units of energy/length^2

  };

  struct psps_data_element{

    int n1,n2,n3,n4;
    double K;
  };



  int find_data(std::string,std::string);  // in data
  int find_data(std::string,std::string,std::string,std::string); // in PP_data
  int find_data(int, int, int, int); // in PSPS_data
 

public:

  vector<data_element> data;  
  vector<pp_data_element> pp_data;  
  vector<psps_data_element> psps_data;


  EHT_K(){ ;; }

  void set_PSPS_value(int, int, int, int, double);

  void set_PPa_value(std::string, std::string, double ); 
  void set_PP0_value(std::string, std::string, double ); 
  void set_PP1_value(std::string, std::string, double ); 
  void set_PP2_value(std::string, std::string, double ); 

  double get_PPa_value(std::string, std::string ); 
  double get_PP0_value(std::string, std::string ); 
  double get_PP1_value(std::string, std::string ); 
  double get_PP2_value(std::string, std::string ); 



  void set_K_value(std::string,std::string,std::string,std::string, double);
  double get_K_value(std::string,std::string,std::string,std::string);

  void set_K1_value(std::string,std::string,std::string,std::string, double);
  double get_K1_value(std::string,std::string,std::string,std::string);

  void set_K2_value(std::string,std::string,std::string,std::string, double);
  double get_K2_value(std::string,std::string,std::string,std::string);

  void set_K3_value(std::string,std::string,std::string,std::string, double);
  double get_K3_value(std::string,std::string,std::string,std::string);

  void set_K4_value(std::string,std::string,std::string,std::string, double);
  double get_K4_value(std::string,std::string,std::string,std::string);



  void set_C0_value(std::string,std::string,std::string,std::string, double);
  double get_C0_value(std::string,std::string,std::string,std::string);

  void set_C1_value(std::string,std::string,std::string,std::string, double);
  double get_C1_value(std::string,std::string,std::string,std::string);

  void set_C2_value(std::string,std::string,std::string,std::string, double);
  double get_C2_value(std::string,std::string,std::string,std::string);

  void set_C3_value(std::string,std::string,std::string,std::string, double);
  double get_C3_value(std::string,std::string,std::string,std::string);

  void set_C4_value(std::string,std::string,std::string,std::string, double);
  double get_C4_value(std::string,std::string,std::string,std::string);





};


class mEHT_K{

public:
  int size; // is the size of the matrix (N below - number of AOs in the basis)
  vector<double> eht_K;
  vector<double> eht_K1;   // eht_K1[I*N+J] = is the parameter for pair of orbitals I and J - each of special type and on special atom
  vector<double> eht_K2;
  vector<double> eht_K3;
  vector<double> eht_K4;

  vector<double> eht_C0;
  vector<double> eht_C1;   // eht_C1[I*N+J] = is the parameter for pair of orbitals I and J - each of special type and on special atom
  vector<double> eht_C2;
  vector<double> eht_C3;
  vector<double> eht_C4;

    
  // Atomic-orbital pseudopotential variables
  vector<vector<double> >  eht_PPa;    // Nnucl x Ntyp(this depends on nucleus)
  vector<vector<double> >  eht_PP0;
  vector<vector<double> >  eht_PP1;
  vector<vector<double> >  eht_PP2;




  mEHT_K(){  size = -1; }
  void set_mapping(EHT_K& k, const vector<AO>& basis);
  void set_mapping1(EHT_K& k, Nuclear& mol);

  inline double get_K_value( int I,int J){ return eht_K[I*size+J];  }
  inline double get_K1_value(int I,int J){ return eht_K1[I*size+J]; }
  inline double get_K2_value(int I,int J){ return eht_K2[I*size+J]; }
  inline double get_K3_value(int I,int J){ return eht_K3[I*size+J]; }
  inline double get_K4_value(int I,int J){ return eht_K4[I*size+J]; }

  inline double get_C0_value(int I,int J){ return eht_C0[I*size+J];  }
  inline double get_C1_value(int I,int J){ return eht_C1[I*size+J]; }
  inline double get_C2_value(int I,int J){ return eht_C2[I*size+J]; }
  inline double get_C3_value(int I,int J){ return eht_C3[I*size+J]; }
  inline double get_C4_value(int I,int J){ return eht_C4[I*size+J]; }


};


class Element{

  public:

  // General atomic properties
  std::string elt_name;   // element symbol  
  int Z;                  // nuclear charge
  int PQN;                // principal quantum number
  int Nval;               // effective number of electrons(e.g. valence electrons)  
  double Zeff;            // effective(screened) nuclear charge
  double mass;            // mass


  // EHT, CNDO, CNDO/2, INDO
  map<std::string,double> IP;   // Ionization potential for a given shell (first = 1s, 2s, 2p, ... 3d, so on)
  map<std::string,double> EA;   // Electron affinity for a given shell (first = 1s, 2s, 2p, ... 3d, so on)  


  map<std::string, int > Nquant;    // principal quantum number for each shell
  map<std::string, int> Nzeta;   // number of different zetas (exponents for given shell)
  map<std::string, vector<double> > zetas;  // zetas for each shell 
  map<std::string, vector<double> > coeffs; // coefficients for each shell

  map<std::string,double> J_param1; // in .../sqrt(r^2 + f(J_param1))
  map<std::string,double> J_param2; // in f'(J_param2)/sqrt(...)
  map<std::string,double> J_param3; // in .../sqrt(r^2 + f(J_param1))  - for el-nucl interaction
  map<std::string,double> J_param4; // in f'(J_param2)/sqrt(...)       - for el-nucl interaction



  // CNDO, CNDO/2, INDO
  map<std::string,double> G1,F2;    // Slater-Condon factors
  map<std::string,double> beta0;    // proportionality constants for off-diagonal matrix elements
  double Zeta;                  // exponent of the radial part of the STOs

  
  Element(){}
  Element(std::string en,int z,int nv,map<std::string, double> ip){ elt_name = en; Z = z; Nval = nv; IP = ip; }

  ~Element(){ 
     if(IP.size()>0) { IP.clear();}
     if(EA.size()>0) { EA.clear();}
   }

  void _set(std::string en,int z){ elt_name = en; Z = z;  }
  void _set(std::string en,int z,int nv){ elt_name = en; Z = z; Nval = nv; }
  void _set(std::string en,int z,int nv,double zeff){ elt_name = en; Z = z; Nval = nv; Zeff = zeff; }
  void _set(std::string en,int z,int nv,double zeff,map<std::string, double> ip){ elt_name = en; Z = z; Nval = nv; IP = ip; Zeff = zeff;}

  void set_mass(double m_){ mass = m_; }

};



class OrbParams{

// Here we just keep parameters for a single orbital type (e.g. H(1s), Si(3p), etc.)

  public:

  // General atomic properties
  double IP, EA;
  int Nquant;             // principal quantum number for each shell
  int Nzeta;              // number of different zetas (exponents for given shell)
  vector<double> zetas;   // zetas for each shell 
  vector<double> coeffs;  // coefficients for each shell

  double J_param1; 
  double J_param2; 
  double J_param3; 
  double J_param4; 


  // CNDO, CNDO/2, INDO
  double G1,F2;          // Slater-Condon factors
  double beta0;          // proportionality constants for off-diagonal matrix elements

  
  OrbParams(){}

  ~OrbParams(){ 
     if(zetas.size()>0) { zetas.clear();}
     if(coeffs.size()>0) { coeffs.clear();}
  }


};




void set_default_elements(map<std::string,Element>&);


class Model_Parameters{

public:
//---------- Members --------

  // Parameters for semiempirical methods
  // Atomic parameters
  map<std::string,Element> PT;    // General and specific (for a given method) atomic parameters
  vector<OrbParams> orb_params;   // properties all orbitals!

  // Pair parameters
  EHT_K  eht_k;                   // K values for EHT

  // Mapped variables - for greater efficiency
  mEHT_K meht_k;                  // mapped EHT parameters

  


  HF_integrals  hf_int;           // precomputed J and K integrals in given basis

//  vector<double> eri;    // precomputed electron repulsion integrals - for all atoms
//  vector<> V_A

  
  //-------------- Constructor --------------
  Model_Parameters(){  
    set_default_elements(PT);
  }

  void set_PT_mapping(const vector<AO>&);

};


void set_parameters_eht(Control_Parameters&, Model_Parameters&);
void set_parameters_hf(Control_Parameters&, Model_Parameters&, vector<AO>&);
void set_parameters_indo(Control_Parameters&, Model_Parameters&);
void set_parameters_geht1(Control_Parameters& prms, Model_Parameters& modprms); // Model_Parameters_GEHT.cpp
void set_parameters_geht2(Control_Parameters& prms, Model_Parameters& modprms); // Model_Parameters_GEHT2.cpp

void set_parameters_eht_mapping(Model_Parameters& modprms,const vector<AO>& basis_ao);
void set_parameters_eht_mapping1(Model_Parameters& modprms, Nuclear& mol);




#endif // MODEL_PARAMETERS_H