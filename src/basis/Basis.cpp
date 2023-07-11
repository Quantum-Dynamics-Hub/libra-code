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
  \file Basis.cpp
  \brief The file implements functions for creating basis atomic orbitals from molecular structure information
    
*/

#include "Basis.h"

/// liblibra namespace
namespace liblibra{

/// libbasis namespace
namespace libbasis{


void basis_params_s(int Nquant, vector<double>& alp, vector<double>& coeff){
/**
  \brief An auxiliary function that returns parameters of STO-3G s-type orbitals
  \param[in] Nquant The principal quantum number of the atom for which we construct the AO
  \param[out] alp The exponents of the primitive Gaussians in the STO-3G expansion
  \param[out] coeff The contraction coefficients the primitive Gaussians enter with in the STO-3G expansion

  The function takes expansion coefficients to be equal to the expansion coefficients
  fron the STO-3G basis in the valence shell. The alphas are obtained from a scaling law.
  Data are taken from EMSL website (mostly) or extrapolated

  Nquant = 1 - H exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  Nquant = 2 - C exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  Nquant = 3 - Si exp is calibrated to zeta = 1.383
  Nquant = 4 - Ge exp is calibrated to zeta = 2.160
  Nquant = 5 - Sn exp is calibrated to zeta = 2.120
  There is no STO-3G bases for 6-row elements, so will need to extrapolate
  so far - just keep it the same as for 5-row

*/


  if(Nquant==1){  // H exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
    alp[0] = 2.227660;   coeff[0] = 0.154329;
    alp[1] = 0.405771;   coeff[1] = 0.535328;
    alp[2] = 0.109818;   coeff[2] = 0.444635;
  }// n = 1
  else if(Nquant==2){ // C  exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
    alp[0] = 0.994203 ;   coeff[0] = -0.09996723;
    alp[1] = 0.231031 ;   coeff[1] =  0.39951283;
    alp[2] = 0.0751386;   coeff[2] =  0.70011547;

  }
  else if(Nquant==3){ // Si  exp is calibrated to zeta = 1.383
    alp[0] = 0.773121308;   coeff[0] = -0.2196203690;
    alp[1] = 0.215698883;   coeff[1] =  0.2255954336;
    alp[2] = 0.084423081;   coeff[2] =  0.9003984260;

  }
  else if(Nquant==4){ // Ge  exp is calibrated to zeta = 2.160
    alp[0] = 0.211298131;   coeff[0] = -0.3088441215;
    alp[1] = 0.077982299;   coeff[1] =  0.0196064117;
    alp[2] = 0.034437806;   coeff[2] =  1.1310344420;

  }
  else if(Nquant==5){  // Sn exp is calibrated to zeta = 2.120
    alp[0] = 0.138746360;   coeff[0] = -0.3842642607;
    alp[1] = 0.074706337;   coeff[1] = -0.1972567438;
    alp[2] = 0.032999103;   coeff[2] =  1.3754955120;

  }
  else if(Nquant==6){  // Pb
    // There is no STO-3G bases for 6-row elements, so will need to extrapolate
    // so far - just keep it the same as for 5-row

    alp[0] = 0.138746360;   coeff[0] = -0.3842642607;
    alp[1] = 0.074706337;   coeff[1] = -0.1972567438;
    alp[2] = 0.032999103;   coeff[2] =  1.3754955120;

  }


}

void basis_params_p(int Nquant, vector<double>& alp, vector<double>& coeff){
/**
  \brief An auxiliary function that returns parameters of STO-3G p-type orbitals
  \param[in] Nquant The principal quantum number of the atom for which we construct the AO
  \param[out] alp The exponents of the primitive Gaussians in the STO-3G expansion
  \param[out] coeff The contraction coefficients the primitive Gaussians enter with in the STO-3G expansion

  The function takes expansion coefficients to be equal to the expansion coefficients
  fron the STO-3G basis in the valence shell. The alphas are obtained from a scaling law.
  Data are taken from EMSL website (mostly) or extrapolated

  Nquant = 1 - H exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  Nquant = 2 - C exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  Nquant = 3 - Si exp is calibrated to zeta = 1.383
  Nquant = 4 - Ge exp is calibrated to zeta = 2.160
  Nquant = 5 - Sn exp is calibrated to zeta = 2.120
  There is no STO-3G bases for 6-row elements, so will need to extrapolate
  so far - just keep it the same as for 5-row

  The exponents are the same as in basis_params_s, only the contraction coefficients are different
*/

  if(Nquant==1){  // H exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
/*
    alp[0] = 2.227660;   coeff[0] = 0.154329;
    alp[1] = 0.405771;   coeff[1] = 0.535328;
    alp[2] = 0.109818;   coeff[2] = 0.444635;
*/
  }// n = 1
  else if(Nquant==2){ // C  exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
    alp[0] = 0.994203 ;   coeff[0] = 0.15591627;
    alp[1] = 0.231031 ;   coeff[1] = 0.60768372;
    alp[2] = 0.0751386;   coeff[2] = 0.39195739;

  }
  else if(Nquant==3){ // Si  exp is calibrated to zeta = 1.383
    alp[0] = 0.773121308;   coeff[0] = 0.01058760429;
    alp[1] = 0.215698883;   coeff[1] = 0.59516700530;
    alp[2] = 0.084423081;   coeff[2] = 0.46200101200;

  }
  else if(Nquant==4){ // Ge  exp is calibrated to zeta = 2.160
    alp[0] = 0.211298131;   coeff[0] = -0.1215468600;
    alp[1] = 0.077982299;   coeff[1] =  0.5715227604;
    alp[2] = 0.034437806;   coeff[2] =  0.5498949471;

  }
  else if(Nquant==5){  // Sn exp is calibrated to zeta = 2.120
    alp[0] = 0.138746360;   coeff[0] = -0.3481691526;
    alp[1] = 0.074706337;   coeff[1] =  0.6290323690;
    alp[2] = 0.032999103;   coeff[2] =  0.6662832743;

  }
  else if(Nquant==6){  // Pb
    // There is no STO-3G bases for 6-row elements, so will need to extrapolate
    // so far - just keep it the same as for 5-row

    alp[0] = 0.138746360;   coeff[0] = -0.3481691526;
    alp[1] = 0.074706337;   coeff[1] =  0.6290323690;
    alp[2] = 0.032999103;   coeff[2] =  0.6662832743;

  }

}

void basis_params_d(int Nquant, vector<double>& alp, vector<double>& coeff){
/**
  \brief An auxiliary function that returns parameters of STO-3G d-type orbitals
  \param[in] Nquant The principal quantum number of the atom for which we construct the AO
  \param[out] alp The exponents of the primitive Gaussians in the STO-3G expansion
  \param[out] coeff The contraction coefficients the primitive Gaussians enter with in the STO-3G expansion

  The function takes expansion coefficients to be equal to the expansion coefficients
  fron the STO-3G basis in the valence shell. The alphas are obtained from a scaling law.
  Data are taken from EMSL website (mostly) or extrapolated

  Nquant = 1 - H exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  Nquant = 2 - C exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  Nquant = 3 - taken to be the same as for Nquant = 4
  Nquant = 4 - Sc  exp is calibrated to zeta = 1.7 (exponent with the largest weight)
  Nquant = 5 - Nb exp is calibrated to zeta = 1.64
  There is no STO-3G bases for 6-row elements, so will need to extrapolate
  so far - just keep it the same as for 5-row

*/


  if(Nquant==1){  // H exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  }// n = 1
  else if(Nquant==2){ // C  exp is taken directly from Pople: Hehre, Stewart, Pople  JCP 1969, 51, 2657
  }
  else if(Nquant==3){ // taken to be the same as for Nquant = 4
    alp[0] = 0.190899677;   coeff[0] = 0.2197679508;
    alp[1] = 0.058230486;   coeff[1] = 0.6555473627;
    alp[2] = 0.022467132;   coeff[2] = 0.2865732590;

  }
  else if(Nquant==4){ // Sc  exp is calibrated to zeta = 1.7 (exponent with the largest weight)
    alp[0] = 0.190899677;   coeff[0] = 0.2197679508;
    alp[1] = 0.058230486;   coeff[1] = 0.6555473627;
    alp[2] = 0.022467132;   coeff[2] = 0.2865732590;

  }
  else if(Nquant==5){  // Nb exp is calibrated to zeta = 1.64
    alp[0] = 0.500029323;   coeff[0] = 0.1250662138;
    alp[1] = 0.194708826;   coeff[1] = 0.6686785577;
    alp[2] = 0.085711305;   coeff[2] = 0.3052468245;

  }
  else if(Nquant==6){  // 
    // There is no STO-3G bases for 6-row elements, so will need to extrapolate
    // so far - just keep it the same as for 5-row
    alp[0] = 0.500029323;   coeff[0] = 0.1250662138;
    alp[1] = 0.194708826;   coeff[1] = 0.6686785577;
    alp[2] = 0.085711305;   coeff[2] = 0.3052468245;

  }



}


void add_basis_ao(std::string Atom_name, VECTOR& R, std::string Atom_shell, int Nzeta, int Nquant,
                    double  IP, double exp1, double exp2, double coeff1, double coeff2,
                    vector<AO>& basis_ao){
/**
  \brief Create a new AO (STO-3G) object and add it to existing (including empty) list of orbitals (basis)
  \param[in] Atom_name The name of the atom (element) for which we want to create AO 
  \param[in] R The coordinate at which to create the atomic orbital. 
  \param[in] Atom_shell The name of the AO shell e.g. 1s, 2p, 3p, 3d, etc.
  \param[in] Nzeta The number of different Slater-type zetas - in most cases 1, for d-orbitals can be 2
  \param[in] Nquant The principal quantum number of the orbital
  \param[in] IP The state-specific ionization potential (not actually used in this function!)
  \param[in] exp1 The first Slater-type exponent
  \param[in] exp2 The second Slater-type exponent - is not used when Nzeta = 1
  \param[in] coeff1 The contraction coefficient of the first Slater zeta orbital
  \param[in] coeff2 The contraction coefficient of the second Slater zeta orbital (is effectively = 0 if Nzeta = 1)
  \param[in,out] basis_ao The list of AO objects - the basis we are creating

  Sequential application of this function with proper parameters will generate the atomic basis - the list of AO objects

  For conversion of the spherical hamonics to Cartesian coordinates see e.g.:
  http://csi.chemie.tu-darmstadt.de/ak/immel/script/redirect.cgi?filename=http://csi.chemie.tu-darmstadt.de/ak/immel/tutorials/orbitals/hydrogenic.html
  So the z2 orbital is actually: (1/(?)) *[ 2*z2 - (x2+x2) ]  
  and the x2-y2 orbital is (1/sqrt(2)) * [x2 - y2]

*/

  int basis_size = 0;


  double a1, a2, a3;
  PrimitiveG s;
  AO ao;  
  int Norb = 0;
//  boost::python::list res;

  stringstream ss(stringstream::in | stringstream::out);
  std::string out;
  (ss << Nquant);  ss >> out; 

  vector<double> s_alp(3,0.0);       vector<double> s_coeff(3,0.0);
  vector<double> p_alp(3,0.0);       vector<double> p_coeff(3,0.0);
  vector<double> d_alp(3,0.0);       vector<double> d_coeff(3,0.0);

  basis_params_s(Nquant, s_alp, s_coeff);
  basis_params_p(Nquant, p_alp, p_coeff);
  basis_params_d(Nquant, d_alp, d_coeff);


  if(Nzeta==1){ // single-zeta basis function

    if(Atom_shell=="s"){
      s.init(0,0,0, s_alp[0]*exp1*exp1, R);  ao.add_primitive( s_coeff[0]*coeff1, s);
      s.init(0,0,0, s_alp[1]*exp1*exp1, R);  ao.add_primitive( s_coeff[1]*coeff1, s);
      s.init(0,0,0, s_alp[2]*exp1*exp1, R);  ao.add_primitive( s_coeff[2]*coeff1, s);
                            
      ao.ao_name = out+Atom_shell;          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 0;                         ao.is_x_exp = 1;
      ao.y_exp = 0;                         ao.is_y_exp = 1;
      ao.z_exp = 0;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;

    }// Atom_shell=="s"
    else if(Atom_shell=="p"){
      s.init(1,0,0, p_alp[0]*exp1*exp1, R);  ao.add_primitive( p_coeff[0]*coeff1, s);
      s.init(1,0,0, p_alp[1]*exp1*exp1, R);  ao.add_primitive( p_coeff[1]*coeff1, s);
      s.init(1,0,0, p_alp[2]*exp1*exp1, R);  ao.add_primitive( p_coeff[2]*coeff1, s);
      
      ao.ao_name = out+Atom_shell+"x";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 1;                         ao.is_x_exp = 1;
      ao.y_exp = 0;                         ao.is_y_exp = 1;
      ao.z_exp = 0;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(0,1,0, p_alp[0]*exp1*exp1, R);  ao.add_primitive( p_coeff[0]*coeff1, s);
      s.init(0,1,0, p_alp[1]*exp1*exp1, R);  ao.add_primitive( p_coeff[1]*coeff1, s);
      s.init(0,1,0, p_alp[2]*exp1*exp1, R);  ao.add_primitive( p_coeff[2]*coeff1, s);
      
      ao.ao_name = out+Atom_shell+"y";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 0;                         ao.is_x_exp = 1;
      ao.y_exp = 1;                         ao.is_y_exp = 1;
      ao.z_exp = 0;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(0,0,1, p_alp[0]*exp1*exp1, R);  ao.add_primitive( p_coeff[0]*coeff1, s);
      s.init(0,0,1, p_alp[1]*exp1*exp1, R);  ao.add_primitive( p_coeff[1]*coeff1, s);
      s.init(0,0,1, p_alp[2]*exp1*exp1, R);  ao.add_primitive( p_coeff[2]*coeff1, s);
      
      ao.ao_name = out+Atom_shell+"z";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 0;                         ao.is_x_exp = 1;
      ao.y_exp = 0;                         ao.is_y_exp = 1;
      ao.z_exp = 1;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;

    }// Atom_shell=="p"

    else if(Atom_shell=="d"){

      double nrm = 1.0/sqrt(2.0);
  
      s.init(2,0,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive( nrm*d_coeff[0]*coeff1, s);
      s.init(2,0,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive( nrm*d_coeff[1]*coeff1, s);
      s.init(2,0,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive( nrm*d_coeff[2]*coeff1, s);
      s.init(0,2,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[0]*coeff1, s);
      s.init(0,2,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[1]*coeff1, s);
      s.init(0,2,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[2]*coeff1, s);

      
      ao.ao_name = out+Atom_shell+"x2-y2";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;



      // 2*z2 - (x2+y2)
      s.init(0,0,2, d_alp[0]*exp1*exp1, R);  ao.add_primitive(2.0*d_coeff[0]*coeff1, s);
      s.init(0,0,2, d_alp[1]*exp1*exp1, R);  ao.add_primitive(2.0*d_coeff[1]*coeff1, s);
      s.init(0,0,2, d_alp[2]*exp1*exp1, R);  ao.add_primitive(2.0*d_coeff[2]*coeff1, s);

      s.init(2,0,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive(-d_coeff[0]*coeff1, s);
      s.init(2,0,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive(-d_coeff[1]*coeff1, s);
      s.init(2,0,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive(-d_coeff[2]*coeff1, s);

      s.init(0,2,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive(-d_coeff[0]*coeff1, s);
      s.init(0,2,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive(-d_coeff[1]*coeff1, s);
      s.init(0,2,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive(-d_coeff[2]*coeff1, s);

      
      ao.ao_name = out+Atom_shell+"z2";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(1,1,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive( d_coeff[0]*coeff1, s);
      s.init(1,1,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive( d_coeff[1]*coeff1, s);
      s.init(1,1,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive( d_coeff[2]*coeff1, s);
      
      ao.ao_name = out+Atom_shell+"xy";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(1,0,1, d_alp[0]*exp1*exp1, R);  ao.add_primitive( d_coeff[0]*coeff1, s);
      s.init(1,0,1, d_alp[1]*exp1*exp1, R);  ao.add_primitive( d_coeff[1]*coeff1, s);
      s.init(1,0,1, d_alp[2]*exp1*exp1, R);  ao.add_primitive( d_coeff[2]*coeff1, s);
      
      ao.ao_name = out+Atom_shell+"xz";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(0,1,1, d_alp[0]*exp1*exp1, R);  ao.add_primitive( d_coeff[0]*coeff1, s);
      s.init(0,1,1, d_alp[1]*exp1*exp1, R);  ao.add_primitive( d_coeff[1]*coeff1, s);
      s.init(0,1,1, d_alp[2]*exp1*exp1, R);  ao.add_primitive( d_coeff[2]*coeff1, s);
      
      ao.ao_name = out+Atom_shell+"yz";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


    }// Atom_shell=="d"



  }// Nzeta == 1

  else if(Nzeta==2){ // double-zeta basis function

    if(Atom_shell=="s"){
      s.init(0,0,0, s_alp[0]*exp1*exp1, R);  ao.add_primitive( s_coeff[0]*coeff1, s);
      s.init(0,0,0, s_alp[1]*exp1*exp1, R);  ao.add_primitive( s_coeff[1]*coeff1, s);
      s.init(0,0,0, s_alp[2]*exp1*exp1, R);  ao.add_primitive( s_coeff[2]*coeff1, s);
      s.init(0,0,0, s_alp[0]*exp2*exp2, R);  ao.add_primitive( s_coeff[0]*coeff2, s);
      s.init(0,0,0, s_alp[1]*exp2*exp2, R);  ao.add_primitive( s_coeff[1]*coeff2, s);
      s.init(0,0,0, s_alp[2]*exp2*exp2, R);  ao.add_primitive( s_coeff[2]*coeff2, s);
      
      ao.ao_name = out+Atom_shell;          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 0;                         ao.is_x_exp = 1;
      ao.y_exp = 0;                         ao.is_y_exp = 1;
      ao.z_exp = 0;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;

    }// Atom_shell=="s"
    else if(Atom_shell=="p"){
      s.init(1,0,0, p_alp[0]*exp1*exp1, R);  ao.add_primitive( p_coeff[0]*coeff1, s);
      s.init(1,0,0, p_alp[1]*exp1*exp1, R);  ao.add_primitive( p_coeff[1]*coeff1, s);
      s.init(1,0,0, p_alp[2]*exp1*exp1, R);  ao.add_primitive( p_coeff[2]*coeff1, s);
      s.init(1,0,0, p_alp[0]*exp2*exp2, R);  ao.add_primitive( p_coeff[0]*coeff2, s);
      s.init(1,0,0, p_alp[1]*exp2*exp2, R);  ao.add_primitive( p_coeff[1]*coeff2, s);
      s.init(1,0,0, p_alp[2]*exp2*exp2, R);  ao.add_primitive( p_coeff[2]*coeff2, s);
      
      ao.ao_name = out+Atom_shell+"x";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 1;                         ao.is_x_exp = 1;
      ao.y_exp = 0;                         ao.is_y_exp = 1;
      ao.z_exp = 0;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(0,1,0, p_alp[0]*exp1*exp1, R);  ao.add_primitive( p_coeff[0]*coeff1, s);
      s.init(0,1,0, p_alp[1]*exp1*exp1, R);  ao.add_primitive( p_coeff[1]*coeff1, s);
      s.init(0,1,0, p_alp[2]*exp1*exp1, R);  ao.add_primitive( p_coeff[2]*coeff1, s);
      s.init(0,1,0, p_alp[0]*exp2*exp2, R);  ao.add_primitive( p_coeff[0]*coeff2, s);
      s.init(0,1,0, p_alp[1]*exp2*exp2, R);  ao.add_primitive( p_coeff[1]*coeff2, s);
      s.init(0,1,0, p_alp[2]*exp2*exp2, R);  ao.add_primitive( p_coeff[2]*coeff2, s);
      
      ao.ao_name = out+Atom_shell+"y";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 0;                         ao.is_x_exp = 1;
      ao.y_exp = 1;                         ao.is_y_exp = 1;
      ao.z_exp = 0;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(0,0,1, p_alp[0]*exp1*exp1, R);  ao.add_primitive( p_coeff[0]*coeff1, s);
      s.init(0,0,1, p_alp[1]*exp1*exp1, R);  ao.add_primitive( p_coeff[1]*coeff1, s);
      s.init(0,0,1, p_alp[2]*exp1*exp1, R);  ao.add_primitive( p_coeff[2]*coeff1, s);
      s.init(0,0,1, p_alp[0]*exp2*exp2, R);  ao.add_primitive( p_coeff[0]*coeff2, s);
      s.init(0,0,1, p_alp[1]*exp2*exp2, R);  ao.add_primitive( p_coeff[1]*coeff2, s);
      s.init(0,0,1, p_alp[2]*exp2*exp2, R);  ao.add_primitive( p_coeff[2]*coeff2, s);
      
      ao.ao_name = out+Atom_shell+"z";          ao.is_ao_name = 1;      
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.
      ao.x_exp = 0;                         ao.is_x_exp = 1;
      ao.y_exp = 0;                         ao.is_y_exp = 1;
      ao.z_exp = 1;                         ao.is_z_exp = 1;
      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;

    }// Atom_shell=="p"
    else if(Atom_shell=="d"){

/*  For conversion of the spherical hamonics to Cartesian coordinates see e.g.:

   http://csi.chemie.tu-darmstadt.de/ak/immel/script/redirect.cgi?filename=http://csi.chemie.tu-darmstadt.de/ak/immel/tutorials/orbitals/hydrogenic.html

   So the z2 orbital is actually: (1/(?)) *[ 2*z2 - (x2+x2) ]  ?=6??
   and the x2-y2 orbital is (1/sqrt(2)) * [x2 - y2]
*/


      double nrm = 1.0/sqrt(2.0);
  
      s.init(2,0,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive( nrm*d_coeff[0]*coeff1, s);
      s.init(2,0,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive( nrm*d_coeff[1]*coeff1, s);
      s.init(2,0,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive( nrm*d_coeff[2]*coeff1, s);
      s.init(0,2,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[0]*coeff1, s);
      s.init(0,2,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[1]*coeff1, s);
      s.init(0,2,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[2]*coeff1, s);

      s.init(2,0,0, d_alp[0]*exp2*exp2, R);  ao.add_primitive( nrm*d_coeff[0]*coeff2, s);
      s.init(2,0,0, d_alp[1]*exp2*exp2, R);  ao.add_primitive( nrm*d_coeff[1]*coeff2, s);
      s.init(2,0,0, d_alp[2]*exp2*exp2, R);  ao.add_primitive( nrm*d_coeff[2]*coeff2, s);
      s.init(0,2,0, d_alp[0]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[0]*coeff2, s);
      s.init(0,2,0, d_alp[1]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[1]*coeff2, s);
      s.init(0,2,0, d_alp[2]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[2]*coeff2, s);

      
      ao.ao_name = out+Atom_shell+"x2-y2";          ao.is_ao_name = 1;      
//      ao.at_indx = Atom_index;              ao.is_at_indx = 1;       // 1, 2, ...
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;




      nrm = 1.0; // in fact, this normalization constant does not matter, since the AO will be normalized 
                 // once all primitives are collected


      s.init(0,0,2, d_alp[0]*exp1*exp1, R);  ao.add_primitive(2.0*nrm*d_coeff[0]*coeff1, s);
      s.init(0,0,2, d_alp[1]*exp1*exp1, R);  ao.add_primitive(2.0*nrm*d_coeff[1]*coeff1, s);
      s.init(0,0,2, d_alp[2]*exp1*exp1, R);  ao.add_primitive(2.0*nrm*d_coeff[2]*coeff1, s);
      s.init(0,0,2, d_alp[0]*exp2*exp2, R);  ao.add_primitive(2.0*nrm*d_coeff[0]*coeff2, s);
      s.init(0,0,2, d_alp[1]*exp2*exp2, R);  ao.add_primitive(2.0*nrm*d_coeff[1]*coeff2, s);
      s.init(0,0,2, d_alp[2]*exp2*exp2, R);  ao.add_primitive(2.0*nrm*d_coeff[2]*coeff2, s);

      s.init(2,0,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[0]*coeff1, s);
      s.init(2,0,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[1]*coeff1, s);
      s.init(2,0,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[2]*coeff1, s);
      s.init(2,0,0, d_alp[0]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[0]*coeff2, s);
      s.init(2,0,0, d_alp[1]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[1]*coeff2, s);
      s.init(2,0,0, d_alp[2]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[2]*coeff2, s);

      s.init(0,2,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[0]*coeff1, s);
      s.init(0,2,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[1]*coeff1, s);
      s.init(0,2,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive(-nrm*d_coeff[2]*coeff1, s);
      s.init(0,2,0, d_alp[0]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[0]*coeff2, s);
      s.init(0,2,0, d_alp[1]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[1]*coeff2, s);
      s.init(0,2,0, d_alp[2]*exp2*exp2, R);  ao.add_primitive(-nrm*d_coeff[2]*coeff2, s);

      
      ao.ao_name = out+Atom_shell+"z2";          ao.is_ao_name = 1;
//      ao.at_indx = Atom_index;              ao.is_at_indx = 1;       // 1, 2, ...
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(1,1,0, d_alp[0]*exp1*exp1, R);  ao.add_primitive( d_coeff[0]*coeff1, s);
      s.init(1,1,0, d_alp[1]*exp1*exp1, R);  ao.add_primitive( d_coeff[1]*coeff1, s);
      s.init(1,1,0, d_alp[2]*exp1*exp1, R);  ao.add_primitive( d_coeff[2]*coeff1, s);
      s.init(1,1,0, d_alp[0]*exp2*exp2, R);  ao.add_primitive( d_coeff[0]*coeff2, s);
      s.init(1,1,0, d_alp[1]*exp2*exp2, R);  ao.add_primitive( d_coeff[1]*coeff2, s);
      s.init(1,1,0, d_alp[2]*exp2*exp2, R);  ao.add_primitive( d_coeff[2]*coeff2, s);
      
      ao.ao_name = out+Atom_shell+"xy";          ao.is_ao_name = 1;      
//      ao.at_indx = Atom_index;              ao.is_at_indx = 1;       // 1, 2, ...
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(1,0,1, d_alp[0]*exp1*exp1, R);  ao.add_primitive( d_coeff[0]*coeff1, s);
      s.init(1,0,1, d_alp[1]*exp1*exp1, R);  ao.add_primitive( d_coeff[1]*coeff1, s);
      s.init(1,0,1, d_alp[2]*exp1*exp1, R);  ao.add_primitive( d_coeff[2]*coeff1, s);
      s.init(1,0,1, d_alp[0]*exp2*exp2, R);  ao.add_primitive( d_coeff[0]*coeff2, s);
      s.init(1,0,1, d_alp[1]*exp2*exp2, R);  ao.add_primitive( d_coeff[1]*coeff2, s);
      s.init(1,0,1, d_alp[2]*exp2*exp2, R);  ao.add_primitive( d_coeff[2]*coeff2, s);
      
      ao.ao_name = out+Atom_shell+"xz";          ao.is_ao_name = 1;      
//      ao.at_indx = Atom_index;              ao.is_at_indx = 1;       // 1, 2, ...
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


      s.init(0,1,1, d_alp[0]*exp1*exp1, R);  ao.add_primitive( d_coeff[0]*coeff1, s);
      s.init(0,1,1, d_alp[1]*exp1*exp1, R);  ao.add_primitive( d_coeff[1]*coeff1, s);
      s.init(0,1,1, d_alp[2]*exp1*exp1, R);  ao.add_primitive( d_coeff[2]*coeff1, s);
      s.init(0,1,1, d_alp[0]*exp2*exp2, R);  ao.add_primitive( d_coeff[0]*coeff2, s);
      s.init(0,1,1, d_alp[1]*exp2*exp2, R);  ao.add_primitive( d_coeff[1]*coeff2, s);
      s.init(0,1,1, d_alp[2]*exp2*exp2, R);  ao.add_primitive( d_coeff[2]*coeff2, s);
      
      ao.ao_name = out+Atom_shell+"yz";          ao.is_ao_name = 1;      
//      ao.at_indx = Atom_index;              ao.is_at_indx = 1;       // 1, 2, ...
      ao.element = Atom_name;               ao.is_element = 1;       // Li, Be, B, etc.
      ao.ao_shell = out+Atom_shell;         ao.is_ao_shell = 1;      // 3s, 3p, etc.
      ao.ao_shell_type = Atom_shell;        ao.is_ao_shell_type = 1; // s, p, d, etc.

      ao.normalize(); basis_ao.push_back(ao); ao.clear(); basis_size++; Norb++;


    }// Atom_shell=="d"

  }// Nzeta == 2

//  if(verbose_level>1){
//    cout<<"in add_ao: ao_name= "<<ao.ao_name<<" at_indx= "<<ao.at_indx<<" element= "<<ao.element<<endl;
//  }

}  // void add_basis_ao(...)


void add_basis_ao(std::string Atom_name, VECTOR& R, std::string Atom_shell, int Nzeta, int Nquant,
                    double  IP, double exp1, double exp2, double coeff1, double coeff2,
                    boost::python::list basis_ao){
/**
  \brief Create a new AO (STO-3G) object and add it to existing (including empty) list of orbitals (basis) - 
  this version produces the Python list of AO objects.

  \param[in] Atom_name The name of the atom (element) for which we want to create AO 
  \param[in] R The coordinate at which to create the atomic orbital. 
  \param[in] Atom_shell The name of the AO shell e.g. 1s, 2p, 3p, 3d, etc.
  \param[in] Nzeta The number of different Slater-type zetas - in most cases 1, for d-orbitals can be 2
  \param[in] Nquant The principal quantum number of the orbital
  \param[in] IP The state-specific ionization potential (not actually used in this function!)
  \param[in] exp1 The first Slater-type exponent
  \param[in] exp2 The second Slater-type exponent - is not used when Nzeta = 1
  \param[in] coeff1 The contraction coefficient of the first Slater zeta orbital
  \param[in] coeff2 The contraction coefficient of the second Slater zeta orbital (is effectively = 0 if Nzeta = 1)
  \param[in,out] basis_ao The list of AO objects - the basis we are creating

  Sequential application of this function with proper parameters will generate the atomic basis - the list of AO objects

  For conversion of the spherical hamonics to Cartesian coordinates see e.g.:
  http://csi.chemie.tu-darmstadt.de/ak/immel/script/redirect.cgi?filename=http://csi.chemie.tu-darmstadt.de/ak/immel/tutorials/orbitals/hydrogenic.html
  So the z2 orbital is actually: (1/(?)) *[ 2*z2 - (x2+x2) ]  
  and the x2-y2 orbital is (1/sqrt(2)) * [x2 - y2]

*/


  vector<AO> tmp;
  add_basis_ao(Atom_name,R,Atom_shell, Nzeta, Nquant, IP, exp1, exp2, coeff1, coeff2, tmp);

  for(int i=0;i<tmp.size();i++){
    basis_ao.append(tmp[i]);
  }

}



int num_valence_elec(int Z){
/**
  \brief The function to compute the number of valence electrons in the atom with given nucleus charge.

  param[in] Z The charge of the atomic nucleus
*/

  // Compute the number of valence electrons
  int dn = 0; 
       if(Z<=2)          {  dn = Z;          }  //  H,      ... He
  else if(Z>= 3 && Z<=10){  dn = Z - 2;      }  //  Li, Be, ... Ne
  else if(Z>=11 && Z<=18){  dn = Z - 10;     }  //  Na, Mg, ... Ar
  
  else if(Z>=19 && Z<=36){  dn = Z - 18;     }  //  K,  Ca, ... Kr  - here we include all d electrons as well
  else if(Z>=37 && Z<=54){  dn = Z - 36;     }  //  Rb, Sr, ... Xe  - here we include all d electrons as well

  else if(Z>=55 && Z<=56){  dn = Z - 54;     }  //  Cs, Ba
  else if(Z>=71 && Z<=86){  dn = Z - 70 + 2; }  //  Lu, Hf, ... Rn  - skip f-electrons

  else{
     cout<<"Calculation of # of electrons for Z>86 is not implemented\n"; exit(0);
  }

  return dn;

}




}// namespace libbasis
}// namespace liblibra


