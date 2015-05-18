/*********************************************************************************
* Copyright (C) 2014 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
 \file Dipole.cpp
 \brief Implementation of the transition dipole matrix computations in AO basis

*/
/****************************************************************************
  This file contains following functions:

  void update_dipole_matrix(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
                            vector<int>& basis_fo,vector<AO>& basis_ao,
                            MATRIX* Mux_ao,MATRIX* Muy_ao,MATRIX* Muz_ao,vector<double*>& aux, int n_aux, Nuclear& mol)

  void compute_dipole_moments(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                              vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                              Electronic* el, Memory* mem,
                              VECTOR& mu_nucl,MATRIX* mux,MATRIX* muy,MATRIX* muz)



****************************************************************************/
#include "Dipole.h"

void update_dipole_matrix(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
                          vector<int>& basis_fo,vector<AO>& basis_ao,
                          MATRIX* Mux_ao,MATRIX* Muy_ao,MATRIX* Muz_ao,vector<double*>& aux, int n_aux, Nuclear& mol){

  int i,j,n,I,J;
  VECTOR dIdA,dIdB,TV, Rij,Rij0;
  double dist, dist_min;  
  int opt_x,opt_y,opt_z;

  int Norb = basis_fo.size();
//  cout<<"In update_dipole_matrix\n"; 

  for(i=0;i<Norb;i++){
    I = basis_fo[i];
    for(j=0;j<Norb;j++){      
      J = basis_fo[j];

      TV = mol.get_tv(basis_ao[I].at_indx,basis_ao[J].at_indx,x_period,y_period,z_period,t1,t2,t3);

      VECTOR muij = DIPOLE_INTEGRAL(basis_ao[I],basis_ao[J],0,dIdA,dIdB,aux,n_aux,TV); // <i|r|j>

      // Electronic contribution
      Mux_ao->M[i*Norb+j] = -muij.x;
      Muy_ao->M[i*Norb+j] = -muij.y;
      Muz_ao->M[i*Norb+j] = -muij.z;


    }// j
  }// i

}


void compute_dipole_moments(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                            vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                            Electronic* el, Memory* mem,
                            VECTOR& mu_nucl,MATRIX* mux,MATRIX* muy,MATRIX* muz){

  int n;

  // Compute center of mass of a given fragment (subset of nuclei in molecular structure mol)
  VECTOR com; com = 0.0;
  double mass = 0.0;

  for(n=0;n<fragment.size();n++){
    com += mol.R[fragment[n]] * mol.mass[fragment[n]];
    mass += mol.mass[fragment[n]];
  }
  com = com/mass;


  cout<<"Total mass of the fragment is: "<<mass<<endl;
  cout<<"Coordinates of the center of mass of the fragment are: "<<com<<endl;


  // Compute nuclear dipole moment in MO
  // Nuclear contribution to dipole moment
  mu_nucl = 0.0;
  for(n=0;n<fragment.size();n++){
    mu_nucl += mol.Zeff[fragment[n]] * (mol.R[fragment[n]] - com);
  }

  cout<<"Nuclear contribution to dipole momentum (in MO basis), a.u. = "<<mu_nucl<<endl;


  // Precompute transition dipole moments in AO basis
  // Bare dipole momenta <xi_i|r|xi_j> - not translationally invariant
  update_dipole_matrix(prms.x_period,prms.y_period,prms.z_period,prms.t1,prms.t2,prms.t3,
                         basis_fo, basis_ao, mux, muy, muz, mem->aux, mem->n_aux, mol);

  // To make dipole momentum translationally-invariant, subtract component due to COM
  // working in AO,so need overlam Sao
  // mu = -<i|r-COM|j> 
  //                   
  *mux = *mux + com.x * *el->Sao;// + mu_nucl.x * *el[i]->Sao;
  *muy = *muy + com.y * *el->Sao;// + mu_nucl.y * *el[i]->Sao;
  *muz = *muz + com.z * *el->Sao;// + mu_nucl.z * *el[i]->Sao;

  // Note, total dipole moment (in MO) will include the contrimution from nuclei
  //+ <i| sum[ Z_a * (R_a - COM)] |j>  - add nuclear contribution at the very end !!! 
  //       a                                                                          


}


