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
/**
  \file Hamiltonian_INDO.cpp
  \brief The file implements functions for INDO calculations
*/

#include "Hamiltonian_INDO.h"


/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{


vector<int> compute_sorb_indices
( int sz, vector<AO>& basis_ao, vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
){
/** 
  Compute the global indices of the last s-type orbitals on each atom

  \param[in] sz The number of atoms in the system
  \param[in] basis_ao The AO basis for the system, including s-type and other functions
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized

  The function returns the vector of the indices of the first s-type AOs in the global array of AOs, sorb_indx, so
  that is sorb_indx[i] is the index of the last found s-type orbital localized on the atom with index i.
  We assume, that when the basis functions are created, the order of generation of orbitals is
  1s, 2s, ....,  3s, ...., etc so the returned indices refer to the valence shell s-type orbital
    
*/


  vector<int> sorb_indx(sz,0); // global index of s-type orbital on i-th atom

  for(int a=0;a<sz;a++){  // for all atoms
    for(int i=0;i<atom_to_ao_map[a].size();i++){  // all orbitals on given atom (i-dummy)

      int I = atom_to_ao_map[a][i];  // i-th AO on atom a, I - is the global index of this AO in the given basis

      if(basis_ao[I].ao_shell_type=="s" ){  sorb_indx[a] = I; }
          
    }// for j
  }// for i

  return sorb_indx;
}


void compute_indo_core_parameters
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, double& eri, double& V_AB){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] sorb_indx The vector of global indices of the last s-type orbitals on each atom
  \param[in] opt Option for computing V_AB terms: this controlls the distinction between INDO (opt = 1) and CNDO2 (opt = 0)
  \param[in] a the index of one of the atoms, for which the V_ab term is computed
  \param[in] b the index of one of the atoms, for which the V_ab term is computed
  \param[out] eri The electron repulsion integral between a and b cores
  \param[out] V_AB The core repulsion integral
  
  Computes ERIs and V_AB parameters for given pair of atoms, and for given geometry.
*/

  // Compute ERIs and V_AB

  int I = sorb_indx[a];
  int J = sorb_indx[b];
     
  // ERI
  eri = electron_repulsion_integral(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]); // eri[a][b]

  // V_AB
  int B = b;
  double Zeff = modprms.PT[syst.Atoms[B].Atom_element].Zeff; 

  if(opt==0){
    V_AB = Zeff*nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm);// V_AB[a][b]
  }
  else if(opt==1){
    V_AB = Zeff*eri; //[a*sz+b];
  }

}

void compute_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, int c, VECTOR& deri, VECTOR& dV_AB){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] sorb_indx The vector of global indices of the last s-type orbitals on each atom
  \param[in] opt Option for computing V_AB terms: this controlls the distinction between INDO (opt = 1) and CNDO2 (opt = 0)
  \param[in] a the index of one of the atoms, for which the V_ab term is computed
  \param[in] b the index of one of the atoms, for which the V_ab term is computed
  \param[in] c the index pf the atom w.r.t. which the derivatives are computed
  \param[out] deri The derivative of the electron repulsion integral between a and b cores
  \param[out] dV_AB The derivative of the core repulsion integral
  
  Compute derivatives of ERI and V_AB parameters:  d ERI[a][b] / dR[c]  and d V[a][b] / dR[c]
*/




  int I = sorb_indx[a];
  int J = sorb_indx[b];
  int K = sorb_indx[c]; 
     
  // ERI
  //eri = electron_repulsion_integral(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]); // eri[a][b]


  VECTOR DA,DB,DC,DD;
  double eri = electron_repulsion_integral(&basis_ao[I],&basis_ao[I],&basis_ao[J],&basis_ao[J],1,1,DA,DB,DC,DD);

  VECTOR deri_dc; deri_dc = 0.0;
  if(a==b){ ;; }
  else{
    if(c==a){ deri += (DA + DB); }
    if(c==b){ deri += (DC + DD); }
  }


  // V_AB
  double V_AB = 0.0;
  int B = b; //ao_to_atom_map[b]; // global index of atom on orbital b
  double Zeff = modprms.PT[syst.Atoms[B].Atom_element].Zeff; 

  if(opt==0){
    V_AB = Zeff*nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm);// V_AB[a][b]

    double nai = nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm, 1, 1, DA, DB, DC);// V_AB[a][b]

    // Alright - i'm not really sure about this, so skip this part for now
    //if(K==J){ dV_AB += (DA + DB); }
    //if(K==J){ dV_AB += (DA + DB); }

  }
  else if(opt==1){
    V_AB = Zeff*eri; //[a*sz+b];
    dV_AB = Zeff * deri;
  }

}



void compute_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, int c, VECTOR& deri, VECTOR& dV_AB,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] sorb_indx The vector of global indices of the last s-type orbitals on each atom
  \param[in] opt Option for computing V_AB terms: this controlls the distinction between INDO (opt = 1) and CNDO2 (opt = 0)
  \param[in] a the index of one of the atoms, for which the V_ab term is computed
  \param[in] b the index of one of the atoms, for which the V_ab term is computed
  \param[in] c the index pf the atom w.r.t. which the derivatives are computed
  \param[out] deri The derivative of the electron repulsion integral between a and b cores
  \param[out] dV_AB The derivative of the core repulsion integral
  \param[in,out] aux The auxiliary memory allocated for double values
  \param[in] n_aux The length of each of the allocated double array
  \param[in,out] auxv The auxiliary memory allocated for VECTOR values
  \param[in] n_auxv The length of each of the allocated VECTOR array 
  
  Compute derivatives of ERI and V_AB parameters:  d ERI[a][b] / dR[c]  and d V[a][b] / dR[c] for given pair of atoms
  and for the selected gradient component
  This is supposed to be an accelerated version, since no memory allocation/deallocation is necessary
*/


  int I = sorb_indx[a];
  int J = sorb_indx[b];
  int K = sorb_indx[c]; 
     
  // ERI
  //eri = electron_repulsion_integral(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]); // eri[a][b]


  VECTOR DA,DB,DC,DD;

  /// This version doesn't do memory re-allocation every time

  double eri = electron_repulsion_integral(&basis_ao[I],&basis_ao[I],&basis_ao[J],&basis_ao[J],1,1,DA,DB,DC,DD, aux, n_aux, auxv, n_auxv);


  VECTOR deri_dc; deri_dc = 0.0;
  if(a==b){ ;; }
  else{
    if(c==a){ deri += (DA + DB); }
    if(c==b){ deri += (DC + DD); }
  }


  // V_AB
  double V_AB = 0.0;
  int B = b; //ao_to_atom_map[b]; // global index of atom on orbital b
  double Zeff = modprms.PT[syst.Atoms[B].Atom_element].Zeff; 

  if(opt==0){
    V_AB = Zeff*nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm);// V_AB[a][b]

    double nai = nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm, 1, 1, DA, DB, DC);// V_AB[a][b]
    // Alright - i'm not really sure about this, so skip this part for now
    //if(K==J){ dV_AB += (DA + DB); }
    //if(K==J){ dV_AB += (DA + DB); }

  }
  else if(opt==1){
    V_AB = Zeff*eri; //[a*sz+b];
    dV_AB = Zeff * deri;
  }

}

void compute_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  vector<int>& sorb_indx,
  int opt, int a, int b, vector<VECTOR>& deri, vector<VECTOR>& dV_AB,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] sorb_indx The vector of global indices of the last s-type orbitals on each atom
  \param[in] opt Option for computing V_AB terms: this controlls the distinction between INDO (opt = 1) and CNDO2 (opt = 0)
  \param[in] a the index of one of the atoms, for which the V_ab term is computed
  \param[in] b the index of one of the atoms, for which the V_ab term is computed
  \param[out] deri Thevector of derivatives of the electron repulsion integral between a and b cores w.r.t to each nuclear DOF
  \param[out] dV_AB The derivative of the core repulsion integral between a and b cores w.r.t to each nuclear DOF
  \param[in,out] aux The auxiliary memory allocated for double values
  \param[in] n_aux The length of each of the allocated double array
  \param[in,out] auxv The auxiliary memory allocated for VECTOR values
  \param[in] n_auxv The length of each of the allocated VECTOR array 
  
  Compute derivatives of ERI and V_AB parameters:  d ERI[a][b] / dR[c]  and d V[a][b] / dR[c] for given pair of atoms
  The derivatives w.r.t. all atoms are computed at once, in a single swipe - this is much more efficient approach than 
  when we call this computations one by one.
  This is supposed to be an accelerated version, since no memory allocation/deallocation is necessary
*/



  // compute derivatives only once - for a fixed pair of a and b
  VECTOR DA,DB,DC,DD;  

  int I = sorb_indx[a];
  int J = sorb_indx[b];

  double eri = electron_repulsion_integral(&basis_ao[I],&basis_ao[I],&basis_ao[J],&basis_ao[J],1,1,DA,DB,DC,DD, aux, n_aux, auxv, n_auxv);



  if(deri.size()!=syst.Number_of_atoms){  deri = vector<VECTOR>(syst.Number_of_atoms, VECTOR(0.0, 0.0, 0.0)); }
  if(dV_AB.size()!=syst.Number_of_atoms){  dV_AB = vector<VECTOR>(syst.Number_of_atoms, VECTOR(0.0, 0.0, 0.0)); }


  for(int c=0;c<syst.Number_of_atoms;c++){  // all atoms

    // Now set up derivatives
    deri[c] = 0.0;
    dV_AB[c] = 0.0;

    if(c==a){ deri[c] += (DA + DB); }
    if(c==b){ deri[c] += (DC + DD); }

    double Zeff = modprms.PT[syst.Atoms[b].Atom_element].Zeff; 

    if(opt==0){
    }
    else if(opt==1){
      dV_AB[c] = Zeff * deri[c];
    }


  }// for c


}




void compute_all_indo_core_parameters_derivs
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, 
  vector<int>& sorb_indx, int opt
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] sorb_indx The vector of global indices of the last s-type orbitals on each atom
  \param[in] opt Option for computing V_AB terms: this controlls the distinction between INDO (opt = 1) and CNDO2 (opt = 0)
  
  Compute the ERI and V_AB parameters for all pairs of atoms, a and b.
  Also, compute all derivatives of each ERI and V_AB parameter:  d ERI[a][b] / dR[c]  and d V[a][b] / dR[c]  w.r.t. all atoms 

  The ERI and V_AB matrices are stored internally, while the derivatives are printed out in file in binary format
  So, the following files will be created in the calling directory:
  "deri.x_1.bin", "deri.x_2.bin", ... , "deri.x_N.bin"  (where N - is the number of atoms)
  "deri.y_1.bin", "deri.y_2.bin", ... , "deri.y_N.bin"  (where N - is the number of atoms)
  "deri.z_1.bin", "deri.z_2.bin", ... , "deri.z_N.bin"  (where N - is the number of atoms)
  "dV_AB.x_1.bin", "dV_AB.x_2.bin", ... , "dV_AB.x_N.bin"  (where N - is the number of atoms)
  "dV_AB.y_1.bin", "dV_AB.y_2.bin", ... , "dV_AB.y_N.bin"  (where N - is the number of atoms)
  "dV_AB.z_1.bin", "dV_AB.z_2.bin", ... , "dV_AB.z_N.bin"  (where N - is the number of atoms)

  These files will be accessed by the SCF procedures, so they are needed

  This is the older (MUCH LESS EFFICIENT) version which takes O(N^3) computations because it computes derivatives one by one for
  each nuclear DOF.
*/


  int N = syst.Number_of_atoms;


  MATRIX* aux1;  aux1 = new MATRIX(N,N);
  MATRIX* aux2;  aux2 = new MATRIX(N,N);
  MATRIX* aux3;  aux3 = new MATRIX(N,N);

  MATRIX* aux4;  aux4 = new MATRIX(N,N);
  MATRIX* aux5;  aux5 = new MATRIX(N,N);
  MATRIX* aux6;  aux6 = new MATRIX(N,N);

  // Memory for ERI computations
  int i;
  int n_auxd = 40;
  int n_auxv = 40;
  vector<double*> auxd(30);
  for(i=0;i<30;i++){ auxd[i] = new double[n_auxd]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }


  for(int c=0;c<N;c++){

    for(int a=0;a<N;a++){
      for(int b=0;b<N;b++){

        VECTOR deri, dV_AB;
        deri = 0.0; dV_AB = 0.0;

//        compute_indo_core_parameters_derivs(syst,basis_ao,modprms, atom_to_ao_map, ao_to_atom_map, sorb_indx, opt, a, b, c, deri, dV_AB);

        // Memory-efficient version
        compute_indo_core_parameters_derivs(syst,basis_ao,modprms, atom_to_ao_map, ao_to_atom_map, sorb_indx, opt, a, b, c, deri, dV_AB, auxd, n_auxd, auxv, n_auxv);
    
        aux1->set(a, b, deri.x);
        aux2->set(a, b, deri.y);
        aux3->set(a, b, deri.z);

        aux4->set(a, b, dV_AB.x);
        aux5->set(a, b, dV_AB.y);
        aux6->set(a, b, dV_AB.z);

      }// for b
    }// for a

    stringstream ss(stringstream::in | stringstream::out);
    std::string out;
    (ss << c);  ss >> out;

    aux1->bin_dump("deri.x_"+out+".bin");
    aux2->bin_dump("deri.y_"+out+".bin");
    aux3->bin_dump("deri.z_"+out+".bin");

    aux4->bin_dump("dV_AB.x_"+out+".bin");
    aux5->bin_dump("dV_AB.y_"+out+".bin");
    aux6->bin_dump("dV_AB.z_"+out+".bin");

  }// for c

  // Clean working memory
  for(i=0;i<30;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();
 


  delete aux1; delete aux2; delete aux3;
  delete aux4; delete aux5; delete aux6;
  
}



void compute_all_indo_core_parameters_derivs1
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, 
  vector<int>& sorb_indx, int opt
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] sorb_indx The vector of global indices of the last s-type orbitals on each atom
  \param[in] opt Option for computing V_AB terms: this controlls the distinction between INDO (opt = 1) and CNDO2 (opt = 0)
  
  Compute the ERI and V_AB parameters for all pairs of atoms, a and b.
  Also, compute all derivatives of each ERI and V_AB parameter:  d ERI[a][b] / dR[c]  and d V[a][b] / dR[c]  w.r.t. all atoms 

  The ERI and V_AB matrices are stored internally, while the derivatives are printed out in file in binary format
  So, the following files will be created in the calling directory:
  "deri.x_1.bin", "deri.x_2.bin", ... , "deri.x_N.bin"  (where N - is the number of atoms)
  "deri.y_1.bin", "deri.y_2.bin", ... , "deri.y_N.bin"  (where N - is the number of atoms)
  "deri.z_1.bin", "deri.z_2.bin", ... , "deri.z_N.bin"  (where N - is the number of atoms)
  "dV_AB.x_1.bin", "dV_AB.x_2.bin", ... , "dV_AB.x_N.bin"  (where N - is the number of atoms)
  "dV_AB.y_1.bin", "dV_AB.y_2.bin", ... , "dV_AB.y_N.bin"  (where N - is the number of atoms)
  "dV_AB.z_1.bin", "dV_AB.z_2.bin", ... , "dV_AB.z_N.bin"  (where N - is the number of atoms)

  These files will be accessed by the SCF procedures, so they are needed

  This is the new (MUCH MORE EFFICIENT) version which takes O(N^2) computations because it computes derivatives in batches.
*/


  int N = syst.Number_of_atoms;


  MATRIX* aux1;  aux1 = new MATRIX(N,N);
  MATRIX* aux2;  aux2 = new MATRIX(N,N);
  MATRIX* aux3;  aux3 = new MATRIX(N,N);

  MATRIX* aux4;  aux4 = new MATRIX(N,N);
  MATRIX* aux5;  aux5 = new MATRIX(N,N);
  MATRIX* aux6;  aux6 = new MATRIX(N,N);

  // Memory for ERI computations
  int i;
  int n_auxd = 40;
  int n_auxv = 40;
  vector<double*> auxd(30);
  for(i=0;i<30;i++){ auxd[i] = new double[n_auxd]; }
  vector<VECTOR*> auxv(5);
  for(i=0;i<5;i++){ auxv[i] = new VECTOR[n_auxv]; }


  // This is bad memory scaling, i know, but should be more-or-less fine for now
  vector< vector<vector<VECTOR> > > deri;
  vector< vector<vector<VECTOR> > > dV_AB;
  deri = vector< vector<vector<VECTOR> > >(N, vector<vector<VECTOR> >(N, vector<VECTOR>(N, VECTOR(0.0, 0.0, 0.0) ) ) );
  dV_AB = vector< vector<vector<VECTOR> > >(N, vector<vector<VECTOR> >(N, vector<VECTOR>(N, VECTOR(0.0, 0.0, 0.0) ) ) );


  for(int a=0;a<N;a++){
    for(int b=0;b<N;b++){

      // Memory-efficient version
      compute_indo_core_parameters_derivs(syst,basis_ao,modprms, atom_to_ao_map, ao_to_atom_map, sorb_indx, opt, a, b , deri[a][b], dV_AB[a][b], auxd, n_auxd, auxv, n_auxv);
    
    }// for b
  }// for a



  for(int c=0;c<N;c++){

    // Form matrices
    for(int a=0;a<N;a++){
      for(int b=0;b<N;b++){

        aux1->set(a, b, deri[a][b][c].x);
        aux2->set(a, b, deri[a][b][c].y);
        aux3->set(a, b, deri[a][b][c].z);

        aux4->set(a, b, dV_AB[a][b][c].x);
        aux5->set(a, b, dV_AB[a][b][c].y);
        aux6->set(a, b, dV_AB[a][b][c].z);

      }// for b
    }// for a

    // Print matrices
    stringstream ss(stringstream::in | stringstream::out);
    std::string out;
    (ss << c);  ss >> out;

    aux1->bin_dump("deri.x_"+out+".bin");
    aux2->bin_dump("deri.y_"+out+".bin");
    aux3->bin_dump("deri.z_"+out+".bin");

    aux4->bin_dump("dV_AB.x_"+out+".bin");
    aux5->bin_dump("dV_AB.y_"+out+".bin");
    aux6->bin_dump("dV_AB.z_"+out+".bin");

  }// for c

  // Clean working memory
  for(i=0;i<30;i++){ delete [] auxd[i]; }  
  auxd.clear();
  for(i=0;i<5;i++){ delete [] auxv[i]; }  
  auxv.clear();
 


  delete aux1; delete aux2; delete aux3;
  delete aux4; delete aux5; delete aux6;
  
}





void indo_core_parameters
( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int opt, int DF){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] opt Option for computing V_AB terms: this controlls the distinction between INDO (opt = 1) and CNDO2 (opt = 0)
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  
  The upper-level function for initializing INDO parameters and their derivatives
*/


// opt == 0 - cndo
// opt == 1 - indo
  modprms.indo_opt = opt;

//  int DF = 0;
  int i,j,a,b,i1,a1,I,J,A;
  VECTOR da,db,dc;

  if(DF){ cout<<"in indo_core_parameters\n"; }

  /// Allocates memory, if needed

  int sz = syst.Number_of_atoms; // number of atoms in the system

  if(modprms.eri.size()!=sz*sz){  cout<<"In indo_core_parameters: eri array is not allocated\nDo allocation...\n"; 
    modprms.eri.clear(); modprms.eri = vector<double>(sz*sz,0.0); 
  }
  if(modprms.V_AB.size()!=sz*sz){  cout<<"In indo_core_parameters: V_AB array is not allocated\nDo allocation...\n"; 
    modprms.V_AB.clear(); modprms.V_AB = vector<double>(sz*sz,0.0); 
  }


  /// Compute the indices of the valence shell s-type orbital

  vector<int> sorb_indx;
  sorb_indx = compute_sorb_indices(sz,basis_ao,atom_to_ao_map,ao_to_atom_map);


  // Printing mapping
  if(DF){ 
    cout<<"i - runs over indices of atoms in given system\n";
    cout<<"sorb_indx[i] - is the global index of s-type orbital (assuming only one) centered on atom i\n";
    for(i=0;i<sorb_indx.size();i++){
      cout<<"i= "<<i<<" sorb_indx[i]= "<<sorb_indx[i]<<endl;
    }
  }
    
/*
  // Compute ERIs and V_AB
  for(a=0;a<sz;a++){
    for(b=0;b<sz;b++){


      I = sorb_indx[a];
      J = sorb_indx[b];
     
      // ERI
      modprms.eri[a*sz+b] = electron_repulsion_integral(basis_ao[I],basis_ao[I],basis_ao[J],basis_ao[J]);

      // V_AB
      int B = b; //ao_to_atom_map[b]; // global index of atom on orbital b
      if(opt==0){
        modprms.V_AB[a*sz+b] = modprms.PT[syst.Atoms[B].Atom_element].Zeff*nuclear_attraction_integral(basis_ao[J],basis_ao[J], syst.Atoms[B].Atom_RB.rb_cm);// V_AB[a][b]
      }
      else if(opt==1){
        modprms.V_AB[a*sz+b] = modprms.PT[syst.Atoms[B].Atom_element].Zeff*modprms.eri[a*sz+b];
      }

      if(DF){ cout<<"a= "<<a<<" b= "<<b<<" I= "<<I<<" J= "<<J<<" eri= "<<modprms.eri[a*sz+b]<<" V_AB= "<<modprms.V_AB[a*sz+b]<<endl; }

    }// for j
  }// for i
*/

  /// Compute ERIs and V_AB

  for(a=0;a<sz;a++){
    for(b=0;b<sz;b++){

    compute_indo_core_parameters
    (syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, sorb_indx, opt, a, b, modprms.eri[a*sz+b], modprms.V_AB[a*sz+b] ); 

    }// for b
  }// for a

  /// Compute their derivatives and store on the disk

  //compute_all_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, sorb_indx, opt);
  compute_all_indo_core_parameters_derivs1(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, sorb_indx, opt);


}




void Hamiltonian_core_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The pointer to the matrix object in which the core Hamiltonian will be stored
  \param Sao The pointer to the AO overla matrix (not actually used here)
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  
  Compute the core INDO Hamiltonian
*/


  int i,j,k,a,b,I,J,A,B;
  VECTOR dIdA,dIdB;
  //cout<<"in Hamiltonian_core_indo\n";

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_indo is called\n";
    cout<<"In Hamiltonian_core_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  int sz = syst.Number_of_atoms; // number of atoms in this fragment

  if(modprms.eri.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_indo: size of auxiliary eri array is not right\n"; exit(0);}
  if(modprms.V_AB.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_indo: size of auxiliary V_AB array is not right\n"; exit(0);}



  //----------- Compute core terms of the Hamiltonian ------------
  *Hao = 0.0;

  for(i=0;i<Norb;i++){  // global orbital indices
    // values of IP are different for cndo (just IPs) and cndo2 and indo ( 0.5*(IP + EA) )
    a = ao_to_atom_map[i];
    Hao->M[i*Norb+i] += modprms.PT[basis_ao[i].element].IP[basis_ao[i].ao_shell];

    if(DF){
      cout<<"Setting diagonal element i = "<<i<<endl;
      cout<<"Contribution from IP = "<<Hao->M[i*Norb+i]<<endl;
    }

    /// The code below is the same for CNDO, CNDO2 and INDO - but the difference comes in use of different G1 and F2 parameters
    /// for CNDO and CNDO2 they are zero
    double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
    double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
    
    if(DF){
      cout<<"a= "<<a<<" G1= "<<G1<<" F2= "<<F2<<" Atom[a].Atom_Z= "<<syst.Atoms[a].Atom_Z
          <<" basis_ao[i].ao_shell_type= "<<basis_ao[i].ao_shell_type<<endl;
    }
   
    /// Eqs. 3.17 - 3.23 from Pople, Beveridge, Dobosh, JCP 47, 2026 (1967)
    ///
    int Z = syst.Atoms[a].Atom_Z;  // modprms.PT[elt].Z - e.g. 6 for C

         if(Z==1){ Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a]; }  // H
    else if(Z==3 || Z==11){ 

      if(basis_ao[i].ao_shell_type=="s"){   Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a]; } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (0.5*modprms.eri[a*sz+a] - G1/12.0); } // p

    }  // Li or Na

    else if(Z==4 || Z==12){ 

      if(basis_ao[i].ao_shell_type=="s"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.5*G1); } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.25*G1); } // p

    }  // Be or Mg

    else if( (Z>=5 && Z<=9) || (Z>=13 && Z<=18)){ // B - F or Al - Cl

      double Z_core = modprms.PT[basis_ao[i].element].Nval; // core charge of atom 

      if(basis_ao[i].ao_shell_type=="s"){      
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/6.0)*( Z_core - 1.5 )*G1); // s
      }
      else if(basis_ao[i].ao_shell_type=="p"){
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/3.0)*G1 - 0.08*( Z_core - 2.5 )*F2); // p
      }

    }
    else{  cout<<"Error: INDO is not implemented for elements beyond Cl\n"; exit(0);  }
    if(DF){ cout<<" + Contribution from Frank-Condon factors = "<<Hao->M[i*Norb+i]<<endl;   }

    
    //----------------- Coulombic terms --------------

    for(b=0;b<sz;b++){
      if(b!=a){
        Hao->M[i*Norb+i] -= modprms.V_AB[a*sz+b];  //  = V_AB = Z_B * eri[A][B] -in INDO
      }
    }// for b
    if(DF){  cout<<" + Contribution from Coulombic terms = "<<Hao->M[i*Norb+i]<<endl;   }


    //-------------- Off-diagonal terms of the core matrix ---------
    for(j=0;j<Norb;j++){

      if(j!=i){
        b = ao_to_atom_map[j];

        if(b==a){ ;; }  // different orbitals centered on the same atom - give zero (not true for hybrid orbitals)
        else{           // centered on different atoms - use overlap formula

          // Overlap is set to identity in INDO, so need to recompute it explicitly
          double sao_ij = gaussian_overlap(basis_ao[i],basis_ao[j]); // 0, dIdA,dIdB), mem->aux, mem->n_aux);

          Hao->M[i*Norb+j] += 0.5*(modprms.PT[basis_ao[i].element].beta0[basis_ao[i].ao_shell] 
                                 + modprms.PT[basis_ao[j].element].beta0[basis_ao[j].ao_shell]) * sao_ij;
        }
      }// j!=i

    }// for j    

  }// for i

}


void Hamiltonian_core_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The matrix object in which the core Hamiltonian will be stored
  \param Sao The AO overla matrix (not actually used here)
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  
  Compute the core INDO Hamiltonian - Python-friendly version
*/


  Hamiltonian_core_indo( syst, basis_ao, prms, modprms,  atom_to_ao_map, ao_to_atom_map, &Hao, &Sao, DF);

}


void Hamiltonian_core_deriv_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX* Hao, MATRIX* Sao, int DF,
  int c,
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz, 
  MATRIX* dSao_dx, MATRIX* dSao_dy, MATRIX* dSao_dz
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The pointer to the matrix object in which the core Hamiltonian will be stored
  \param[out] Sao The pointer to the AO overlap matrix computed here
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  \param[in] c The index of the atom w.r.t. which coordinates we take the derivatives
  \param[out] dHao_dx The derivative of the Hamiltonian w.r.t. the x-coordinate of the selected atom
  \param[out] dHao_dy The derivative of the Hamiltonian w.r.t. the y-coordinate of the selected atom
  \param[out] dHao_dz The derivative of the Hamiltonian w.r.t. the z-coordinate of the selected atom
  \param[out] dSao_dx The derivative of the AO overlap matrix w.r.t. the x-coordinate of the selected atom
  \param[out] dSao_dy The derivative of the AO overlap matrix w.r.t. the y-coordinate of the selected atom
  \param[out] dSao_dz The derivative of the AO overlap matrix w.r.t. the z-coordinate of the selected atom
  
  Compute the core INDO Hamiltonian and its derivatives w.r.t. specified nuclear DOFs.
*/


  //================ Basically, here we compute derivatives of the core Hamiltonian ========================

  int i,j,k,a,b,I,J,A,B;
  VECTOR dIdA,dIdB;
//  cout<<"in Hamiltonian_core_deriv_indo\n";

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=Hao->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dx->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dy->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }
  if(Norb!=dHao_dz->n_cols){  
    cout<<"Hao matrix is not allocated\n Must be allocated before Hamiltonian_core_deriv_indo is called\n";
    cout<<"In Hamiltonian_core_deriv_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }

  int sz = syst.Number_of_atoms; // number of atoms in this fragment


  if(modprms.eri.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_deriv_indo: size of auxiliary eri array is not right\n"; exit(0);}
  if(modprms.V_AB.size()!=sz*sz){  cout<<"Error in Hamiltonian_core_deriv_indo: size of auxiliary V_AB array is not right\n"; exit(0);}


  stringstream ss(stringstream::in | stringstream::out);
  std::string out;
  (ss << c);  ss >> out;


  int use_disk = 1;

  

  MATRIX* aux1;  
  MATRIX* aux2;  
  MATRIX* aux3;  
  MATRIX* aux4;  
  MATRIX* aux5;  
  MATRIX* aux6;  

  if(use_disk){

    aux1 = new MATRIX(sz,sz);
    aux2 = new MATRIX(sz,sz);
    aux3 = new MATRIX(sz,sz);
    aux4 = new MATRIX(sz,sz);
    aux5 = new MATRIX(sz,sz);
    aux6 = new MATRIX(sz,sz);

    // Load matrices with derivatives of the core parameters
    aux1->bin_load("deri.x_"+out+".bin");
    aux2->bin_load("deri.y_"+out+".bin");
    aux3->bin_load("deri.z_"+out+".bin");

    aux4->bin_load("dV_AB.x_"+out+".bin");
    aux5->bin_load("dV_AB.y_"+out+".bin");
    aux6->bin_load("dV_AB.z_"+out+".bin");

  }

  //----------- Compute core terms of the Hamiltonian ------------
  *Hao = 0.0;
  *dHao_dx = 0.0;
  *dHao_dy = 0.0;
  *dHao_dz = 0.0;

  *dSao_dx = 0.0;
  *dSao_dy = 0.0;
  *dSao_dz = 0.0;
  


  vector<int> sorb_indx;
  sorb_indx = compute_sorb_indices(sz, basis_ao, atom_to_ao_map, ao_to_atom_map);


  for(i=0;i<Norb;i++){  // global orbital indices
    // values of IP are different for cndo (just IPs) and cndo2 and indo ( 0.5*(IP + EA) )
    a = ao_to_atom_map[i];    
    Hao->M[i*Norb+i] += modprms.PT[basis_ao[i].element].IP[basis_ao[i].ao_shell];

    if(DF){
      cout<<"Setting diagonal element i = "<<i<<endl;
      cout<<"Contribution from IP = "<<Hao->M[i*Norb+i]<<endl;
    }

    /// The code below is the same for CNDO, CNDO2 and INDO - but the difference comes in use of different G1 and F2 parameters
    /// for CNDO and CNDO2 they are zero

    double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
    double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
    
    if(DF){
      cout<<"a= "<<a<<" G1= "<<G1<<" F2= "<<F2<<" Atom[a].Atom_Z= "<<syst.Atoms[a].Atom_Z
          <<" basis_ao[i].ao_shell_type= "<<basis_ao[i].ao_shell_type<<endl;
    }
   
    /// Eqs. 3.17 - 3.23 from Pople, Beveridge, Dobosh, JCP 47, 2026 (1967)
    ///

    int Z = syst.Atoms[a].Atom_Z;  // modprms.PT[elt].Z - e.g. 6 for C


    //================== The portion below does not contribute to the derivatives ====================

         if(Z==1){      Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a];    }  // H
    else if(Z==3 || Z==11){ 

      if(basis_ao[i].ao_shell_type=="s"){   Hao->M[i*Norb+i] -= 0.5*modprms.eri[a*sz+a]; } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (0.5*modprms.eri[a*sz+a] - G1/12.0); } // p

    }  // Li or Na

    else if(Z==4 || Z==12){ 

      if(basis_ao[i].ao_shell_type=="s"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.5*G1); } // s
      else if(basis_ao[i].ao_shell_type=="p"){  Hao->M[i*Norb+i] -= (1.5*modprms.eri[a*sz+a] - 0.25*G1); } // p

    }  // Be or Mg

    else if( (Z>=5 && Z<=9) || (Z>=13 && Z<=18)){ // B - F or Al - Cl

      double Z_core = modprms.PT[basis_ao[i].element].Nval; // core charge of atom 

      if(basis_ao[i].ao_shell_type=="s"){      
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/6.0)*( Z_core - 1.5 )*G1); // s
      }
      else if(basis_ao[i].ao_shell_type=="p"){
        Hao->M[i*Norb+i] -= ( ( Z_core - 0.5 )*modprms.eri[a*sz+a] - (1.0/3.0)*G1 - 0.08*( Z_core - 2.5 )*F2); // p
      }

    }
    else{  cout<<"Error: INDO is not implemented for elements beyond Cl\n"; exit(0);  }
    if(DF){ cout<<" + Contribution from Frank-Condon factors = "<<Hao->M[i*Norb+i]<<endl;   }

    
    //----------------- Coulombic terms: Contribution to diagonal elements --------------


    for(b=0;b<sz;b++){
//      cout<<"i= "<<i<<" a= "<<a<<" b= "<<b<<endl;

      if(b!=a){
        Hao->M[i*Norb+i] -= modprms.V_AB[a*sz+b];  //  = V_AB = Z_B * eri[A][B] -in INDO


        VECTOR deri, dV_AB;    
        if(use_disk){

          dV_AB.x = aux4->get(a,b);
          dV_AB.y = aux5->get(a,b);
          dV_AB.z = aux6->get(a,b);

        }
        else{
          compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                              sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);
        }

        dHao_dx->M[i*Norb+i] -= dV_AB.x;
        dHao_dy->M[i*Norb+i] -= dV_AB.y;
        dHao_dz->M[i*Norb+i] -= dV_AB.z;
        


      }
    }// for b
    if(DF){  cout<<" + Contribution from Coulombic terms = "<<Hao->M[i*Norb+i]<<endl;   }


    //-------------- Off-diagonal terms of the core matrix ---------
    for(j=0;j<Norb;j++){

      if(j!=i){
        b = ao_to_atom_map[j];

        if(b==a){ ;; }  // different orbitals centered on the same atom - give zero (not true for hybrid orbitals)
        else{           // centered on different atoms - use overlap formula

          /// Overlap is set to identity in INDO, so need to recompute it explicitly

          //double sao_ij = gaussian_overlap(basis_ao[i],basis_ao[j]); // 0, dIdA,dIdB), mem->aux, mem->n_aux);

          VECTOR dSda, dSdb, dSdc;
          double sao_ij = gaussian_overlap(&basis_ao[i],&basis_ao[j], 1, 1, dSda, dSdb);
          // i - on atom a
          // j - on atom b

          dSdc = 0.0;
          if(c==a){ dSdc += dSda; }
          if(c==b){ dSdc += dSdb; }

          double beta_ij = 0.5*(modprms.PT[basis_ao[i].element].beta0[basis_ao[i].ao_shell] + modprms.PT[basis_ao[j].element].beta0[basis_ao[j].ao_shell]);

          Hao->M[i*Norb+j] += beta_ij * sao_ij;

          dHao_dx->M[i*Norb+j] += beta_ij * dSdc.x;
          dHao_dy->M[i*Norb+j] += beta_ij * dSdc.y;
          dHao_dz->M[i*Norb+j] += beta_ij * dSdc.z;

          dSao_dx->M[i*Norb+j] += dSdc.x;
          dSao_dy->M[i*Norb+j] += dSdc.y;
          dSao_dz->M[i*Norb+j] += dSdc.z;


        }
      }// j!=i

    }// for j    

  }// for i

  if(use_disk){
    delete aux1; delete aux2; delete aux3;
    delete aux4; delete aux5; delete aux6;
  }

}

void Hamiltonian_core_deriv_indo
( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF,
  int c,
  MATRIX& dHao_dx, MATRIX& dHao_dy, MATRIX& dHao_dz, 
  MATRIX& dSao_dx, MATRIX& dSao_dy, MATRIX& dSao_dz
){
/**
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[out] Hao The pointer to the matrix object in which the core Hamiltonian will be stored
  \param[out] Sao The pointer to the AO overlap matrix computed here
  \param[in] DF Debug flag - controlls how much of extra info is printed out
  \param[in] c The index of the atom w.r.t. which coordinates we take the derivatives
  \param[out] dHao_dx The derivative of the Hamiltonian w.r.t. the x-coordinate of the selected atom
  \param[out] dHao_dy The derivative of the Hamiltonian w.r.t. the y-coordinate of the selected atom
  \param[out] dHao_dz The derivative of the Hamiltonian w.r.t. the z-coordinate of the selected atom
  \param[out] dSao_dx The derivative of the AO overlap matrix w.r.t. the x-coordinate of the selected atom
  \param[out] dSao_dy The derivative of the AO overlap matrix w.r.t. the y-coordinate of the selected atom
  \param[out] dSao_dz The derivative of the AO overlap matrix w.r.t. the z-coordinate of the selected atom
  
  Compute the core INDO Hamiltonian and its derivatives w.r.t. specified nuclear DOFs.- Python-friendly version
*/


  Hamiltonian_core_deriv_indo
  ( syst, basis_ao, prms, modprms,  atom_to_ao_map, ao_to_atom_map,  &Hao, &Sao, DF, c,
    &dHao_dx, &dHao_dy, &dHao_dz,   &dSao_dx, &dSao_dy, &dSao_dz);

}





void get_integrals(int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij){
/**
  An auxiliary function: Compute Coulomb and exchange integrals for the orbitals i and j (global indices), both of which 
  are centered on the same atom a (global index)
  eri[a][a] is taken as input argument
  parameters G1 and F2 (Slater-Condon) are taken as input

  \param[in] i Index of one AO   
  \param[in] j Index of another AO   
  \param[in] basis_ao The list of the atomic orbitals - the AO basis
  \param[in] eri_aa On-site electron repulsion integral taken as a parameter
  \param[in] G1 Slater-Condon parameter
  \param[in] F2 Slater-Condon parameter
  \param[out] ii_jj The Coulomb intergal of the AOs i and j
  \param[out] ij_ij The exchange intergal of the AOs i and j
*/

  //=====================================================================================================
  // Integrals:
  ij_ij = ii_jj = 0.0; 

  if( basis_ao[i].ao_shell_type=="s" && basis_ao[j].ao_shell_type=="s"){ 
    ij_ij = ii_jj = eri_aa;                  // ss_ss 
  }

  else if( basis_ao[i].ao_shell_type=="s" && basis_ao[j].ao_shell_type=="p"){
    ij_ij = G1/3.0;                          // sx_sx = sy_sy = sz_sz
    ii_jj = eri_aa;                          // ss_xx = ss_yy = ss_zz
  }
  else if( basis_ao[j].ao_shell_type=="s" && basis_ao[i].ao_shell_type=="p"){
    ij_ij = G1/3.0;                          // xs_xs = ys_ys = zs_zs
    ii_jj = eri_aa;                          // xx_ss = yy_ss = zz_ss
  }
  else if( basis_ao[i].ao_shell_type=="p" && basis_ao[j].ao_shell_type=="p" ){ 
    if( (basis_ao[i].x_exp == basis_ao[j].x_exp) && 
        (basis_ao[i].y_exp == basis_ao[j].y_exp) && 
        (basis_ao[i].z_exp == basis_ao[j].z_exp)){
      ij_ij = ii_jj = eri_aa + 0.16*F2;      // xx_xx = yy_yy = zz_zz
    }
    else{
      ij_ij = 0.12*F2;                         // xy_xy = xz_xz = ...
      ii_jj = eri_aa - 0.08*F2;                // xx_yy = xx_zz = ...
    }
  }

  //=====================================================================================================

}


void Hamiltonian_Fock_indo(Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms, Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          ){
/**
  \param[in,out] el The electronic structre properties of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  
  Compute the INDO (or CNDO/2) Fock Hamiltonian. Unrestricted formulation
*/


  int i,j,k,n,I,J,K,a,b,A,B;

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->n_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }


  // Formation of the Fock matrix: Core part
  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;


  // Update charges
  *el->P = *el->P_alp + *el->P_bet;


  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);


  vector<double> Zeff(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_gross(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_net(syst.Number_of_atoms, 0.0);

  for(a=0;a<syst.Number_of_atoms;a++){ Zeff[a] = modprms.PT[syst.Atoms[a].Atom_element].Zeff; } // e.g. 4 for STO-3G C

  update_Mull_charges(ao_to_atom_map, Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net, Mull_charges_gross, Mull_charges_net);

  for(a=0;a<syst.Number_of_atoms;a++){ 
    syst.Atoms[a].Atom_mull_charge_gross = Mull_charges_gross[a]; 
    syst.Atoms[a].Atom_mull_charge_net = Mull_charges_net[a]; 
  }



    
  // Formation of the Fock matrix: add Coulomb and Exchange parts    
  for(i=0;i<Norb;i++){
    a = ao_to_atom_map[i];

    for(j=0;j<Norb;j++){
      b = ao_to_atom_map[j];

      if(i==j){  // Diagonal terms

//        cout<<"Diagonal terms, i= "<<i<<endl;
//        cout<<"el->Fao_alp->M[i*Norb+i] = "<<el->Fao_alp->M[i*Norb+i]<<endl;

        for(int kk=0;kk<atom_to_ao_map[a].size();kk++){    // for all orbitals on atom a
          k = atom_to_ao_map[a][kk];                       // global orbital index of AO kk on atom a


          double ii_kk, ik_ik; ii_kk = ik_ik = 0.0;
          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms + a];

          get_integrals(i,k,basis_ao,eri_aa,G1,F2,ii_kk,ik_ik);

/*
          cout<<"kk= "<<kk<<" k= "<<k
              <<" basis_ao[i].ao_shell_type= "<<basis_ao[i].ao_shell_type
              <<" basis_ao[k].ao_shell_type= "<<basis_ao[k].ao_shell_type
              <<"basis_ao[i].exponents = ("<<basis_ao[i].x_exp<<","<<basis_ao[i].y_exp<<","<<basis_ao[i].z_exp
              <<"basis_ao[k].exponents = ("<<basis_ao[k].x_exp<<","<<basis_ao[k].y_exp<<","<<basis_ao[k].z_exp
              <<" el->P->M[k*Norb+k]= "<<el->P->M[k*Norb+k]
              <<" el->P_alp->M[k*Norb+k]= "<<el->P_alp->M[k*Norb+k]
              <<" ii_kk= "<<ii_kk<<" ik_ik= "<<ik_ik
              <<endl;
*/

          if(prms.use_rosh){ // Restricted open-shell
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);

          }
          else{ // unrestricted
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_alp->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_bet->M[k*Norb+k]*ik_ik);

          }

//          cout<<"at this point el->Fao_bet->M[i*Norb+i]= "<<el->Fao_bet->M[i*Norb+i]<<endl;
  
        }// for kk - all orbitals on atom A
  

        // Contributions from all other atoms to the diagonal terms
        // don't worry that b determined above will be rewritten - this is ok
        for(b=0;b<syst.Number_of_atoms;b++){

          if(b!=a){

            // Compute density matrix due to all orbitals on atom b 
            el->Fao_alp->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 
          }
        }// for b
  
  
      }// i==j
      else{      // Off-diagonal terms
  
        if(a==b){ // different orbitals are on the same atom
          double ij_ij,ii_jj; ij_ij = ii_jj = 0.0;

          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms+a];
          get_integrals(i,j,basis_ao,eri_aa,G1,F2,ii_jj,ij_ij);

          //cout<<"i= "<<i<<" j= "<<j<<" ii_jj= "<<ii_jj<<", ij_ij= "<<ij_ij<<endl;

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
          }
          else{  
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_alp->M[i*Norb+j])*ij_ij - el->P_alp->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_bet->M[i*Norb+j])*ij_ij - el->P_bet->M[i*Norb+j]*ii_jj );
          }
  
        }// a==b - different orbitals are on the same atom
  
        else{ // different orbitals are on different atoms

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
          }

        }
  
  
      }// i!=j
  
  
    }// for j
  }// for i

}


void Hamiltonian_Fock_indo(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
                           Control_Parameters& prms, Model_Parameters& modprms,
                           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
                          ){
/**
  \param[in,out] el The electronic structre properties of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  
  Compute the INDO (or CNDO/2) Fock Hamiltonian. Unrestricted formulation - Python-friendly version
*/


  Hamiltonian_Fock_indo(&el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map);

}




void Hamiltonian_Fock_derivs_indo
( Electronic_Structure* el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX* dHao_dx, MATRIX* dHao_dy, MATRIX* dHao_dz,
  MATRIX* dFao_alp_dx, MATRIX* dFao_alp_dy, MATRIX* dFao_alp_dz,
  MATRIX* dFao_bet_dx, MATRIX* dFao_bet_dy, MATRIX* dFao_bet_dz
){
/**
  \param[in,out] el The electronic structre properties of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] c The index of the atom w.r.t. which coordinates we take the derivatives
  \param[in] dHao_dx The derivative of the core Hamiltonian w.r.t. the x-coordinate of the selected atom
  \param[in] dHao_dy The derivative of the core Hamiltonian w.r.t. the y-coordinate of the selected atom
  \param[in] dHao_dz The derivative of the core Hamiltonian w.r.t. the z-coordinate of the selected atom
  \param[out] dFao_alp_dx The derivative of the Fock Hamiltonian (alpha-component) w.r.t. the x-coordinate of the selected atom
  \param[out] dFao_alp_dy The derivative of the Fock Hamiltonian (alpha-component) w.r.t. the y-coordinate of the selected atom
  \param[out] dFao_alp_dz The derivative of the Fock Hamiltonian (alpha-component) w.r.t. the z-coordinate of the selected atom
  \param[out] dFao_bet_dx The derivative of the Fock Hamiltonian (beta-component) w.r.t. the x-coordinate of the selected atom
  \param[out] dFao_bet_dy The derivative of the Fock Hamiltonian (beta-component) w.r.t. the y-coordinate of the selected atom
  \param[out] dFao_bet_dz The derivative of the Fock Hamiltonian (beta-component) w.r.t. the z-coordinate of the selected atom
  
  Compute the INDO (or CNDO/2) Fock Hamiltonian as well as the gradients of the Fock matrix. Unrestricted formulation
*/

  int i,j,k,n,I,J,K,a,b,A,B;

  Timer tim1;

  int Norb = basis_ao.size(); // how many AOs are included in this fragment
  if(Norb!=el->Hao->n_cols){  
    cout<<"In Hamiltonian_Fock_indo: Dimension of input/output matrix is not compatible whith the number of the fragment-localized orbitals\n";
    exit(0);
  }

  stringstream ss(stringstream::in | stringstream::out);
  std::string out;
  (ss << c);  ss >> out;

  int use_disk = 1;  

  MATRIX* aux1;  
  MATRIX* aux2;  
  MATRIX* aux3;  
  MATRIX* aux4;  
  MATRIX* aux5;  
  MATRIX* aux6;  

  int nat = syst.Number_of_atoms;
  if(use_disk){

    aux1 = new MATRIX(nat,nat);
    aux2 = new MATRIX(nat,nat);
    aux3 = new MATRIX(nat,nat);
    aux4 = new MATRIX(nat,nat);
    aux5 = new MATRIX(nat,nat);
    aux6 = new MATRIX(nat,nat);

    // Load matrices with derivatives of the core parameters
    aux1->bin_load("deri.x_"+out+".bin");
    aux2->bin_load("deri.y_"+out+".bin");
    aux3->bin_load("deri.z_"+out+".bin");

    aux4->bin_load("dV_AB.x_"+out+".bin");
    aux5->bin_load("dV_AB.y_"+out+".bin");
    aux6->bin_load("dV_AB.z_"+out+".bin");

  }




  // Formation of the Fock matrix: Core part
  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;

  *dFao_alp_dx = *dHao_dx;
  *dFao_alp_dy = *dHao_dy;
  *dFao_alp_dz = *dHao_dz;

  *dFao_bet_dx = *dHao_dx;
  *dFao_bet_dy = *dHao_dy;
  *dFao_bet_dz = *dHao_dz;


  vector<int> sorb_indx;
  sorb_indx = compute_sorb_indices(syst.Number_of_atoms, basis_ao, atom_to_ao_map, ao_to_atom_map);


  // Update charges
  *el->P = *el->P_alp + *el->P_bet;


  update_Mull_orb_pop(el->P, el->Sao, el->Mull_orb_pop_gross, el->Mull_orb_pop_net);


  vector<double> Zeff(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_gross(syst.Number_of_atoms, 0.0);
  vector<double> Mull_charges_net(syst.Number_of_atoms, 0.0);

  for(a=0;a<syst.Number_of_atoms;a++){ Zeff[a] = modprms.PT[syst.Atoms[a].Atom_element].Zeff; } // e.g. 4 for STO-3G C

  update_Mull_charges(ao_to_atom_map, Zeff, el->Mull_orb_pop_gross, el->Mull_orb_pop_net, Mull_charges_gross, Mull_charges_net);

  for(a=0;a<syst.Number_of_atoms;a++){ 
    syst.Atoms[a].Atom_mull_charge_gross = Mull_charges_gross[a]; 
    syst.Atoms[a].Atom_mull_charge_net = Mull_charges_net[a]; 
  }



    
  // Formation of the Fock matrix: add Coulomb and Exchange parts    
  for(i=0;i<Norb;i++){
    a = ao_to_atom_map[i];

    for(j=0;j<Norb;j++){
      b = ao_to_atom_map[j];

      if(i==j){  // Diagonal terms

        for(int kk=0;kk<atom_to_ao_map[a].size();kk++){    // for all orbitals on atom a
          k = atom_to_ao_map[a][kk];                       // global orbital index of AO kk on atom a


          double ii_kk, ik_ik; ii_kk = ik_ik = 0.0;
          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms + a];

          get_integrals(i,k,basis_ao,eri_aa,G1,F2,ii_kk,ik_ik);


          if(prms.use_rosh){ // Restricted open-shell
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - 0.5*el->P->M[k*Norb+k]*ik_ik);

          }
          else{ // unrestricted
            el->Fao_alp->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_alp->M[k*Norb+k]*ik_ik);
            el->Fao_bet->M[i*Norb+i] += (el->P->M[k*Norb+k]*ii_kk - el->P_bet->M[k*Norb+k]*ik_ik);

          }
  
  
        }// for kk - all orbitals on atom A
  

        // Contributions from all other atoms to the diagonal terms
        // don't worry that b determined above will be rewritten - this is ok
        for(b=0;b<syst.Number_of_atoms;b++){

          if(b!=a){

            // Compute density matrix due to all orbitals on atom b 
            el->Fao_alp->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*modprms.eri[a*syst.Number_of_atoms+b]; 


            // Derivatives

            VECTOR deri, dV_AB; 

            tim1.start();
            if(use_disk){
              deri.x = aux1->get(a,b);
              deri.y = aux2->get(a,b);
              deri.z = aux3->get(a,b);
            }
            else{
              compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                                  sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);
            }

            tim1.stop();

            dFao_alp_dx->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.x;
            dFao_alp_dy->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.y;
            dFao_alp_dz->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.z;

            dFao_bet_dx->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.x;
            dFao_bet_dy->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.y;
            dFao_bet_dz->M[i*Norb+i] += (Zeff[b] - syst.Atoms[b].Atom_mull_charge_net)*deri.z;


          }
        }// for b
  
  
      }// i==j
      else{      // Off-diagonal terms
  
        if(a==b){ // different orbitals are on the same atom
          double ij_ij,ii_jj; ij_ij = ii_jj = 0.0;

          double G1 = modprms.PT[basis_ao[i].element].G1[basis_ao[i].ao_shell];
          double F2 = modprms.PT[basis_ao[i].element].F2[basis_ao[i].ao_shell];
          double eri_aa = modprms.eri[a*syst.Number_of_atoms+a];

          get_integrals(i,j,basis_ao,eri_aa,G1,F2,ii_jj,ij_ij);

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - 0.5*el->P->M[i*Norb+j])*ij_ij - 0.5*el->P->M[i*Norb+j]*ii_jj );
          }
          else{  
            el->Fao_alp->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_alp->M[i*Norb+j])*ij_ij - el->P_alp->M[i*Norb+j]*ii_jj );
            el->Fao_bet->M[i*Norb+j] += ( (2.0*el->P->M[i*Norb+j] - el->P_bet->M[i*Norb+j])*ij_ij - el->P_bet->M[i*Norb+j]*ii_jj );
          }
  
        }// a==b - different orbitals are on the same atom
  
        else{ // different orbitals are on different atoms

          if(prms.use_rosh){
            el->Fao_alp->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 

            tim1.start();
            VECTOR deri, dV_AB;    

            if(use_disk){
              deri.x = aux1->get(a,b);
              deri.y = aux2->get(a,b);
              deri.z = aux3->get(a,b);
            }
            else{
              compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                                  sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);
            }
            tim1.stop();

            dFao_alp_dx->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.x; 
            dFao_alp_dy->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.y; 
            dFao_alp_dz->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.z; 

            dFao_bet_dx->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.x; 
            dFao_bet_dy->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.y; 
            dFao_bet_dz->M[i*Norb+j] -= 0.5*el->P->M[i*Norb+j]*deri.z; 

          }
          else{
            el->Fao_alp->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 
            el->Fao_bet->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*modprms.eri[a*syst.Number_of_atoms+b]; 


            tim1.start();
            VECTOR deri, dV_AB;    
            if(use_disk){
              deri.x = aux1->get(a,b);
              deri.y = aux2->get(a,b);
              deri.z = aux3->get(a,b);
            }
            else{
              compute_indo_core_parameters_derivs(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map,
                                                  sorb_indx, modprms.indo_opt, a, b, c, deri, dV_AB);
            }
            tim1.stop();

            dFao_alp_dx->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*deri.x; 
            dFao_alp_dy->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*deri.y; 
            dFao_alp_dz->M[i*Norb+j] -= el->P_alp->M[i*Norb+j]*deri.z; 

            dFao_bet_dx->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*deri.x; 
            dFao_bet_dy->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*deri.y; 
            dFao_bet_dz->M[i*Norb+j] -= el->P_bet->M[i*Norb+j]*deri.z; 

          }

        }
    
      }// i!=j
  
  
    }// for j
  }// for i

  if(use_disk){
    delete aux1; delete aux2; delete aux3;
    delete aux4; delete aux5; delete aux6;
  }


//  cout<<"End of Hamiltonian_Fock_derivs_indo: Time to compute core_parameters_derivs = "<<tim1.show()<<endl;

}

void Hamiltonian_Fock_derivs_indo
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int c, 
  MATRIX& dHao_dx,     MATRIX& dHao_dy,     MATRIX& dHao_dz,
  MATRIX& dFao_alp_dx, MATRIX& dFao_alp_dy, MATRIX& dFao_alp_dz,
  MATRIX& dFao_bet_dx, MATRIX& dFao_bet_dy, MATRIX& dFao_bet_dz
){
/**
  \param[in,out] el The electronic structre properties of the system
  \param[in] syst The object defining molecular structure of the chemical system
  \param[in] basis_ao The vector of AO objects - it constitutes the atomic basis of the system
  \param[in] prms The parameters controlling the quantum mechanical calculations
  \param[in,out] modprms The parameters of the atomistic Hamiltonian
  \param[in] atom_to_ao_map The mapping from the atomic indices to the lists of the indices of AOs localized on given atom
  \param[in] ao_to_atom_map The mapping from the AO index to the index of atoms on which given AO is localized
  \param[in] c The index of the atom w.r.t. which coordinates we take the derivatives
  \param[in] dHao_dx The derivative of the core Hamiltonian w.r.t. the x-coordinate of the selected atom
  \param[in] dHao_dy The derivative of the core Hamiltonian w.r.t. the y-coordinate of the selected atom
  \param[in] dHao_dz The derivative of the core Hamiltonian w.r.t. the z-coordinate of the selected atom
  \param[out] dFao_alp_dx The derivative of the Fock Hamiltonian (alpha-component) w.r.t. the x-coordinate of the selected atom
  \param[out] dFao_alp_dy The derivative of the Fock Hamiltonian (alpha-component) w.r.t. the y-coordinate of the selected atom
  \param[out] dFao_alp_dz The derivative of the Fock Hamiltonian (alpha-component) w.r.t. the z-coordinate of the selected atom
  \param[out] dFao_bet_dx The derivative of the Fock Hamiltonian (beta-component) w.r.t. the x-coordinate of the selected atom
  \param[out] dFao_bet_dy The derivative of the Fock Hamiltonian (beta-component) w.r.t. the y-coordinate of the selected atom
  \param[out] dFao_bet_dz The derivative of the Fock Hamiltonian (beta-component) w.r.t. the z-coordinate of the selected atom
  
  Compute the INDO (or CNDO/2) Fock Hamiltonian as well as the gradients of the Fock matrix. Unrestricted formulation - Python-friendly version
*/


  Hamiltonian_Fock_derivs_indo( &el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map,
  c, &dHao_dx, &dHao_dy, &dHao_dz,  &dFao_alp_dx, &dFao_alp_dy, &dFao_alp_dz,  &dFao_bet_dx, &dFao_bet_dy, &dFao_bet_dz);

}





}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra
