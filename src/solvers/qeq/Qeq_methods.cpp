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

#include "Qeq.h"


void Qeq_Init(std::string idempotential,map<std::string,_GMP>& GMP){

// Set up (or calculate) some required atomic properties
  _GMP gmp;

//=============== Setting up properties of atoms =================
/*  Literature references:

[1] - Rappe, A. K., Goddard III W. A. "Charge Equilibration for Molecular Dynamics
      Simulations" J. Phys. Chem. 1991 V. 95, P. 3358-3363

[2] - Parr, R. G., Pearson, R. G. "Absolute Hardness: Companion Parameter to Absolute
      Electronegativity" J. Am. Chem. Soc. 1983 V. 105, P. 7512-7516

[3] - Bakowies, D., Thiel, W. "Semiempirical Treatment of Electrostatic Potentials
      and Partial Charges in Combined Quantum Mechanical and Molecular Mechanical
      Approiaches" J. Comp. Chem. 1996 V. 17, P. 87-108

[4] - Hans Martin Senn "Transition-Metal-Catalysed Hydroamination of Alkenes:
      Theorethical Studies Using Ab Initio Molecular Dynamics" Dipl. ETH No. 13972
      (he gives values obtained from Rappe&Goddard III).
*/

  if(idempotential=="rappe_goddard"){ // Based on [1]

  // atomic radii - in Angstrems from Rappe, Goddard
  // ksi          - in inverse a.u. from Rappe, Goddard
  // J            - in eV from Rappe, Goddard
  // xi           - in eV from Rappe, Goddard

  gmp.xi  = 4.528; gmp.type = "H"; gmp.r = 0.371; gmp.n = 1; gmp.ksi = 1.0698; gmp.J = 13.8904; GMP["H"] = gmp;
  gmp.xi  = 5.110; gmp.type = "B"; gmp.r = 0.838; gmp.n = 2; gmp.ksi = 0.7458; gmp.J = 8.592  ; GMP["B"] = gmp;
  gmp.xi  = 5.343; gmp.type = "C"; gmp.r = 0.759; gmp.n = 2; gmp.ksi = 0.8563; gmp.J = 10.126 ; GMP["C"] = gmp;
  gmp.xi  = 6.899; gmp.type = "N"; gmp.r = 0.715; gmp.n = 2; gmp.ksi = 0.9089; gmp.J = 11.760 ; GMP["N"] = gmp;
  gmp.xi  = 8.741; gmp.type = "O"; gmp.r = 0.669; gmp.n = 2; gmp.ksi = 0.9745; gmp.J = 13.364 ; GMP["O"] = gmp;
  gmp.xi  = 4.168; gmp.type = "Si";gmp.r = 1.176; gmp.n = 3; gmp.ksi = 0.7737; gmp.J =  6.974 ; GMP["Si"]= gmp;
  gmp.xi  = 5.463; gmp.type = "P"; gmp.r = 1.102; gmp.n = 3; gmp.ksi = 0.8257; gmp.J =  8.000 ; GMP["P"] = gmp;
  gmp.xi  = 6.928; gmp.type = "S"; gmp.r = 1.047; gmp.n = 3; gmp.ksi = 0.8690; gmp.J =  8.972  ; GMP["S"] = gmp;
  gmp.xi  =10.874; gmp.type = "F"; gmp.r = 0.706; gmp.n = 2; gmp.ksi = 0.9206; gmp.J = 14.948 ; GMP["F"] = gmp;
  gmp.xi  = 8.564; gmp.type = "Cl";gmp.r = 0.994; gmp.n = 3; gmp.ksi = 0.9154; gmp.J =  9.892 ; GMP["Cl"]= gmp;
  gmp.xi  = 7.790; gmp.type = "Br";gmp.r = 1.141; gmp.n = 4; gmp.ksi = 1.0253; gmp.J =  8.850 ; GMP["Br"]= gmp;
  gmp.xi  = 6.822; gmp.type = "I"; gmp.r = 1.333; gmp.n = 5; gmp.ksi = 1.0726; gmp.J =  7.524 ; GMP["I"] = gmp;
  gmp.xi  = 3.006; gmp.type = "Li";gmp.r = 1.557; gmp.n = 2; gmp.ksi = 0.4174; gmp.J =  4.772 ; GMP["Li"]= gmp;
  gmp.xi  = 2.843; gmp.type = "Na";gmp.r = 2.085; gmp.n = 3; gmp.ksi = 0.4364; gmp.J =  4.592 ; GMP["Na"]= gmp;
  gmp.xi  = 2.421; gmp.type = "K"; gmp.r = 2.586; gmp.n = 4; gmp.ksi = 0.4524; gmp.J =  3.840 ; GMP["K"] = gmp;
  gmp.xi  = 2.331; gmp.type = "Rb";gmp.r = 2.770; gmp.n = 5; gmp.ksi = 0.5162; gmp.J =  3.692 ; GMP["Rb"]= gmp;
  gmp.xi  = 2.183; gmp.type = "Cs";gmp.r = 2.984; gmp.n = 6; gmp.ksi = 0.5663; gmp.J =  3.422 ; GMP["Cs"]= gmp;

  }
  else if(idempotential=="parr_pearson"){// based on [2]
  // J            - in eV
  // xi           - in eV

  gmp.type  = "H";  gmp.xi = 7.18; gmp.J = 2.0*6.43; GMP["H"]  = gmp;
  gmp.type  = "B";  gmp.xi = 4.29; gmp.J = 2.0*4.01; GMP["B"]  = gmp;
  gmp.type  = "C";  gmp.xi = 6.27; gmp.J = 2.0*5.00; GMP["C"]  = gmp;
  gmp.type  = "N";  gmp.xi = 7.30; gmp.J = 2.0*7.23; GMP["N"]  = gmp;
  gmp.type  = "O";  gmp.xi = 7.54; gmp.J = 2.0*6.08; GMP["O"]  = gmp;
  gmp.type  = "Si"; gmp.xi = 4.77; gmp.J = 2.0*3.38; GMP["Si"] = gmp;
  gmp.type  = "P";  gmp.xi = 5.62; gmp.J = 2.0*4.88; GMP["P"]  = gmp;
  gmp.type  = "S";  gmp.xi = 6.22; gmp.J = 2.0*4.14; GMP["S"]  = gmp;
  gmp.type  = "F";  gmp.xi =10.41; gmp.J = 2.0*7.01; GMP["F"]  = gmp;
  gmp.type  = "Cl"; gmp.xi = 8.30; gmp.J = 2.0*4.68; GMP["Cl"] = gmp;
  gmp.type  = "Br"; gmp.xi = 7.59; gmp.J = 2.0*4.22; GMP["Br"] = gmp;
  gmp.type  = "I";  gmp.xi = 6.76; gmp.J = 2.0*3.69; GMP["I"]  = gmp;
  gmp.type  = "Li"; gmp.xi = 3.01; gmp.J = 2.0*2.39; GMP["Li"] = gmp;
  gmp.type  = "Na"; gmp.xi = 2.85; gmp.J = 2.0*2.30; GMP["Na"] = gmp;
  gmp.type  = "K";  gmp.xi = 2.42; gmp.J = 2.0*1.92; GMP["K"]  = gmp;
  gmp.type  = "Rb"; gmp.xi = 2.34; gmp.J = 2.0*1.85; GMP["Rb"] = gmp;
  gmp.type  = "Cs"; gmp.xi = 2.18; gmp.J = 2.0*1.71; GMP["Cs"] = gmp;

  gmp.type  = "Al"; gmp.xi = 3.21; gmp.J = 2.0*2.77; GMP["Al"] = gmp;
  gmp.type  = "V";  gmp.xi = 3.64; gmp.J = 2.0*3.11; GMP["V"]  = gmp;
  gmp.type  = "Cr"; gmp.xi = 3.76; gmp.J = 2.0*3.05; GMP["Cr"] = gmp;
  gmp.type  = "Fe"; gmp.xi = 4.03; gmp.J = 2.0*3.87; GMP["Fe"] = gmp;
  gmp.type  = "Co"; gmp.xi = 4.26; gmp.J = 2.0*3.60; GMP["Co"] = gmp;
  gmp.type  = "Ni"; gmp.xi = 4.44; gmp.J = 2.0*3.24; GMP["Ni"] = gmp;
  gmp.type  = "Cu"; gmp.xi = 4.48; gmp.J = 2.0*3.25; GMP["Cu"] = gmp;
  gmp.type  = "Se"; gmp.xi = 3.89; gmp.J = 2.0*3.86; GMP["Se"] = gmp;
  gmp.type  = "Zr"; gmp.xi = 3.63; gmp.J = 2.0*3.21; GMP["Zr"] = gmp;
  gmp.type  = "Nb"; gmp.xi = 3.88; gmp.J = 2.0*2.99; GMP["Nb"] = gmp;
  gmp.type  = "Mo"; gmp.xi = 3.92; gmp.J = 2.0*3.17; GMP["Mo"] = gmp;
  gmp.type  = "Rh"; gmp.xi = 4.30; gmp.J = 2.0*3.16; GMP["Rh"] = gmp;
  gmp.type  = "Pd"; gmp.xi = 4.44; gmp.J = 2.0*3.88; GMP["Pd"] = gmp;
  gmp.type  = "Ag"; gmp.xi = 4.44; gmp.J = 2.0*3.14; GMP["Ag"] = gmp;
  gmp.type  = "Sn"; gmp.xi = 4.30; gmp.J = 2.0*3.05; GMP["Sn"] = gmp;
  gmp.type  = "Sb"; gmp.xi = 4.84; gmp.J = 2.0*3.79; GMP["Sb"] = gmp;
  gmp.type  = "Te"; gmp.xi = 5.49; gmp.J = 2.0*3.52; GMP["Te"] = gmp;
  gmp.type  = "Ba"; gmp.xi = 2.60; gmp.J = 2.0*2.60; GMP["Ba"] = gmp;
  gmp.type  = "Pt"; gmp.xi = 5.60; gmp.J = 2.0*3.50; GMP["Pt"] = gmp;
  gmp.type  = "Au"; gmp.xi = 5.80; gmp.J = 2.0*3.50; GMP["Au"] = gmp;

  }
  else if(idempotential=="bakowies_thiel"){ // Potential Derived - [3]

  gmp.type = "H"; gmp.xi  = 4.42211; gmp.J = 13.84036; GMP["H"] = gmp;
  gmp.type = "C"; gmp.xi  = 5.07305; gmp.J = 10.06444; GMP["C"] = gmp;
  gmp.type = "N"; gmp.xi  = 7.73699; gmp.J = 12.96908; GMP["N"] = gmp;
  gmp.type = "O"; gmp.xi  = 8.27885; gmp.J = 14.93241; GMP["O"] = gmp;

  }

  else{
  std::cout<<"Error: Possible options for getting idempotentials are: rappe_goddard, parr_pearson, bakowies_thiel "<<std::endl;
  exit(77);
  }

}

void Solve(_MOL& mol,double epsilon,int mode){

  // This is version with iterative Gauss-Seidel solution

  int sz = mol.Q.size();

  double Iij,I0j,p1,pj;
  double a,b,c,Rij,R0j;
  double ksi0_i,ksi0_j;
  double totQ = 0.0;
  int st,st1,st2,i,j;

  //-----------------------------
  MATRIX C(sz,sz);
  for(i=0;i<sz;i++){C.M[0*sz+i] = 1.0;}
  // Calculate matrix C elements
  for(j=0;j<sz;j++){
    R0j = (mol.R[0]-mol.R[j]).length();
    I0j = Coulomb_Integral(R0j,0,j,mol.Q[0],mol.Q[j],mol.GMP,epsilon,mode);
    for(i=1;i<sz;i++){
      Rij = (mol.R[i]-mol.R[j]).length();
      Iij = Coulomb_Integral(Rij,i,j,mol.Q[i],mol.Q[j],mol.GMP,epsilon,mode);
      C.M[i*sz+j]= I0j - Iij;
    }// for i
  }// for j

  //--------------------------------------------------------
  // Calculate matrix D elements (differences of electronegativities)
  MATRIX D(sz,1);  D.M[0] = 0.0;
  for(i=1;i<sz;i++){   D.M[i] = (mol.GMP[i].xi - mol.GMP[0].xi);  }
  //----------------------------------------------------------
  MATRIX q(sz,1);
  // Iterative solution
  solve_linsys(C,D,q,1e-10,100);
  // Inserting the result into mol structure
  for(i=0;i<sz;i++){mol.Q[i] = q.M[i]; }

}

/*****************************

 TO BE MOVED TO OBJECT_SPACE METHODS !!!

int CQeq_Calculate(ObjectSpace& os){

  _MOL mol;
  map<std::string,_GMP>::iterator it;

  // Initialize _MOL structure: coordinates of atoms, atom elements, initial atomic charges and GMP properties

  for(int i=0;i<os.ObjectSpace_Number_of_atoms;i++){
    if(os.ObjectSpace_Atoms[i].is_Atom_element){
      mol.type.push_back(os.ObjectSpace_Atoms[i].Atom_element);
      it = os.GMP.find(os.ObjectSpace_Atoms[i].Atom_element);
      if(it!=os.GMP.end()){
        mol.GMP.push_back(it->second);
      }else{
        std::cout<<"Error: GMP properties are not defined for atom "<<os.ObjectSpace_Atoms[i].Atom_element<<std::endl;
        exit(77);
      }
    }
    else{std::cout<<"Error: Atom element is not defined for atom indx = "<<os.ObjectSpace_Atoms[i].globAtom_Index<<std::endl;
      exit(77);
    }

    // Note! Convert Angstroms to atomic units
    if(os.ObjectSpace_Atoms[i].is_Atom_coords_cart){mol.R.push_back(os.ObjectSpace_Atoms[i].Atom_coords_cart/0.52917); }
    else{std::cout<<"Error: Atom coordinates are not defined for atom indx = "<<os.ObjectSpace_Atoms[i].globAtom_Index<<std::endl;
      exit(77);
    }
    if(os.ObjectSpace_Atoms[i].is_Atom_charge){mol.Q.push_back(os.ObjectSpace_Atoms[i].Atom_charge); }
    else{                                      mol.Q.push_back(0.0);    }

  }// for i

  //===================== Variables ============================
  int MaxCount;
  double threshold;
  int solution_method; // controls which interative algorithm to use
  int integral_method; // controls how Jab is calculated

  if(os.is_ObjectSpace_CQeq_parameters){
    if(os.ObjectSpace_CQeq_parameters.is_threshold){
      threshold = os.ObjectSpace_CQeq_parameters.threshold;
    }else{ std::cout<<"CQeq threshold is not defined"<<std::endl; exit(69);  }
    if(os.ObjectSpace_CQeq_parameters.is_MaxCount){
      MaxCount = os.ObjectSpace_CQeq_parameters.MaxCount;
    }else{ std::cout<<"CQeq MaxCount is not defined"<<std::endl; exit(69);  }
    if(os.ObjectSpace_CQeq_parameters.is_solution_method){
      solution_method = os.ObjectSpace_CQeq_parameters.solution_method;
    }else{ std::cout<<"CQeq solution_method is not defined"<<std::endl; exit(69);  }
    if(os.ObjectSpace_CQeq_parameters.is_integral_method){
      integral_method = os.ObjectSpace_CQeq_parameters.integral_method;
    }else{ std::cout<<"CQeq integral_method is not defined"<<std::endl; exit(69);  }
  }else{
    std::cout<<"Error: CQeq parameters are not defined!"<<std::endl;
    exit(69);
  }

  double eps;
  double* q1; q1 = NULL;
  double* q2; q2 = NULL;
  int count = 0;
  int sz = mol.R.size();

  q1 = new double[sz];
  q2 = new double[sz];

  //========================= Performing calculations ============================
  for(i=0;i<sz;i++){q1[i] = mol.Q[i];}
  do{
    Solve(mol,os.ObjectSpace_CQeq_parameters.epsilon,integral_method);
    for(i=0;i<sz;i++){ q2[i] = mol.Q[i];}
    eps = 0.0;
    for(i=0;i<sz;i++){eps += fabs(q2[i]-q1[i]); q1[i]=q2[i];}
    count++;
    } while((eps>threshold)&&(count<MaxCount));
  //============================================================
  // Update calculated charges
  for(i=0;i<os.ObjectSpace_Number_of_atoms;i++){
    os.ObjectSpace_Atoms[i].Atom_charge = mol.Q[i];
    os.ObjectSpace_Atoms[i].is_Atom_charge = 1; // just in case
  }
  //---- Update pair parameters --------------
  // Default charges
  for(int p=0;p<os.ObjectSpace_Number_of_frag_pairs;p++){
    int act_pair_indx = p;
    int at1,at2;
    at1 = os.ObjectSpace_Frag_pairs[act_pair_indx].globAtom_Index[0];
    at2 = os.ObjectSpace_Frag_pairs[act_pair_indx].globAtom_Index[1];
    double  ku1,ku2; ku1 = ku2 = 0.0;

    if(os.ObjectSpace_Atoms[at1].is_Atom_charge){ku1 = os.ObjectSpace_Atoms[at1].Atom_charge;}
    if(os.ObjectSpace_Atoms[at2].is_Atom_charge){ku2 = os.ObjectSpace_Atoms[at2].Atom_charge;}

    // Default dielectric constant
    double Dielectric = 1.0;
    if(os.ObjectSpace_Frag_pairs[act_pair_indx].is_Group_pair_dielectric){
      Dielectric = os.ObjectSpace_Frag_pairs[act_pair_indx].Group_pair_dielectric;
    }
    if(os.is_ObjectSpace_dielectric){
      Dielectric = os.ObjectSpace_dielectric;
    }
    os.ObjectSpace_Frag_pairs[act_pair_indx].Group_pair_q1q2 = electric*(ku1*ku2/Dielectric);
    os.ObjectSpace_Frag_pairs[act_pair_indx].is_Group_pair_q1q2 = 1;
  }// for p - all pairs

  if(q1!=NULL){ delete [] q1; }
  if(q2!=NULL){ delete [] q2; }
  return 0;
}

***********************************************/








