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
 \file main.cpp
 \brief Main module of the program

*/


#include "Engine.h"
#include "random.h"
#include "MOAO.h"
#include "classes.h"
#include "Control_Parameters.h"
#include "Model_Parameters.h"
#include "Nuclear.h"
#include "Basis.h"
#include "Electronic.h"
#include "Hamiltonian.h"
#include "SCF.h"
#include "DOS.h"
#include "Print.h"
#include "Memory.h"
#include "Guess.h"
#include "Dipole.h"
#include "Excitations.h"

#include <stdlib.h>
#include <time.h>
#include <sstream>

#include <sys/utsname.h>
#include <unistd.h>
#include <cstdlib>
#include <fstream>



using namespace std;



int main(int argc, char* argv[]){

// 4 main types of obkects/classes
// Parameters: Model (EHT, INDO, CNDO, etc.) and simulation (convergence, options, choice of the model, etc.) parameters
// Molecule: nuclear information, geometry, charges, etc.
// AO basis: AOs - global array of AOs
// Subsystem:  (F)MOs, mapping of global AOs to particular nuclear sysbsystem

// Algorithms - work on subsystems, molecules, etc.

  // Identify platform
  struct utsname un;
  uname (&un);
  cout << "Running on system: <sysname> " << un.sysname << " </sysname>" << endl;
  cout << "Running on node: <nodename> " << un.nodename << " </nodename>" << endl;




  //==============================================================================================
  //============================== Global variables ==============================================  
  // Parameters
  Control_Parameters prms;  // This already loads some default settings
  Model_Parameters modprms; // This already loads some default settings

  // Dynamical variables
  Nuclear mol;             // 0 atoms
  vector<AO> basis_ao;
  vector<vector<int> > at_orbitals; // Mapping

  int Nfrag;
  vector< vector<int> > basis_fo_sad;
  vector< vector<int> > basis_fo;
  vector<Electronic*> el;
  vector<Electronic*> el0;
  
  // transition dipole moments in AO basis
  VECTOR mu_nucl; mu_nucl = 0.0;
  vector<MATRIX*> mux;
  vector<MATRIX*> muy;
  vector<MATRIX*> muz;  
  
 
  // Descriptive variables
  int Nelec;
  int Norb;
  int Nocc;
  int Nocc_alp;
  int Nocc_bet;

  // Auxiliary variables
  vector<Atom> atoms;
  Memory* mem; mem = new Memory(20,40, 20,40); // auxiliary variables for temporary memory storage


  // Debug and verbosity
  int DF = 0;

  // Arbitrary variables
  int i,j,n,a,b;
  time_t ti1,ti2;

  // Timing
  Timer t; 

  //==============================================================================================
  //==============================================================================================



  //==============================================================================================
  // ===================================== Execution part  =======================================


  srand(time(0));

  //---------------------------- Reading input file ----------------------------------
  if(argc<2){ cout<<"Wrong number of command-line arguments\n Must be at least 1\n"; exit(0); }
  std::string filename = argv[1];  


  std::string ftmp;
  // Make all potentially produced files empty, to spot program termination
  ftmp = filename+".xyz";   ofstream out(ftmp.c_str(),ios::out); out.close();
  ftmp = filename+".xyzq";  out.open(ftmp.c_str(),ios::out); out.close();
  ftmp = filename+"0.el_struct";  out.open(ftmp.c_str(),ios::out); out.close();
  ftmp = filename+"0.dipole";  out.open(ftmp.c_str(),ios::out); out.close();

  //----------------------------------------------------------------------------------
  //-------------------------- Build system and set up variables ---------------------
  // Initialize/read control variables
  get_parameters_from_file(filename, prms, atoms);
  DF = prms.DF;

  // Initialize periodic cell (if any)
  Cell cell(prms.t1,prms.t2,prms.t3);

//exit(0);

  // Initialize/read model parameters (need basis info)
  if(prms.hamiltonian=="eht"||prms.hamiltonian=="geht"){   set_parameters_eht(prms, modprms);  }
  else if(prms.hamiltonian=="indo"){  set_parameters_indo(prms, modprms); }
  else if(prms.hamiltonian=="geht1"){ set_parameters_geht1(prms, modprms); }
  else if(prms.hamiltonian=="geht2"){ set_parameters_geht2(prms, modprms); }

//exit(0);

  // Initialize nuclear subsystem
  for(n=0;n<atoms.size();n++){

    std::string elt = atoms[n].atom_type;
    double zeff = modprms.PT[elt].Zeff;

    mol.add_atom(elt, elt, modprms.PT[elt].Z, zeff, modprms.PT[elt].mass, atoms[n].R);

  }// for n



  //----------------------------------------------------------------------------------
  //------------ Initialize basis of AOs and set atom-orbital mapping ----------------
  set_basis_STO_3G_DZ(mol, modprms, basis_ao, Nelec, Norb);

//exit(0);

  if(prms.hamiltonian=="eht"||prms.hamiltonian=="geht"||prms.hamiltonian=="geht1"||prms.hamiltonian=="geht2"){  
    set_parameters_eht_mapping(modprms,basis_ao); 
    set_parameters_eht_mapping1(modprms,mol); 

/*
    cout<<"K matrix:\n";
    for(i=0;i<basis_ao.size();i++){
      for(j=0;j<basis_ao.size();j++){
        cout<<modprms.meht_k.get_K_value(i,j)<<" ";
      }
      cout<<"\n";
    }
*/

  }// if

//exit(0);

  cout<<"Showing basis functions\n";
  for(i=0;i<basis_ao.size();i++){
    cout<<"i= "<<i<<" elt= "<<basis_ao[i].element<<" at_indx= "<<basis_ao[i].at_indx<<"  ao_name= "<<basis_ao[i].ao_name<<endl;
  }

  map_atoms_and_orbitals(mol.Nnucl, basis_ao, at_orbitals);

  cout<<"Atom:  {Orbitals}\n";
  show_mapping(at_orbitals);



  //----------------------------------------------------------------------------------
  //------------ Initialize subsystems and set fragment-orbital mapping --------------
  vector<vector<int> > fragments;
  vector<double> charges;

  if(prms.fragments.size()<1){   // Whole system is a single fragment
    vector<int> x;
    for(i=0;i<mol.Nnucl;i++){ x.push_back(i); }
    fragments.push_back(x);
    charges.push_back(prms.charge);
  }
  else{ // Read fregments from input file
    fragments = prms.fragments;
    charges = prms.frag_charge;

    double sum = 0.0;
    for(j=0;j<charges.size();j++){ sum += charges[j]; }
    if(fabs(sum-prms.charge)>0.01){ cout<<"Error: sum of fragment charges is not equal to total charge\n"; exit(0); }
  }

  if(basis_fo.size()>1){ basis_fo.clear(); }
  Nfrag = init_electronic_subsystems(fragments,at_orbitals, basis_fo, el, basis_ao, modprms,mol,charges);


  // Reserve a copy
  el0 = vector<Electronic*>(Nfrag);



  // Memory for INDO
  if(mem->eri.size()>0){  mem->eri.clear(); }
  if(mem->V_AB.size()>0){  mem->V_AB.clear(); }


  cout<<"Fragment:  {Orbitals}\n";
  show_mapping(basis_fo);

  // Allocate memory for transition dipole moments
  mux = vector<MATRIX*>(Nfrag);
  muy = vector<MATRIX*>(Nfrag);
  muz = vector<MATRIX*>(Nfrag);
  for(i=0;i<Nfrag;i++){
    mux[i] = new MATRIX(el[i]->Norb,el[i]->Norb);
    muy[i] = new MATRIX(el[i]->Norb,el[i]->Norb);
    muz[i] = new MATRIX(el[i]->Norb,el[i]->Norb);
  }

//exit(0);

  //----------------------------------------------------------------------------------
  //-------------------- Compute core Hamiltonian and transition dipole in AO basis for all fragments ----------------------
  for(i=0;i<Nfrag;i++){

    // Prepare additional memory for INDO
    if(prms.hamiltonian=="indo"){
      cout<<"Generating INDO core parameters for the fragment i = "<<i<<endl;
      indo_core_parameters(fragments[i], basis_fo[i], basis_ao, mol, mem->eri, mem->V_AB, mem, 1);
    }// if "indo"


    t.start(); cout<<"Computation of Guess density for fragment i = "<<i<<endl;     
    guess(prms, modprms, mol, fragments[i], basis_fo[i], basis_ao, at_orbitals, el[i], mem);
    t.stop();  cout<<"Guess computation for fragment i = "<<i<<" takes "<<t.show()<<" seconds\n";


    if(DF){ 
      cout<<"P[i] = \n"<<*el[i]->P<<endl;
      cout<<"H[i] = \n"<<*el[i]->Hao<<endl;
      cout<<"S[i] = \n"<<*el[i]->Sao<<endl;
    }


    el0[i] = new Electronic(el[i]); // make copy at guess

    // Overlap matrix



    // Compute center of mass of a given fragment
    VECTOR com; com = 0.0;
    double mass = 0.0;
    for(n=0;n<fragments[i].size();n++){
      com += mol.R[fragments[i][n]] * mol.mass[fragments[i][n]];
      mass += mol.mass[fragments[i][n]];
    }
    com = com/mass;
    cout<<"Total mass of the fragment "<<i<<" is: "<<mass<<endl;
    cout<<"Coordinates of the center of mass of the fragment "<<i<<" are: "<<com<<endl;

    // Nuclear contribution to dipole moment
    for(n=0;n<fragments[i].size();n++){
      mu_nucl += mol.Zeff[fragments[i][n]] * (mol.R[fragments[i][n]] - com);
    }

    cout<<"Nuclear contribution to dipole moment = "<<mu_nucl<<endl;

//exit(0);

    // Precompute transition dipole moments in AO basis
    t.start();
    if(prms.compute_dipole){
      update_dipole_matrix(prms.x_period,prms.y_period,prms.z_period,prms.t1,prms.t2,prms.t3,
                           basis_fo[i],basis_ao, mux[i], muy[i], muz[i], mem->aux, mem->n_aux, mol);
    }
    t.stop();
    cout<<"Computation of transition dipole matrix for fragment i = "<<i<<" takes "<<t.show()<<" seconds\n";

    if(DF){
      cout<<"Bare dipole momenta <xi_i|r|xi_j> - not translationally invariant\n";
      cout<<"AO basis: Transition dipole moment x-component for fragment i= "<<i<<"\n"<<*mux[i]<<endl;
      cout<<"AO basis: Transition dipole moment y-component for fragment i= "<<i<<"\n"<<*muy[i]<<endl;
      cout<<"AO basis: Transition dipole moment z-component for fragment i= "<<i<<"\n"<<*muz[i]<<endl;
    }

    // mu = -<i|r-COM|j> + <i| sum[ Z_a * (R_a - COM)] |j>  - add nuclear contribution at the very end !!!
    //                          a
    *mux[i] = *mux[i] + com.x * *el[i]->Sao;// + mu_nucl.x * *el[i]->Sao;
    *muy[i] = *muy[i] + com.y * *el[i]->Sao;// + mu_nucl.y * *el[i]->Sao;
    *muz[i] = *muz[i] + com.z * *el[i]->Sao;// + mu_nucl.z * *el[i]->Sao;

    if(DF){
      cout<<"Dipole momenta <xi_i|r-Rcom(frag)|xi_j> - translationally invariant, but not rotationally invariant\n";
      cout<<"AO basis: Transition dipole moment x-component for fragment i= "<<i<<"\n"<<*mux[i]<<endl;
      cout<<"AO basis: Transition dipole moment y-component for fragment i= "<<i<<"\n"<<*muy[i]<<endl;
      cout<<"AO basis: Transition dipole moment z-component for fragment i= "<<i<<"\n"<<*muz[i]<<endl;
    }


  }// for i - all fragments


  //======================== The heart of the code =======================


  //----------------------------------------------------------------------------------
  //-------------------------- Self-consistent field (SCF) ---------------------------
  if(prms.hamiltonian=="hf"){

    set_parameters_hf(prms,modprms,basis_ao);

  }


  if(prms.runtype=="scf"){

    vector<double> Eelec(Nfrag,0.0);
    vector<double> Enucl(Nfrag,0.0);
 
    for(i=0;i<Nfrag;i++){

      if(prms.hamiltonian=="indo" || prms.hamiltonian=="geht1"|| prms.hamiltonian=="geht2" ){

        if(prms.hamiltonian=="indo"){
          // Update core parameters and core Hamiltonian for this fragment
          indo_core_parameters(fragments[i], basis_fo[i], basis_ao, mol, mem->eri, mem->V_AB, mem, 1);
          Hamiltonian_core_indo(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals, el[i]->Hao, el[i]->Sao, mem, mem->eri, mem->V_AB);
        }
        else if(prms.hamiltonian=="geht1"){
          // Update core parameters and core Hamiltonian for this fragment
          //geht1_core_parameters(fragments[i], basis_fo[i], basis_ao, mol, mem->eri, mem->V_AB, mem, 1);
          Hamiltonian_core_geht1(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals, el[i]->Hao, el[i]->Sao, mem);
//          Hamiltonian_core_geht1(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals, el[i]->Hao, el[i]->Sao, mem, mem->eri, mem->V_AB);
        }
        else if(prms.hamiltonian=="geht2"){
          // Update core parameters and core Hamiltonian for this fragment
          Hamiltonian_core_geht2(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals, el[i]->Hao, el[i]->Sao, mem);
        }



        if(DF){
          // Printing some intermediate results
          MATRIX* eri; eri = new MATRIX(fragments[i].size(),fragments[i].size());
          MATRIX* V_AB; V_AB = new MATRIX(fragments[i].size(),fragments[i].size());

          for(a=0;a<fragments[i].size();a++){
            for(b=0;b<fragments[i].size();b++){
              eri->M[a*fragments[i].size()+b] = mem->eri[a*fragments.size()+b]; 
              V_AB->M[a*fragments[i].size()+b] = mem->V_AB[a*fragments.size()+b]; 
            }
          }// for a
          cout<<"eri = \n"<<*eri<<endl;
          cout<<"V_AB = \n"<<*V_AB<<endl;
          delete eri;
          delete V_AB;
        }// if DF 

      }// if "indo" || "geht1" || "geht2"



      // Initialize Fock matrix
      Hamiltonian_Fock(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);

      if(DF){
        cout<<"After Fock initialization...\n";
        cout<<"Hao = \n"<<*el[i]->Hao<<endl;
        cout<<"P_alp = \n"<<*el[i]->P_alp<<endl;
        cout<<"P_bet = \n"<<*el[i]->P_bet<<endl;
        cout<<"Fao_alp = \n"<<*el[i]->Fao_alp<<endl;
        cout<<"Fao_bet = \n"<<*el[i]->Fao_bet<<endl;
      }

      cout<<"|F|= "<<el[i]->Fao_alp->max_elt()<<endl;
      cout<<"|S|= "<<el[i]->Sao->max_elt()<<endl;




      // Compute SCF for all fragments
      if(prms.use_disk){

        Eelec[i] = scf_oda_disk(prms,modprms,mol,fragments[i],basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);

      }
      else{

        Eelec[i] = scf(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);

      }

      // Compute nuclear part of the energy
      Enucl[i] = Energy_nucl(mol,fragments[i]);

      cout<<"Nuclear energy = "<<Enucl[i]<<endl;


      if(prms.compute_excitations==1){

        compute_excitations(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);
        print_excitations(filename+".spectrum", el[i], prms);

      }// if compute_excitations == 1


      if(prms.compute_dos==1){

        compute_dos(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],mem);

      }// if compute_dos == 1

//      if(prms.compute_vertical_ip==1){
//
//        compute_vertical_ip(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);
//      }


    
    }// for i



    // Now do printing    
    print_xyz(filename+".xyz",mol);     // Print geometry of the entire system
    print_xyzq(filename+".xyzq",mol,Eelec[0],Enucl[0],Eelec[0]+Enucl[0]);   // Print geometry and charges of the entire system

    // Print electronic structure info for all fragments, one by one
    for(i=0;i<Nfrag;i++){

      stringstream ss(stringstream::in | stringstream::out);
      std::string out; int inp = i;  (ss << inp);  ss >> out;

      print_el_struct(filename+out+".el_struct",el[i], Eelec[i],prms.degen_tol);

      if(prms.compute_dipole){
        print_dipole(filename+out+".dipole",el[i], Eelec[i],mux[i],muy[i],muz[i],mu_nucl);
      }


    }// for i



//  Debug
    if(DF){ 
      for(i=0;i<Nfrag;i++){

        MATRIX* Ptmpa; Ptmpa = new MATRIX(el[i]->Norb,el[i]->Norb);
        MATRIX* Ptmpb; Ptmpb = new MATRIX(el[i]->Norb,el[i]->Norb);

        compute_density_matrix(el[i]->occ_alp, el[i]->C_alp, Ptmpa);
        compute_density_matrix(el[i]->occ_bet, el[i]->C_bet, Ptmpb);

        cout<<"Composed density matrix (alpha) for fragment i = "<<i<<endl<<*el[i]->P_alp<<endl;
        cout<<"Density matrix recovered from WFC (alpha): \n"<<(*Ptmpa)<<endl; 

        cout<<"Composed density matrix (beta) for fragment i = "<<i<<endl<<*el[i]->P_bet<<endl;
        cout<<"Density matrix recovered from WFC (beta): \n"<<(*Ptmpb)<<endl; 


        cout<<"Composed density matrix for fragment i = "<<i<<endl<<*el[i]->P<<endl;
        cout<<"Density matrix recovered from WFC: \n"<<(*Ptmpa + *Ptmpb)<<endl; 

/*!
      Verification of certain identities:
      http://arxiv.org/pdf/cond-mat/9407115.pdf
      
      C*C^{+} = S^{-1}
      
      This is verified:

      \code{.cpp}
       cout<<"C*C^+ = \n"<<*el[i]->C_alp *  (*el[i]->C_alp).T()<<endl;
       cout<<"Sao = \n"<<*el[i]->Sao<<endl;
       MATRIX* Sinv;  Sinv = new MATRIX(el[i]->Norb,el[i]->Norb);
       el[i]->Sao->Inverse(Sinv);
       cout<<"Sao_inv = \n"<<*Sinv<<endl;

       cout<<"This gives unity matrix\n";
       cout<<"C^{+} * S * C = \n"<<(*el[i]->C_alp).T() * (*el[i]->Sao) *  (*el[i]->C_alp)<<endl;
       cout<<"This gives diagonal matrix\n";
       cout<<"C^{+} * F * C = \n"<<(*el[i]->C_alp).T() * (*el[i]->Fao_alp) *  (*el[i]->C_alp)<<endl;
       cout<<"This gives occupation matrix\n";
       cout<<"C^{+} * S * P * S * C = \n"<<(*el[i]->C_alp).T() * (*el[i]->Sao) * (*el[i]->P_alp) * (*el[i]->Sao) * (*el[i]->C_alp)<<endl;
       cout<<"FON density matrix need not commute with Fock matrix. This may be non-zero matrix.\n";
       cout<<"F * P * S - S * P * F = \n"<<(*el[i]->Fao_alp) * (*el[i]->P_alp) * (*el[i]->Sao) - (*el[i]->Sao) * (*el[i]->P_alp) * (*el[i]->Fao_alp)<<endl;
       cout<<"Density matrix constructed from eigenstates does commute with Fock matrix. This should be zero matrix.\n";
       cout<<"F * P * S - S * P * F = \n"<<(*el[i]->Fao_alp) * (*Ptmpa) * (*el[i]->Sao) - (*el[i]->Sao) * (*Ptmpa) * (*el[i]->Fao_alp)<<endl;


      \endcode
*/
        cout<<"This gives unity matrix\n";
        cout<<"C^{+} * S * C = \n"<<(*el[i]->C_alp).T() * (*el[i]->Sao) *  (*el[i]->C_alp)<<endl;
        cout<<"This gives diagonal matrix (alpha, energies)\n";
        cout<<"C^{+} * F * C = \n"<<(*el[i]->C_alp).T() * (*el[i]->Fao_alp) *  (*el[i]->C_alp)<<endl;
        cout<<"This gives diagonal matrix (beta, energies)\n";
        cout<<"C^{+} * F * C = \n"<<(*el[i]->C_bet).T() * (*el[i]->Fao_bet) *  (*el[i]->C_bet)<<endl;
        cout<<"This gives occupation matrix (alpha)\n";
        cout<<"C^{+} * S * P * S * C = \n"<<(*el[i]->C_alp).T() * (*el[i]->Sao) * (*el[i]->P_alp) * (*el[i]->Sao) * (*el[i]->C_alp)<<endl;
        cout<<"This gives occupation matrix (beta)\n";
        cout<<"C^{+} * S * P * S * C = \n"<<(*el[i]->C_bet).T() * (*el[i]->Sao) * (*el[i]->P_bet) * (*el[i]->Sao) * (*el[i]->C_bet)<<endl;
        cout<<"FON density matrix need not commute with Fock matrix. This may be non-zero matrix.\n";
        cout<<"F * P * S - S * P * F = \n"<<(*el[i]->Fao_alp) * (*el[i]->P_alp) * (*el[i]->Sao) - (*el[i]->Sao) * (*el[i]->P_alp) * (*el[i]->Fao_alp)<<endl;
        cout<<"Density matrix constructed from eigenstates does commute with Fock matrix. This should be zero matrix.\n";
        cout<<"F * P * S - S * P * F = \n"<<(*el[i]->Fao_alp) * (*Ptmpa) * (*el[i]->Sao) - (*el[i]->Sao) * (*Ptmpa) * (*el[i]->Fao_alp)<<endl;       


        // Natural orbitals

        MATRIX* Id; Id = new MATRIX(el[i]->Norb,el[i]->Norb);
        Id->Init_Unit_Matrix(1.0);
        MATRIX* Nat_pop; Nat_pop = new MATRIX(el[i]->Norb,el[i]->Norb);
        MATRIX* Nat_orb; Nat_orb = new MATRIX(el[i]->Norb,el[i]->Norb);

        solve_eigen(el[i]->Norb, el[i]->P, Id, Nat_pop, Nat_orb);

        cout<<"Natural populations = \n"<<*Nat_pop<<endl;
        cout<<"Natural orbitals = \n"<<*Nat_orb<<endl;



        delete Ptmpa; delete Ptmpb;
        delete Id;
        delete Nat_pop; delete Nat_orb;

      }// for i - for all fragments
    }// if in DEBUG mode (DF==1)


  }// runtype == scf

  else if(prms.runtype=="nac"){

    cout<<"prms.runtype == nac\n";

    //----------------------------------------------------------------------------------
    //------------------------------ Read coordinates ----------------------------------
    std::string st;
    vector< vector<std::string> > file;
    cout<<"Reading input MD trajectory XYZ = "<<prms.nac_md_trajectory_filename<<endl;

    ifstream in(prms.nac_md_trajectory_filename.c_str(), ios::in);
    if(in.is_open()){
      while(!in.eof()){
        getline(in,st); 
     
        vector<std::string> line;
        stringstream ss(st,stringstream::in|stringstream::out);
        while(ss>>st){ line.push_back(st);}
     
        file.push_back(line);
     
      }// while
    }else{ cout<<"Error: Can not open file\n";  exit(0); }
    in.close();


    int Natoms = mol.Nnucl;
    cout<<"Natoms = "<<Natoms<<endl;


    //----------------------------------------------------------------------------------
    //----------------------- Now use only needed coordinates --------------------------
    vector<Electronic*> el_old; el_old = vector<Electronic*>(Nfrag);
    vector<AO> basis_ao_old; basis_ao_old = vector<AO>(basis_ao.size());

    // Prepare NAC matrices:
    vector<MATRIX*> nac; nac = vector<MATRIX*>(Nfrag*Nfrag);
    vector<MATRIX*> ene; ene = vector<MATRIX*>(Nfrag*Nfrag);
    vector<MATRIX*> ovlp; ovlp = vector<MATRIX*>(Nfrag*Nfrag);

    i = 0;
    int Xsize = 0;
    for(int f1=0;f1<Nfrag;f1++){     // fragment A
      int indx_min1 = prms.nac_min_orbs[f1];
      int indx_max1 = prms.nac_max_orbs[f1];

      for(int f2=0;f2<Nfrag;f2++){   // fragment B

        int indx_min2 = prms.nac_min_orbs[f2];
        int indx_max2 = prms.nac_max_orbs[f2];

        nac[i] = new MATRIX(indx_max1-indx_min1+1,indx_max2-indx_min2+1);
        ene[i] = new MATRIX(indx_max1-indx_min1+1,indx_max2-indx_min2+1);
        ovlp[i] = new MATRIX(indx_max1-indx_min1+1,indx_max2-indx_min2+1);
        i++;
      }// for f2
      Xsize += (indx_max1-indx_min1+1);
      
    }// for f1

    MATRIX* NAC; NAC = new MATRIX(Xsize,Xsize);
    MATRIX* ENE; ENE = new MATRIX(Xsize,Xsize);
    MATRIX* OVLP; OVLP = new MATRIX(Xsize,Xsize);
    MATRIX* invOVLP; invOVLP = new MATRIX(Xsize,Xsize);

    vector<int> frags(Xsize);
    vector<int> orbit(Xsize);



    

    for(int md_frame=prms.nac_min_frame;md_frame<=prms.nac_max_frame;md_frame++){


      cout<<"frame= "<<md_frame<<endl;

      // Actual coordinates and atom types
      for(i=0;i<Natoms;i++){ 
        VECTOR r; 
        r  = VECTOR(atof(file[(Natoms+2)*md_frame+i+2][1].c_str()),atof(file[(Natoms+2)*md_frame+i+2][2].c_str()),atof(file[(Natoms+2)*md_frame+i+2][3].c_str()));
        r *= Angst;  
        mol.R[i] = r;        
      }// for i

      // Update coordinates of basis functions
      for(i=0;i<basis_ao.size();i++){ basis_ao[i].set_position(mol.R[basis_ao[i].at_indx]); }
//      cout<<"basis_ao.size = "<<basis_ao.size()<<endl;



      // Electronic structure calculations for this geometry
      for(i=0;i<Nfrag;i++){
        guess(prms, modprms, mol, fragments[i], basis_fo[i], basis_ao, at_orbitals, el[i], mem);                // guess for this geometry

        //el0[i] = new Electronic(el[i]); 
        *el0[i] = *el[i];                                                                        // make copy at guess


        // Compute core Hamiltonian ...
        if(prms.hamiltonian=="indo" || prms.hamiltonian=="geht1"|| prms.hamiltonian=="geht2" ){

          if(prms.hamiltonian=="indo"){
            // Update core parameters and core Hamiltonian for this fragment
            indo_core_parameters(fragments[i], basis_fo[i], basis_ao, mol, mem->eri, mem->V_AB, mem, 1);
            Hamiltonian_core_indo(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals, el[i]->Hao, el[i]->Sao, mem, mem->eri, mem->V_AB);
          }
          else if(prms.hamiltonian=="geht1"){
            // Update core parameters and core Hamiltonian for this fragment
            //geht1_core_parameters(fragments[i], basis_fo[i], basis_ao, mol, mem->eri, mem->V_AB, mem, 1);
            Hamiltonian_core_geht1(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals, el[i]->Hao, el[i]->Sao, mem);
          }
          else if(prms.hamiltonian=="geht2"){
            // Update core parameters and core Hamiltonian for this fragment
            Hamiltonian_core_geht2(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals, el[i]->Hao, el[i]->Sao, mem);
          }

        }// if "indo" || "geht1" || "geht2"



        // Initialize Fock matrix
        Hamiltonian_Fock(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);

//      Hamiltonian_Fock_eht(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem); // init Fock


        // Compute SCF for all fragments
        if(prms.use_disk){

          //Eelec[i] = 
          scf_oda_disk(prms,modprms,mol,fragments[i],basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);

        }
        else{

          //Eelec[i] = 
          scf(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);

        }


//        scf(prms,modprms,mol,fragments[i], basis_fo[i],basis_ao,at_orbitals,el[i],el0[i],mem);                  // solve SCF problem

      }


      if(md_frame==prms.nac_min_frame){  // skip first point, but instantiate old copies
        for(i=0;i<Nfrag;i++){ el_old[i] = new Electronic(el[i]);  }
      }
      else{  // can do NAC computations
        cout<<"Computing nac for frame "<<md_frame<<endl;


        i = 0;
        for(int f1=0;f1<Nfrag;f1++){     // fragment A
          for(int f2=0;f2<Nfrag;f2++){   // fragment B

            if(DF){ cout<<"Coupling between fragments "<<f1<<" and "<<f2<<endl; }

            compute_overlap_nac(el_old[f1], el[f1], prms.nac_min_orbs[f1], prms.nac_max_orbs[f1], basis_fo[f1],
                                el_old[f2], el[f2], prms.nac_min_orbs[f2], prms.nac_max_orbs[f2], basis_fo[f2], 
                                basis_ao_old, basis_ao, prms.nac_dt, mem->aux, mem->n_aux, ovlp[i],nac[i],ene[i],prms.nac_opt);

            if(DF){
              cout<<"nac=\n"<<*nac[i]<<endl;
              cout<<"ovlp=\n"<<*ovlp[i]<<endl;
              cout<<"ene=\n"<<*ene[i]<<endl;
              cout<<"el[f1]->Sao = \n"<<*el[f1]->Sao<<endl;
            }

            i++;

          }// for f2
        }// for f1


        collect_matrices(Nfrag, ene, prms.nac_min_orbs, prms.nac_max_orbs, ENE, frags,orbit);
        collect_matrices(Nfrag, nac, prms.nac_min_orbs, prms.nac_max_orbs, NAC, frags,orbit);
        collect_matrices(Nfrag, ovlp, prms.nac_min_orbs, prms.nac_max_orbs, OVLP, frags,orbit);


        // Now output NAC between all fragments
        stringstream ss(stringstream::in | stringstream::out);
        std::string sframe;
        (ss << md_frame);  ss >> sframe;
        ofstream out_im((prms.nac_prefix+sframe+"_im").c_str(),ios::out);
        ofstream out_re((prms.nac_prefix+sframe+"_re").c_str(),ios::out);

        out_im << *NAC<<endl;
        out_re << *ENE<<endl; //0.5*(*el_old[0]->E_alp + *el[0]->E_alp)<<endl;
//        out<<"NAC=\n"<<*NAC<<"\n";
//        out<<"S=\n"<<*OVLP<<"\n";

        // Keep this for NA-MD in diabatic basis
        /*
        OVLP->Inverse(invOVLP);
        out<<"S^{-1} * NAC=\n"<<*invOVLP * *NAC<<"\n";
        */

        out_im.close();
        out_re.close();

      
      }// else - do NAC computations    


      // Make current values to be old
      for(i=0;i<basis_ao.size();i++){ 
        basis_ao_old[i] = basis_ao[i];
      }
      for(i=0;i<Nfrag;i++){ 
        *el_old[i]->E_alp = *el[i]->E_alp; //copy only this info, because we will need only it
        *el_old[i]->E_bet = *el[i]->E_bet; //copy only this info, because we will need only it

        *el_old[i]->C_alp = *el[i]->C_alp; //copy only this info, because we will need only it
        *el_old[i]->C_bet = *el[i]->C_bet; //copy only this info, because we will need only it
        el_old[i]->bands_alp = el[i]->bands_alp;
        el_old[i]->bands_bet = el[i]->bands_bet;
        el_old[i]->occ_alp = el[i]->occ_alp;
        el_old[i]->occ_bet = el[i]->occ_bet;


      }// i


    }// md_frame

  }// nac 



  

/*

//  int md_frame = 1; // for nac

  double en2,en1;

  if(prms.runtype=="info"){
    // Compute and print some volumetric and molecular info

    Mol.print_com(prms);
 
  }

  else if(prms.runtype=="foe"){  // Fermi operator expansion
    Mol.foe(prms,1);
  }

  else if(prms.runtype=="scf"){
    Mol.scf(prms,1);

    if(prms.guess_type=="fmo"){
    }
    else{
      Mol.print_Mulliken_pop(prms, filename+".Mull_pop");
      Mol.print_multipoles(prms, prms.spin_method, filename+".multipoles");
      Mol.print_el_struct(prms.spin_method, filename+".el_struct");
      if(prms.compute_excitations==1){
        Mol.print_excitations(prms,prms.spin_method, filename+".excitations");
      }

//      for(int i=0;i<prms.Natoms;i++){
//        Mol.R[i].x += 5.0;  // this is to test translational invariance
//      }
//      Mol.print_multipoles(prms.spin_method, filename+".multipoles+10");

    }

//    exit(0);
  }
  else if(prms.runtype=="scf_oda"){
    Mol.scf_oda(prms,1);

    if(prms.guess_type=="fmo"){
    }
    else{
      Mol.print_Mulliken_pop(prms, filename+".Mull_pop");
      Mol.print_multipoles1(prms, prms.spin_method, filename+".multipoles");
      Mol.print_el_struct(prms.spin_method, filename+".el_struct");
      if(prms.compute_excitations==1){
        Mol.print_excitations(prms,prms.spin_method, filename+".excitations");
      }

//      for(int i=0;i<prms.Natoms;i++){
//        Mol.R[i].x += 5.0;  // this is to test translational invariance
//      }
//      Mol.print_multipoles(prms.spin_method, filename+".multipoles+10");

    }

//    exit(0);
  }


  else if(prms.runtype=="scf_ls"){
    Mol.scf_ls(prms,1);
  }
  else if(prms.runtype=="scf_feeq"){
    Mol.scf_feeq(prms,1);
  }


  else if(prms.runtype=="dos"){
    Mol.scf_oda(prms,1);
    Mol.dos_dens(prms.dos_prefix,prms.hamiltonian,prms.spin_method,1);

    Mol.print_Mulliken_pop(prms, filename+".Mull_pop");
    Mol.print_multipoles1(prms, prms.spin_method, filename+".multipoles");
    Mol.print_el_struct(prms.spin_method, filename+".el_struct");

    if(prms.compute_excitations==1){
      Mol.print_excitations(prms,prms.spin_method, filename+".excitations");
    }

//    exit(0);
  }

  else if(prms.runtype=="charge_density"){
    Mol.scf(prms,1);
    Mol.charge_density(prms,1);
//    exit(0);
  }



  else if(prms.runtype=="nac"){

    int turn = 0;
    Mol.scf(prms,1);   

//exit(0);

    MATRIX* Sao;  Sao = new MATRIX(Mol.Norb,Mol.Norb);
    
    // Initialize general-purpose auxiliary arrays
    int n_aux = 20;
    vector<double*> aux;
    for(int i=0;i<40;i++){    // 15 auxiliary arrays with double
      double* x; x = new double[n_aux];
      aux.push_back(x);
    }


    // Read in the input file:
    std::string st;
    vector< vector<std::string> > file;
      
    //---------------------------- Reading input file --------------------------------------
   
    cout<<"Reading input MD trajectory XYZ = "<<prms.nac_md_trajectory_filename<<endl;
    ifstream in(prms.nac_md_trajectory_filename.c_str(), ios::in);
    if(in.is_open()){
      while(!in.eof()){
        getline(in,st); 
   
        vector<std::string> line;
        stringstream ss(st,stringstream::in|stringstream::out);
        while(ss>>st){ line.push_back(st);}
   
        file.push_back(line);
   
      }// while
    }else{ cout<<"Error: Can not open file\n";  exit(0); }
    in.close();


    int Natoms = Mol.Nnucl;
    cout<<"Natoms = "<<Natoms<<endl;

//exit(0);

    for(int md_frame=prms.nac_min_frame;md_frame<=prms.nac_max_frame;md_frame++){

      cout<<"frame= "<<md_frame<<endl;

     //!!!!!!!!!!!!!!!!! THIS NEEDS REVISION !!!!!!!!!!!!!!

      // Actual coordinates and atom types
      for(int i=0;i<Natoms;i++){ 
        VECTOR r; 
        r  = VECTOR(atof(file[(Natoms+2)*md_frame+i+2][1].c_str()),atof(file[(Natoms+2)*md_frame+i+2][2].c_str()),atof(file[(Natoms+2)*md_frame+i+2][3].c_str()));
        cout<<i<<"  "<<r<<endl;
        r *= Angst;                


 
        // Convert coordinates
//        if(coordinates=="Cartesian"){ r *= Angst; } // convert to Bohr
//        else if(coordinates=="Crystal"){  // so far we assume it is cubic cell!!!
//          cout<<"Warning: This works only for cubic cell so far. Check that this is what you have.\n";
//          r *= cell.a;
//        }

        if(turn==0){  Mol1.R[i] = r; }
        else if(turn==1){ Mol2.R[i] = r; }
        else if(turn==2){ Mol.R[i] = r; }


      }// for i
      cout<<endl;

      cout<<"Computing nac for frame "<<md_frame<<endl;

      int do_guess = 1;
      if(md_frame-prms.nac_min_frame>1){  do_guess = 0; }  // to save alot!

//exit(0);

      if(turn==0){
        Mol1.scf(prms,do_guess);    
//        compute_nac(Mol,Mol1,dt,nac_prefix,md_frame); 
        compute_nac(Mol,Mol1,prms.nac_dt,prms.nac_prefix,md_frame, prms.nac_min_orb,prms.nac_max_orb,Sao,aux,n_aux); 

      }
      else if(turn==1){
        Mol2.scf(prms,do_guess);    
//        compute_nac(Mol1,Mol2,dt,nac_prefix,md_frame);        
        compute_nac(Mol1,Mol2,prms.nac_dt,prms.nac_prefix,md_frame, prms.nac_min_orb,prms.nac_max_orb,Sao,aux,n_aux); 
      }
      else if(turn==2){
        Mol.scf(prms,do_guess);    
//        compute_nac(Mol2,Mol,dt,nac_prefix,md_frame);        
        compute_nac(Mol2,Mol,prms.nac_dt,prms.nac_prefix,md_frame, prms.nac_min_orb,prms.nac_max_orb,Sao,aux,n_aux); 
      }

      turn++;
      turn = turn%2;

    }// md_frame

    delete Sao;

  }// nac

*/

/*
  else if(runtype=="scan"){


    ofstream out1("scan_energ.txt",ios::out);
    ofstream out2("scan_eig.txt",ios::out);
    ofstream out2a("scan_eig_alp.txt",ios::out);
    ofstream out2b("scan_eig_bet.txt",ios::out);
    ofstream out3("scan_coord.xyz",ios::out);

    for(double x=scan_dxmin;x<scan_dxmax;x+= scan_dx){

      Mol.R[scan_mov_at] += x*scan_dir;

      double dx = 0.001;
      Mol.R[scan_mov_at] -= dx*scan_dir;
      Mol.scf(method,spin_method,1);
      *e1 = *Mol.E;
      *f1 = *Mol.Fao;
      *c1 = *Mol.C;
      en1 = Mol.Etot;

      Mol.R[scan_mov_at] += 2.0*dx*scan_dir;
      Mol.scf(method,spin_method,1);
      *e2 = *Mol.E;
      *f2 = *Mol.Fao;
      *c2 = *Mol.C;
      en2 = Mol.Etot;
  
      Mol.R[scan_mov_at] -= dx*scan_dir;
      Mol.scf(method,spin_method,1);



      out1<<"R(a.u.)= "<<(Mol.R[scan_mov_at]-Mol.R[scan_ref_at]*scan_dir)*scan_dir<<" E(a.u.)= "<<Mol.Etot
          <<" grad(anal)= "<<(Mol.grad[scan_mov_at]*scan_dir)*scan_dir<<"  grad(numer)= "<<((en2-en1)/(2.0*dx))*scan_dir
          <<"Enucl+Epair(a.u.)= "<<Mol.Enucl+Mol.Epair<<endl;  


      if(spin_method=="restricted"){

        out2<<"R(a.u.)= "<<(Mol.R[scan_mov_at]-Mol.R[scan_ref_at]*scan_dir)*scan_dir;
        for(int i=0;i<Mol.Norb;i++){ out2<<" E["<<i<<"],eV= "<<Mol.E->M[i*Mol.Norb+i]/eV<<" ";   }
        out2<<endl;

      }
      else if(spin_method=="unrestricted"){

        out2a<<"R(a.u.)= "<<(Mol.R[scan_mov_at]-Mol.R[scan_ref_at]*scan_dir)*scan_dir;
//        for(int i=0;i<Mol.Norb;i++){ out2a<<" E["<<i<<"],eV= "<<Mol.E_alp->M[i*Mol.Norb+i]/eV<<" ";   }
        for(int i=0;i<Mol.Norb;i++){ out2a<<" E["<<i<<"],eV= "<<Mol.bands_alp[i].second/eV<<" ";   }
        out2a<<endl;

        out2b<<"R(a.u.)= "<<(Mol.R[scan_mov_at]-Mol.R[scan_ref_at]*scan_dir)*scan_dir;
//        for(int i=0;i<Mol.Norb;i++){ out2b<<" E["<<i<<"],eV= "<<Mol.E_bet->M[i*Mol.Norb+i]/eV<<" ";   }
        for(int i=0;i<Mol.Norb;i++){ out2b<<" E["<<i<<"],eV= "<<Mol.bands_bet[i].second/eV<<" ";   }
        out2b<<endl;


      }


 

//      cout<<"E'(numer)= \n"<<(*e2 - *e1)/(2.0*dx)<<"\n E'(anal)= \n"; for(int i=0;i<Mol.Norb;i++){   cout<<Mol.dEdR[i][0]<<endl;  }  cout<<endl;
//      cout<<"Fao'(numer)= \n"<<(*f2 - *f1)/(2.0*dx)<<"\n Fao'(anal)= \n"; cout<<*Mol.dFao_dx[0]<<endl;
//      cout<<"C'(numer)= \n"<<(*c2 - *c1)/(2.0*dx)<<"\n C'(anal)= \n"; 
//      cout<<"en1 = "<<en1<<endl;
//      cout<<"en2 = "<<en2<<endl;
//      cout<<"Eelec'(numer) = "<<(en2-en1)/(2.0*dx)<<"  Eelec'(anal)= "<<Mol.grad[0]<<endl;


      out3<<Mol.Nnucl<<endl;
      out3<<"molecule"<<endl;
      for(int n=0;n<Mol.Nnucl;n++){
        out3<<Mol.at_types[n]<<"   "<<(1.0/Angst)*Mol.R[n]<<endl;
      }


      Mol.R[scan_mov_at] -= x*scan_dir;

    }// for x

    out1.close();
    out2.close();
    out2a.close();
    out2b.close();
    out3.close();

  }// scan
  else if(runtype=="test1"){


    double dx = 0.001;
    Mol.scf(method,spin_method,1);

    Mol1.R[1].x -= dx;
    Mol1.scf(method,spin_method,1);

    Mol2.R[1].x += dx;
    Mol2.scf(method,spin_method,1);

     // Initialize general-purpose auxiliary arrays
  vector<double*> aux;
  int n_aux = 20;
  for(int i=0;i<40;i++){    // 15 auxiliary arrays with double
    double* x; x = new double[n_aux];
    aux.push_back(x);
  }

  vector<VECTOR*> auxv;
  int n_auxv = 20;
  for(int i=0;i<40;i++){    // 5 auxiliary arrays with VECTOR
    VECTOR* x; x = new VECTOR[n_auxv];
    auxv.push_back(x);
  }

  cout<<"Numerical couplings\n";

  MATRIX* nac2;

  int Norb = Mol.Norb;

  nac2 =  new MATRIX(Norb,Norb);
  *nac2 = 0.0;

    for(int i=0;i<Norb;i++){
      for(int j=0;j<Norb;j++){
        VECTOR dIdA, dIdB;

        double ovlp1, ovlp2, nac, nac1,der1,der2;
        ovlp1 = 0.0;
        ovlp2 = 0.0;
        nac = 0.0;
        nac1 = 0.0;

        for(int k1=0;k1<Norb;k1++){
          for(int k2=0;k2<Norb;k2++){

            //if(k1==k2){
            ovlp1 += Mol1.C->M[k1*Norb+i]*Mol2.C->M[k2*Norb+j]*OVERLAP_INTEGRAL(Mol1.basis_ao[k1],Mol2.basis_ao[k2],0,dIdA,dIdB,aux,n_aux);
            ovlp2 += Mol2.C->M[k2*Norb+i]*Mol1.C->M[k1*Norb+j]*OVERLAP_INTEGRAL(Mol2.basis_ao[k2],Mol1.basis_ao[k1],0,dIdA,dIdB,aux,n_aux);        
            

            der1 = DERIVATIVE_COUPLING_INTEGRAL(Mol.basis_ao[k1],Mol.basis_ao[k2],aux,n_aux).x;
            der2 = DERIVATIVE_COUPLING_INTEGRAL(Mol.basis_ao[k2],Mol.basis_ao[k1],aux,n_aux).x;

//            cout<<"der1= "<<der1<<" der2= "<<der2<<endl;

            nac += Mol.C->M[k1*Norb+i]*Mol.C->M[k2*Norb+j]*Mol.dFao_dx[1]->M[k1*Norb+k2];

           // nac2 += Mol.C->M[k1*Norb+i]*Mol.C->M[k2*Norb+j]*(der1-der2);

            nac1 += Mol.C->M[k1*Norb+i]*Mol.C->M[k2*Norb+j]*(Mol.dFao_dx[1]->M[k1*Norb+k2] - (Mol.E->M[i*Norb+i]*der1 + Mol.E->M[j*Norb+j]*der2) );

          }// for k2
        }// for k1

        nac /= (Mol.E->M[j*Norb+j]-Mol.E->M[i*Norb+i]);
        nac1 /= (Mol.E->M[j*Norb+j]-Mol.E->M[i*Norb+i]);

        nac2->M[i*Norb+j] = (ovlp1 - ovlp2)/(4.0*dx);

        cout<<"numer,coupling("<<i<<","<<j<<")= "<<(ovlp1 - ovlp2)/(4.0*dx)<<endl;
        cout<<"anal, coupling("<<i<<","<<j<<")= "<<nac<<endl;
        cout<<"anal1, coupling("<<i<<","<<j<<")= "<<nac1<<endl;


      }
    }  


    cout<<"nac2= \n"<<*nac2<<endl;
    Mol.test1(method);



  }// test1

  else if(runtype=="md"|| runtype=="opt"){

    ofstream out1("md.txt",ios::out);
    ofstream out2("md.xyz",ios::out);

    double Ekin = 0.0;
    double Epot = 0.0;
    double Etot = Epot + Ekin;

    // Set random velocities:
    VECTOR totP; totP = 0.0;
    for(int n=0;n<Mol.Nnucl;n++){
      Mol.P_[n] = VECTOR(normal(),normal(),normal());
      totP += Mol.P_[n];
    }
    // Remove COM motion:
    totP = totP/((double)Mol.Nnucl);

    // Also need to remove total system rotation
   
    for(int n=0;n<Mol.Nnucl;n++){
      Mol.P_[n] -= totP;
      Ekin += 0.5*Mol.P_[n].length2()/Mol.mass[n];
    }
    
    double T = 300.0; // K
    double kb = 3.1668114e-6;  // Ha/K

    double Ekin_target = 0.5*(3.0*Mol.Nnucl-6.0)*kb*T; 
    double scl = sqrt(Ekin_target/Ekin);

    // Rescale velocities:
    for(int n=0;n<Mol.Nnucl;n++){
      Mol.P_[n] *= scl;
      Ekin += 0.5*Mol.P_[n].length2()/Mol.mass[n];
    }


     

    Mol.scf(method,spin_method,1);

    for(int t=0;t<nsteps;t++){
      
      for(int n=0;n<Mol.Nnucl;n++){
        Mol.P_[n] -= 0.5*dt*Mol.grad[n];
        Mol.R[n] += dt*Mol.P_[n]/Mol.mass[n];
      }

      Mol.scf(method,spin_method,0);

      Epot = Mol.Etot;


      Ekin = 0.0;
      for(int n=0;n<Mol.Nnucl;n++){
        Mol.P_[n] -= 0.5*dt*Mol.grad[n];
        Ekin += 0.5*Mol.P_[n].length2()/Mol.mass[n];

        if(runtype=="opt"){
          Mol.P_[n] = 0.0;
        }
      }

      Etot = Epot + Ekin;
      
      double Tcurr = 2.0*(Ekin/((3.0*Mol.Nnucl-6.0)*kb)); 

      out1<<"t= "<<t<<" Etot= "<<Etot<<" Epot= "<<Epot<<" Ekin= "<<Ekin<<" Tcurr= "<<Tcurr<<endl; //(Mol.R[0] - Mol.R[1]).length()<<endl;

      out2<<Mol.Nnucl<<endl;
      out2<<"molecule"<<endl;
      for(int n=0;n<Mol.Nnucl;n++){
        out2<<Mol.at_types[n]<<"   "<<(1.0/Angst)*Mol.R[n]<<endl;
      }
    
    
    }// for t

    out1.close();
    out2.close();

  }// md

  

//    Mol.compute_forces(method);

*/
  ti2 = clock();
  cout<<"Time elapsed: "<<(ti2-ti1)/((double)CLOCKS_PER_SEC)<<" sec"<<endl;


  for(i=0;i<Nfrag;i++){
    el[i]->~Electronic();
    el0[i]->~Electronic();
  }




  return 0;
}
