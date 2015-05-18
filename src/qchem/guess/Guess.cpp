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
 \file Guess.cpp
 \brief Implementation of the methods for guessing electronic structure of fragment

*/
/****************************************************************************
  This file contains following functions:

  void guess(Control_Parameters& prms, Model_Parameters& modprms, Nuclear& mol,
             vector<int>& fragment, vector<int>& basis_fo, vector<AO>& basis_ao, vector<vector<int> >& at_orbitals, 
             Electronic* el,Memory* mem)

****************************************************************************/
#include "Guess.h"

void guess(Control_Parameters& prms, Model_Parameters& modprms, Nuclear& mol,
           vector<int>& fragment, vector<int>& basis_fo, vector<AO>& basis_ao, vector<vector<int> >& at_orbitals, 
           Electronic* el,Memory* mem){

  int i,a,b;
  std::string eigen_method="generalized";

  vector<Timer> bench_t2(4);
  Timer t;

  //----------------------------------------------------------------------------------
  //-------------------- Compute core Hamiltonian for given fragment -----------------
  // Prepare matrices    
  cout<<"Nelec= "<<el->Nelec<<endl;
  cout<<"Nocc_alp = "<<el->Nocc_alp<<endl;
  cout<<"Nocc_bet = "<<el->Nocc_bet<<endl;

  t.start();
  update_overlap_matrix(prms.x_period,prms.y_period,prms.z_period,prms.t1,prms.t2,prms.t3,
                        basis_fo,basis_ao,el->Sao, mem->aux, mem->n_aux, mol);
//  cout<<"Sao = "<<*el->Sao<<endl;
  t.stop();
  cout<<"Computation of overlap matrix takes: "<<t.show()<<" seconds\n";


//  t.start();
//  update_overlap_matrix_new(prms.x_period,prms.y_period,prms.z_period,prms.t1,prms.t2,prms.t3,
//                            basis_fo,basis_ao,el->Sao, mem->aux, mem->n_aux, mol);
//  cout<<"Sao = "<<*el->Sao<<endl;
//  t.stop();
//  cout<<"Computation of overlap matrix takes: "<<t.show()<<" seconds\n";
//  exit(0);

  
  if(prms.hamiltonian=="indo"){    

    el->Sao->Init_Unit_Matrix(1.0);  

    // 1 - for INDO
    // 0 - for CNDO
    indo_core_parameters(fragment, basis_fo, basis_ao, mol, mem->eri, mem->V_AB, mem, 1);

  }// if "indo"
  if(prms.hamiltonian=="geht1"){    

    //el->Sao->Init_Unit_Matrix(1.0);  

    // 1 - for INDO
    // 0 - for CNDO
    //geht1_core_parameters(fragment, basis_fo, basis_ao, mol, mem->eri, mem->V_AB, mem, 1);

  }// if "geht1"



  t.start();      
  Hamiltonian_core(prms,modprms,mol,fragment,basis_fo,basis_ao,at_orbitals,el->Hao,el->Sao,mem);

  t.stop();
  cout<<"Computation of core Hamiltonian matrix takes: "<<t.show()<<" seconds\n";



  *el->Fao_alp = *el->Hao;
  *el->Fao_bet = *el->Hao;


  Fock_to_P(el->Norb,el->Nocc_alp, 1, el->Nelec, eigen_method, prms.pop_opt,
            el->Fao_alp, el->Sao, el->C_alp, el->E_alp, el->bands_alp, el->occ_alp, el->P_alp, bench_t2);

  Fock_to_P(el->Norb,el->Nocc_bet, 1, el->Nelec, eigen_method, prms.pop_opt,
            el->Fao_bet, el->Sao, el->C_bet, el->E_bet, el->bands_bet, el->occ_bet, el->P_bet, bench_t2);

  *el->P = *el->P_alp + *el->P_bet;


//  el0 = new Electronic(el); 
//  cout<<"Composed density matrix for fragment f = "<<i<<endl<<*el->P<<endl;


  //----------------------------------------------------------------------------------
  //--------------------------- Guess fragment MOs -----------------------------------

  if(prms.guess_type=="core"){
    // Nothing special goes here - the "core" density matrix is already computed in the section above
  }
  else if(prms.guess_type=="sad"){ //
    t.start();
    // Initialize
    vector< vector<int> > fragments_sad;  
    vector< vector<int> > basis_fo_sad;
    vector<double> charges(mol.Nnucl,0.0);
    Memory* mem0; mem0 = new Memory(20,40, 20,40); // auxiliary variables for temporary memory storage


    for(i=0;i<mol.Nnucl;i++){ fragments_sad.push_back(vector<int>(1,i));  }  

    vector<Electronic*> el_sad;
    int Nfrag_sad = init_electronic_subsystems(fragments_sad, at_orbitals, basis_fo_sad, el_sad, basis_ao, modprms,mol,charges);  

    cout<<"Fragment:  {Orbitals}\n";
    show_mapping(basis_fo_sad);

  

    // Prepare matrices    
    for(i=0;i<Nfrag_sad;i++){
 
      cout<<"Nelec= "<<el_sad[i]->Nelec<<endl;
      cout<<"Nocc_alp = "<<el_sad[i]->Nocc_alp<<endl;
      cout<<"Nocc_bet = "<<el_sad[i]->Nocc_bet<<endl;



      if(prms.hamiltonian=="indo"){
        indo_core_parameters(fragments_sad[i], basis_fo_sad[i], basis_ao, mol, mem0->eri, mem0->V_AB, mem0, 1);
      }// if "indo"
      else if(prms.hamiltonian=="indo"){
        geht1_core_parameters(fragments_sad[i], basis_fo_sad[i], basis_ao, mol, mem0->eri, mem0->V_AB, mem0, 1);
      }// if "indo"



      update_overlap_matrix(prms.x_period,prms.y_period,prms.z_period,prms.t1,prms.t2,prms.t3,
                            basis_fo_sad[i],basis_ao,el_sad[i]->Sao, mem->aux, mem->n_aux, mol);

      
      Hamiltonian_core(prms,modprms,mol,fragments_sad[i],basis_fo_sad[i],basis_ao,at_orbitals,el_sad[i]->Hao,el_sad[i]->Sao,mem0);
      
      *el_sad[i]->Fao_alp = *el_sad[i]->Hao;
      *el_sad[i]->Fao_bet = *el_sad[i]->Hao;
  

      Fock_to_P(el_sad[i]->Norb,el_sad[i]->Nocc_alp, 1, el_sad[i]->Nelec, eigen_method, prms.pop_opt,
                el_sad[i]->Fao_alp, el_sad[i]->Sao, el_sad[i]->C_alp, el_sad[i]->E_alp,
                el_sad[i]->bands_alp, el_sad[i]->occ_alp, el_sad[i]->P_alp, bench_t2);

  
      Fock_to_P(el_sad[i]->Norb,el_sad[i]->Nocc_bet, 1, el_sad[i]->Nelec, eigen_method, prms.pop_opt,
                el_sad[i]->Fao_bet, el_sad[i]->Sao, el_sad[i]->C_bet, el_sad[i]->E_bet,
                el_sad[i]->bands_bet, el_sad[i]->occ_bet, el_sad[i]->P_bet, bench_t2);

      *el_sad[i]->P = *el_sad[i]->P_alp + *el_sad[i]->P_bet;
      
//      cout<<"P = \n"<<*el_sad[i]->P<<endl;
//      cout<<"Fao_alp = \n"<<*el_sad[i]->Fao_alp<<endl;
//      cout<<"Sao = \n"<<*el_sad[i]->Sao<<endl;

    }// for i


    // Now we can compose the fragment density matrices from the atomics ones
    *el->P_alp = 0.0;
    *el->P_bet = 0.0;
    *el->P     = 0.0;

    int N1 = basis_fo.size(); // how many orbitals is fragment f1

    for(int i1=0;i1<N1;i1++){
      for(int j1=0;j1<N1;j1++){ 

        for(int f2=0;f2<basis_fo_sad.size();f2++){  

          int N2 = basis_fo_sad[f2].size(); // how many orbitals is fragment f2
     
          for(int i2=0;i2<N2;i2++){ 
            for(int j2=0;j2<N2;j2++){ 
 
              if( (basis_fo[i1]==basis_fo_sad[f2][i2]) && (basis_fo[j1]==basis_fo_sad[f2][j2])){

                el->P_alp->M[i1*N1+j1] = el_sad[f2]->P_alp->M[i2*N2+j2];
                el->P_bet->M[i1*N1+j1] = el_sad[f2]->P_bet->M[i2*N2+j2];
                el->P->M[i1*N1+j1]     = el_sad[f2]->P->M[i2*N2+j2];

              }// same orbital

            }// for all orbitals in the fragment f2 - j2
          }// for all orbitals in the fragment f2 - i2
        }// for all sad fragment - f2

      }// for all orbitals in the fragment f1 - j1
    }// for all orbitals in the fragment f1 - i1

//    cout<<"Composed density matrix for given fragment = "<<endl<<*el->P<<endl;

    // Delete temporary SAD objects
    for(i=0;i<Nfrag_sad;i++){
      el_sad[i]->~Electronic();
    }

    t.stop();
    cout<<"SAD computations take "<<t.show()<<" seconds\n";


//    mem0->~Memory();

  }// guess_type == "sad"

//  cout<<"End with guess\n";

}// guess
