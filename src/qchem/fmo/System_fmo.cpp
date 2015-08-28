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

#include "System.h"
#include "Density_Builder.h"


void System::FMO_init_subsets(vector< vector<int> >& fragment){
// fragment[i] - is a set of atoms that belong to the i-th fragment

  int i,j,k,a,b,n;

  //----------- Initialize orbital indices subsets ------------

  if(FMO_subsets.size()>0){ FMO_subsets.clear(); }

  FMO_num_subsets = fragment.size(); // not a fancy way, but it'll work for now


  // Allocate memory for the number of orbitals in each fragment (subset)
  // and the number of occupied (alp, bet, and total) orbitals in each subset ----- 
  FMO_norb = vector<int>(FMO_num_subsets, 0);
  FMO_nocc = vector<int>(FMO_num_subsets, 0);
  FMO_nocc_alp = vector<int>(FMO_num_subsets, 0);
  FMO_nocc_bet = vector<int>(FMO_num_subsets, 0);

  // Because the following objects are only temporary, we don't create alp,bet,total variants, only a generic one
  FMO_Sao = new MATRIX*[FMO_num_subsets];
  FMO_Fao = new MATRIX*[FMO_num_subsets];
  FMO_Fao_alp = new MATRIX*[FMO_num_subsets];
  FMO_Fao_bet = new MATRIX*[FMO_num_subsets];
  FMO_P = new MATRIX*[FMO_num_subsets];
  FMO_P_alp = new MATRIX*[FMO_num_subsets];
  FMO_P_bet = new MATRIX*[FMO_num_subsets];
  FMO_C = new MATRIX*[FMO_num_subsets];
  FMO_C_alp = new MATRIX*[FMO_num_subsets];
  FMO_C_bet = new MATRIX*[FMO_num_subsets];
  FMO_E = new MATRIX*[FMO_num_subsets];
  FMO_E_alp = new MATRIX*[FMO_num_subsets];
  FMO_E_bet = new MATRIX*[FMO_num_subsets];

  vector< pair<int,double> > x; // empty
  FMO_bands = vector< vector< pair<int,double> > >(FMO_num_subsets, x);
  FMO_bands_alp = vector< vector< pair<int,double> > >(FMO_num_subsets, x);
  FMO_bands_bet = vector< vector< pair<int,double> > >(FMO_num_subsets, x);
  FMO_occ = vector< vector< pair<int,double> > >(FMO_num_subsets, x);
  FMO_occ_alp = vector< vector< pair<int,double> > >(FMO_num_subsets, x);
  FMO_occ_bet = vector< vector< pair<int,double> > >(FMO_num_subsets, x);

 
  for(i=0;i<FMO_num_subsets;i++){
    vector<int> subset;
    int nelec = 0;

    for(j=0;j<fragment[i].size();j++){
      n = fragment[i][j];

      for(k=0;k<at_orbitals[n].size();k++){   subset.push_back( at_orbitals[n][k] );  }// k - all orbitals centered on atom n 

      nelec += PT[at_types[n]].Nval;

    }// j (=>n) - all atoms in the fragment i

    FMO_subsets.push_back(subset);

    FMO_norb[i] = subset.size();
    FMO_nocc[i] = (nelec%2==0)?nelec/2:(nelec/2+1); 
    FMO_nocc_alp[i] = FMO_nocc[i];
    FMO_nocc_bet[i] = nelec - FMO_nocc_alp[i];


    for(j=0;j<FMO_norb[i];j++){
      FMO_bands[i].push_back(pair<int,double>(j,0.0));
      FMO_bands_alp[i].push_back(pair<int,double>(j,0.0));
      FMO_bands_bet[i].push_back(pair<int,double>(j,0.0));

      
      FMO_occ[i].push_back(pair<int,double>(j, ((j<FMO_nocc[i])?2.0:0.0) ) );
      FMO_occ_alp[i].push_back(pair<int,double>(j, ((j<FMO_nocc_alp[i])?1.0:0.0) ) );
      FMO_occ_bet[i].push_back(pair<int,double>(j, ((j<FMO_nocc_bet[i])?1.0:0.0) ) );
    }// for j


    FMO_Sao[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_Fao[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_Fao_alp[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_Fao_bet[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_P[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_P_alp[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_P_bet[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_C[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_C_alp[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_C_bet[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_E[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_E_alp[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);
    FMO_E_bet[i] = new MATRIX(FMO_norb[i],FMO_norb[i]);

  }// i - all fragments


}// void System::FMO_init_subsets(vector< vector<int> >& fragment)


void System::FMO_sad(Control_Parameters& prms,std::string eigen_method){
// Long set of instructions and dealing with temporary objects is hidden in this function,
// so one is only focused on preparation of proper Hamiltonian and on resulting atomic densities,
// when this function is called

  // Split generated Hamiltoian into Nnucl tasks - one per each atom
  vector< vector<int> > fragment;
  for(int n=0;n<Nnucl;n++){
    vector<int> fr(1,n);
    fragment.push_back(fr);
  }

  // Allocate memory and init some basis numbers
  FMO_init_subsets(fragment);


  if(prms.spin_method=="unrestricted"){
    *P_alp = 0.0;
    *P_bet = 0.0;

    //---- Handle subsets ------------
    // Alpha
    handle_n_subsets(Fao_alp, Sao, P_alp, fragment.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc_alp, 1, eigen_method, prms.pop_opt);

    handle_n_subsets(Fao_bet, Sao, P_bet, fragment.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc_bet, 1, eigen_method, prms.pop_opt);

    
    *P = *P_alp + *P_bet;

  }// unrestricted

  else if(prms.spin_method=="restricted"){
    *P = 0.0;

    //---- Handle subsets ------------
    // Alpha
    handle_n_subsets(Fao, Sao, P, fragment.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc, 2, eigen_method, prms.pop_opt);

  }// restricted



  if(0){
    cout<<"Total composed density matrix\n";
    cout<<"P_alp=\n"<<*P_alp<<endl;
    cout<<"P_bet=\n"<<*P_bet<<endl;
  }

}//void System::FMO_sad(Control_Parameters& prms,std::string eigen_method)


void System::FMO_gen(Control_Parameters& prms,std::string eigen_method){
// Generic handling of FMOs
  int k;

  // Allocate memory and init some basis numbers
  FMO_init_subsets(prms.fragments);


  if(prms.spin_method=="unrestricted"){
    *P_alp = 0.0;
    *P_bet = 0.0;

      for(k=0;k<FMO_num_subsets;k++){

        MATRIX* s; s = new MATRIX(FMO_norb[k],FMO_norb[k]);
        pop_submatrix(Fao_alp, FMO_Fao_alp[k], FMO_subsets[k]);
        pop_submatrix(Fao_bet, FMO_Fao_bet[k], FMO_subsets[k]);
        pop_submatrix(Sao, s, FMO_subsets[k]);

        Fock_to_P(FMO_norb[k],FMO_nocc_alp[k], 1, FMO_nocc_alp[k], eigen_method, prms.pop_opt,
                 FMO_Fao_alp[k], s, FMO_C_alp[k], FMO_E_alp[k], FMO_bands_alp[k], FMO_occ_alp[k], FMO_P_alp[k]);


        Fock_to_P(FMO_norb[k],FMO_nocc_bet[k], 1, FMO_nocc_bet[k], eigen_method, prms.pop_opt,
                 FMO_Fao_bet[k], s, FMO_C_bet[k], FMO_E_bet[k], FMO_bands_bet[k], FMO_occ_bet[k], FMO_P_bet[k]);

        delete s;
        
        cout<<"Fragment k="<<k<<"\n"<<*FMO_E_alp[k]<<endl;

      }// for k


/*
    //---- Handle subsets ------------
    // Alpha
    handle_n_subsets(Fao_alp, Sao, P_alp, prms.fragments.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc_alp, 1, eigen_method, prms.pop_opt);

    handle_n_subsets(Fao_bet, Sao, P_bet, prms.fragments.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc_bet, 1, eigen_method, prms.pop_opt);

    
    *P = *P_alp + *P_bet;
*/

  }// unrestricted

  else if(prms.spin_method=="restricted"){
    *P = 0.0;
/*

    //---- Handle subsets ------------
    // Alpha
    handle_n_subsets(Fao, Sao, P, prms.fragments.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc, 2, eigen_method, prms.pop_opt);
*/

  }// restricted



  if(0){
    cout<<"Total composed density matrix\n";
    cout<<"P_alp=\n"<<*P_alp<<endl;
    cout<<"P_bet=\n"<<*P_bet<<endl;
  }

}//void System::FMO_gen(Control_Parameters& prms,std::string eigen_method)



void System::FMO_gen1(Control_Parameters& prms,std::string eigen_method){
// Generic handling of FMOs
  int k;

  // Allocate memory and init some basis numbers
  FMO_init_subsets(prms.fragments);


  if(prms.spin_method=="unrestricted"){
    *P_alp = 0.0;
    *P_bet = 0.0;

    vector<double> N_el_alp(FMO_num_subsets,0.0); // Number of electrons in each fragment
    vector<double> N_el_bet(FMO_num_subsets,0.0); // Number of electrons in each fragment
    
    for(k=0;k<FMO_num_subsets;k++){
      N_el_alp[k] = FMO_nocc_alp[k];
      N_el_bet[k] = FMO_nocc_bet[k];

    }// k

    vector<double> ef_alp(FMO_num_subsets,0.0);
    vector<double> ef_bet(FMO_num_subsets,0.0);

      int stop = 0;

      while(!stop){

        double err = 0.0;
        double ef_alp0;
        double ef_bet0;
        

        // Solve eigenvalue problems
        for(k=0;k<FMO_num_subsets;k++){
     
    
          MATRIX* s; s = new MATRIX(FMO_norb[k],FMO_norb[k]);
          pop_submatrix(Fao_alp, FMO_Fao_alp[k], FMO_subsets[k]);
          pop_submatrix(Fao_bet, FMO_Fao_bet[k], FMO_subsets[k]);
          pop_submatrix(Sao, s, FMO_subsets[k]);
    
    
          Fock_to_P(FMO_norb[k],FMO_nocc_alp[k], 1, N_el_alp[k], eigen_method, prms.pop_opt,
                   FMO_Fao_alp[k], s, FMO_C_alp[k], FMO_E_alp[k], FMO_bands_alp[k], FMO_occ_alp[k], FMO_P_alp[k]);

          ef_alp[k] = fermi_energy(FMO_bands_alp[k],N_el_alp[k],1.0);


//!!!          p_alp =   population(FMO_bands_alp[k][i], ef_alp[k], 1.0);
    
    
          Fock_to_P(FMO_norb[k],FMO_nocc_bet[k], 1, N_el_bet[k], eigen_method, prms.pop_opt,
                   FMO_Fao_bet[k], s, FMO_C_bet[k], FMO_E_bet[k], FMO_bands_bet[k], FMO_occ_bet[k], FMO_P_bet[k]);

          ef_bet[k] = fermi_energy(FMO_bands_bet[k],N_el_bet[k],1.0);

    
          delete s;
          
//          cout<<"Fragment k="<<k<<"  E_f(alp)= "<<ef_alp<<"  E_f(bet)= "<<ef_bet<<endl;
//          cout<<"Fragment k="<<k<<"\n"<<*FMO_E_alp[k]<<endl;

//          err += (ef_alp - ef_alp0)*(ef_alp - ef_alp0) + (ef_bet - ef_bet0)*(ef_bet - ef_bet0);   
          
    
        }// for k


        // Update


        if(err<1e-5){ stop = 1; }


      }// while(!stop)


/*
    //---- Handle subsets ------------
    // Alpha
    handle_n_subsets(Fao_alp, Sao, P_alp, prms.fragments.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc_alp, 1, eigen_method, prms.pop_opt);

    handle_n_subsets(Fao_bet, Sao, P_bet, prms.fragments.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc_bet, 1, eigen_method, prms.pop_opt);

    
    *P = *P_alp + *P_bet;
*/

  }// unrestricted

  else if(prms.spin_method=="restricted"){
    *P = 0.0;
/*

    //---- Handle subsets ------------
    // Alpha
    handle_n_subsets(Fao, Sao, P, prms.fragments.size(), FMO_Fao, FMO_Sao, FMO_P, FMO_C, FMO_E,
                     FMO_bands, FMO_occ, FMO_subsets,FMO_norb, FMO_nocc, 2, eigen_method, prms.pop_opt);
*/

  }// restricted



  if(0){
    cout<<"Total composed density matrix\n";
    cout<<"P_alp=\n"<<*P_alp<<endl;
    cout<<"P_bet=\n"<<*P_bet<<endl;
  }

}//void System::FMO_gen1(Control_Parameters& prms,std::string eigen_method)



