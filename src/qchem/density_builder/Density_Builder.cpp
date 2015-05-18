#include "Density_Builder.h"
#include "System.h"   // fow a while, untill i move the generic functions like Fock_to_P out of System.h

/****************************************************************************
  This file contains following functions:
                                                          
  void pop_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset)
  void push_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset)

  void handle_subset(MATRIX* F, MATRIX* S, MATRIX* P, 
                     MATRIX* f, MATRIX* s, MATRIX* p, MATRIX* c, MATRIX* e,
                     vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,                  
                     vector<int>& subset,
                     int norb, int nocc, int degen, std::string eigen_method, int pop_opt)

  void handle_n_subsets(MATRIX* F, MATRIX* S, MATRIX* P, 
                        int n_subs,
                        MATRIX** f, MATRIX** s, MATRIX** p, MATRIX** c, MATRIX** e,
                        vector< vector< pair<int,double> > >& bands, vector< vector< pair<int,double> > >& occ, 
                        vector< vector<int> >& subset,
                        vector<int>& norb, vector<int>& nocc, int degen, std::string eigen_method, int pop_opt)

****************************************************************************/



void pop_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset){
// Extract the submatrix x from the matrix X according to indices given in <subset>
// Assume that memory for x is already allocated and its dimensions are consistent with the map dimensions:
// subset_size == x->num_of_cols = x->num_of_rows = subset.size()
// X->num_of_cols = X->num_of_rows >= subset_size

  int N = X->num_of_cols;
  int n = x->num_of_cols;

  int i,j,a,b;

  for(i=0;i<n;i++){
    a = subset[i];
    for(j=0;j<n;j++){      
      b = subset[j];
      x->M[i*n+j] = X->M[a*N+b];
    }// j
  }// i

}// void pop_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset)

void push_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset){
// Pushes the smaller submatrix x back to the bigger matrix X, according to indices given in <subset>
// Assume that memory for x is already allocated and its dimensions are consistent with the map dimensions:
// subset_size == x->num_of_cols = x->num_of_rows = subset.size()
// X->num_of_cols = X->num_of_rows >= subset_size

  int N = X->num_of_cols;
  int n = x->num_of_cols;

  int i,j,a,b;

  for(i=0;i<n;i++){
    a = subset[i];
    for(j=0;j<n;j++){      
      b = subset[j];
      X->M[a*N+b] = x->M[i*n+j];
    }// j
  }// i

}// void push_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset)




void handle_subset(MATRIX* F, MATRIX* S, MATRIX* P, 
                   MATRIX* f, MATRIX* s, MATRIX* p, MATRIX* c, MATRIX* e,
                   vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,                  
                   vector<int>& subset,
                   int norb, int nocc, int degen, std::string eigen_method, int pop_opt){
// Generic procedure for composing densities for given subset and placing them back
// norb - dimension of smaller matrices: f, s, p, c, e
// nocc - number of occupied orbitals out of the norb
// all temporary matrices f,s,p,c,e are assumed to be well-prepared (right dimensions, memory is allocated)

  //------------- Extract block for AOs from the whole system matrix ------------
  pop_submatrix(F, f, subset);
  pop_submatrix(S, s, subset);


  //------------- Compute density matrix --------------- 
  Fock_to_P(norb, nocc, degen, nocc, eigen_method, pop_opt, f, s, c, e, bands, occ, p);

  //------------- Map sub-set densities back to global matrix ----------------
  // This procedure works only if <subset> sets do not overlap (e.g. for SAD), 
  // otherwise one would need some combination procedure

  if(0){
    cout<<"subset = "; for(int i=0;i<subset.size();i++){ cout<<subset[i]<<"  "; } cout<<endl;
    cout<<"subset F= "<<*f<<endl;
    cout<<"subset S= "<<*s<<endl;
    cout<<"resulting P= "<<*p<<endl;
  }


  push_submatrix(P, p, subset);


}// void handle_subset(MATRIX* F, MATRIX* S, MATRIX* P ...



void handle_n_subsets(MATRIX* F, MATRIX* S, MATRIX* P, 
                      int n_subs,
                      MATRIX** f, MATRIX** s, MATRIX** p, MATRIX** c, MATRIX** e,
                      vector< vector< pair<int,double> > >& bands, vector< vector< pair<int,double> > >& occ, 
                      vector< vector<int> >& subset,
                      vector<int>& norb, vector<int>& nocc, int degen, std::string eigen_method, int pop_opt){

  for(int n=0;n<n_subs;n++){

    handle_subset(F, S, P, 
                  f[n], s[n], p[n], c[n], e[n],  bands[n], occ[n], subset[n], norb[n], nocc[n], 
                  degen, eigen_method, pop_opt);

  }// for n

}// void handle_n_subsets(MATRIX* F, MATRIX* S, MATRIX* P...



void sad(){
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



