#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

using namespace Eigen;


#include "Engine.h"
#include "units.h"
#include "classes.h"


using namespace std;



/****************************************************************************
  This file contains following functions:

  void solve_eigen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C)
  void solve_eigen(int Norb, MATRIX* H, MATRIX* E, MATRIX* C)

  int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out)
  void order_bands(int Norb, MATRIX* E, vector< pair<int,double> >& bands)

  double fermi_integral(vector< pair<int,double> >& bnds, double ef, double degen)
  double fermi_energy(vector< pair<int,double> >& bnds,double Nel,double degen)
  double population(double e,double ef,double degen)
  void populate(int Nocc, int Norb, int degen, double Nel, int pop_opt,vector< pair<int,double> >& bands,vector< pair<int,double> >& occ)

  void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P)

  void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb, MATRIX* Ra, MATRIX* Rb)
  void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb)


  void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt,
               MATRIX* Fao, MATRIX* Sao, MATRIX* C, MATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               MATRIX* P, vector<Timer>&)

  void update_Mull_orb_pop(MATRIX* P, MATRIX* S, vector<double>& Mull_orb_pop_gross,vector<double>& Mull_orb_pop_net)

  void update_Mull_charges(vector<int>& fragment, vector<int>& basis_fo, vector<vector<int> >& at_orbitals,vector<double>& Zeff,
                           vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net,
                           vector<double>& Mull_charges_gross, vector<double>& Mull_charges_net)


  double energy_elec(int Norb,MATRIX* Pao,MATRIX* Hao,MATRIX* Fao)

  void show_bands(int Norb, int Nocc, vector< pair<int,double> >& bands,vector< pair<int,double> >& occ)

****************************************************************************/

void solve_eigen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C){
// Solve H * C = S * C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXd A(Norb,Norb), B(Norb,Norb);
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = 0.5*(H->M[i*Norb+j] + H->M[j*Norb+i]);     
      B(i,j) = 0.5*(S->M[i*Norb+j] + S->M[j*Norb+i]);     
    }// for j
  }// for i


  // Solve eigenvalue problem
  GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B);
  if(solution.info()!=Success){ cout<<"Eigen fails\n";   exit(0); }


  // Copy results into output matrices
  for(i=0;i<Norb;i++){

     E->M[i*Norb+i] = solution.eigenvalues()[i]; 

    for(j=0;j<Norb;j++){  C->M[j*Norb+i] = solution.eigenvectors().col(i)[j];    }// for j

  }// for i

  if(0){
    // Checking solution - for debug purposes, inactive by default
    cout<<"in solve_eigen (with MATRIX):\n";
    cout<<"A = "<<A<<endl;
    cout<<"B = "<<B<<endl;
    cout<<"E = "<<(*E)<<endl;
    cout<<"C = "<<(*C)<<endl;
    cout<<"S = "<<(*S)<<endl;
    cout<<"H = "<<(*H)<<endl;
    cout<<"A*C = "<<(*H) * (*C)<<endl;
    cout<<"S*C*E = "<<(*S) * (*C) * (*E)<<endl;
    cout<<"A*C - S*C*E = "<<*H * (*C) - (*S) * (*C) * (*E)<<endl;
    cout<<"C.T()*S*C= "<<((*C).T()) * (*S) * (*C)<<endl; // This gives unity, indeed.
  }


}//void solve_eigen(int Norb, MATRIX* H, MATRIX* S, MATRIX* E, MATRIX* C)


void solve_eigen(int Norb, MATRIX* H, MATRIX* E, MATRIX* C){
// Solve H * C = C * E
// i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
// C[i][j] - is the weight of j-th AO in i-th MO

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// ACHTUNG:    EigenSolver<MatrixXd> solution(A); - messes up the ordering of eigenvalues
// GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B); does not!, so 
// we eventually would need to abandon the former (special case) solver

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  *E = 0.0;
  *C = 0.0;

  int i,j;

  // Wrapper matrices for Eigen3
  MatrixXd A(Norb,Norb);
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      // Symmetrize, to reduce numerical errors:
      A(i,j) = 0.5*(H->M[i*Norb+j] + H->M[j*Norb+i]);     
    }// for j
  }// for i


  // Solve eigenvalue problem
  EigenSolver<MatrixXd> solution(A);


  // Copy results into output matrices
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

      if(i==j){   E->M[i*Norb+i] = solution.eigenvalues()[i].real();  }// i == j    

      C->M[j*Norb+i] = solution.eigenvectors().col(i)[j].real(); 

    }// for j
  }// for i

  if(0){
    // Checking solution - for debug purposes, inactive by default
    cout<<"in solve_eigen (with MATRIX):\n";
    cout<<"A*C = "<<*H * (*C)<<endl;
    cout<<"C*E = "<<(*C) * (*E)<<endl;
    cout<<"C.T()*C= "<<((*C).T()) * (*C)<<endl; // This gives unity, indeed.
  }

}// void solve_eigen(int Norb, MATRIX* H, MATRIX* E, MATRIX* C)


int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out){

  if(out.size()>0){ out.clear(); }
  int sz = in.size();
  if(sz==0){ }
  else if(sz==1){ out = in; }
  else{
    // Divide in into 2 blocks of approximately same size each
    int half = sz/2;
    vector< pair<int,double> > in1,in2,out1,out2;
    for(int i=0;i<half;i++){ in1.push_back(in[i]); }
    merge_sort(in1,out1);
       
    for(int i=half;i<sz;i++){ in2.push_back(in[i]); }
    merge_sort(in2,out2); 
      
    // Now merge two parts
    int cl,cr; cl = 0; cr = half;
    while((cl<half) && (cr<sz)){
 
    /// EXTREMELY IMPORTANT !!!  This simple, slight difference - the use of < or <= makes HUGE differnece
    /// The "good" version maximally preserves the ordering of orbitals, so one does not run into trouble of alternating
    /// charges - this also leads to symmetric charge distribution in unrestricted formulations even without population smearing
    /// I think it is even more than that - this can lead to convergence (or faster convergence), while the wrong
    /// method may lead to either non-convergent scheme or to sifnificantly slower convergenc.

    if(out1[cl].second<=out2[cr-half].second){ out.push_back(out1[cl]); cl++; }  ///< <-- This is good
//    if(out1[cl].second<out2[cr-half].second){ out.push_back(out1[cl]); cl++; }  ///< <-- Try old

    /// The "bad" version will alternate order of nearby orbitals
    /// It is here only for the purpose of "demonstration of pathological implementation"
//    if(out1[cl].second<out2[cr-half].second){ out.push_back(out1[cl]); cl++; } <-- This is BAD
      else{ out.push_back(out2[cr-half]); cr++; }
    } 
    while(cl<half){ out.push_back(out1[cl]); cl++;}
    while(cr<sz)  { out.push_back(out2[cr-half]); cr++;}
  }
  return 0;

}// int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out)

void order_bands(int Norb, MATRIX* E, vector< pair<int,double> >& bands){

  vector< pair<int,double> > in;
  pair<int,double> x;
  int i;

  for(i=0;i<Norb;i++){
    x = make_pair(i,E->M[i*Norb+i]);
    in.push_back(x);
  }

  if(0){
    cout<<"Unordered bands:\n";
    for(i=0;i<Norb;i++){  cout<<"i= "<<i<<" orb_indx = "<<in[i].first<<" E[i]= "<<in[i].second<<endl; }
  }

  merge_sort(in,bands);  

  if(0){
    cout<<"Ordered bands:\n";
    for(i=0;i<Norb;i++){  cout<<"i= "<<i<<" orb_indx = "<<in[i].first<<" E[i]= "<<in[i].second<<endl; }
  }


  
}// void order_bands(int Norb, MATRIX* E, vector< pair<int,double> >& bands)



double fermi_integral(vector< pair<int,double> >& bnds, double ef, double degen){
// Compute integral(sum):
// sum [degen / (1 + exp(e-ef))]  = N
//  i
// where N - is a number of electrons (valence), 2 - is because we use closed shell formulation - each orbital represents 2 spin-orbitals
// so 2 - is a degeneracy of the energy level

  double kT = 0.00095; // kT in a.u. for T = 300 K
  int Norb = bnds.size();

  double sum = 0.0;
  for(int i=0;i<Norb;i++){

    double argg = ((bnds[i].second-ef)/(1.0*kT));
    double pop = degen;

    if(argg>50.0){   }
    else if(argg<-50.0){   sum += degen;   }
    else{
      pop = degen/(1.0 + exp(argg));
      sum += pop;
    }  


  }// for i

  return sum;

}// double fermi_integral(vector< pair<int,double> >& bnds, double ef, double degen)

double fermi_energy(vector< pair<int,double> >& bnds,double Nel,double degen){
// Computes Fermi energy
  double ef_l,ef_m,ef_r,i_l,i_m,i_r;
  double tol = 1e-5;
  double err = 2.0*tol;
  int Norb = bnds.size();

  ef_l = bnds[0].second - 10.0;
  ef_r = bnds[Norb-1].second + 10.0;

  i_l = fermi_integral(bnds,ef_l,degen) - Nel;
  i_r = fermi_integral(bnds,ef_r,degen) - Nel;

  do{
    ef_m = 0.5*(ef_l + ef_r);
    i_m = fermi_integral(bnds,ef_m,degen) - Nel;

    if(0){
      cout<<"ef_l= "<<ef_l<<" ef_m= "<<ef_m<<" ef_r= "<<ef_r<<endl;
      cout<<"i_l= "<<i_l<<" i_m= "<<i_m<<" i_r= "<<i_r<<endl;
    }

    int var;
    if(i_m*i_r<=0.0 && i_m*i_l>=0.0){ var = 1; }
    else if(i_m*i_r>=0.0 && i_m*i_l<=0.0){ var = 2; }
    else{ cout<<"Error in fermi_energy\n"; exit(0); }

    switch(var){
      case 1: {i_l = i_m; ef_l = ef_m; break; }
      case 2: {i_r = i_m; ef_r = ef_m; break; }
      default: break;
    }

//    err = 0.5*(ef_r-ef_l);
    err = 0.5*(i_r - i_l);
    
  }while(err>tol);
 
  return 0.5*(ef_l+ef_r);

}// double fermi_energy(vector< pair<int,double> >& bnds,double Nel,double degen)


double population(double e,double ef,double degen){
// Compute Fermi-based occupation of doubly-degenerate energy level i

  double kT = 0.00095; // kT in a.u. for T = 300 K

  double argg = ((e-ef)/(1.0*kT));
  double pop = 2.0;

  if(argg>50.0){ pop = 0.0;  }
  else if(argg<-50.0){   pop = degen;   }
  else{ pop = degen/(1.0 + exp(argg));  }  


  return pop;

}// double population(double e,double ef,double degen)


void populate(int Nocc, int Norb, int degen, double Nel, int pop_opt, vector< pair<int,double> >& bands,vector< pair<int,double> >& occ){

  int i;
  occ = bands;  // maybe not the best way, but it allocates memory

  if(pop_opt==0){  // integer populations

    for(i=0;i<Norb;i++){
      occ[i].first = bands[i].first;

      if(i<Nocc){ occ[i].second = (double)degen; }
      else{ occ[i].second = 0.0; }

    }// for i
  }// pop_opt = 0

  else if(pop_opt==1){  // Fermi distribution (fractional populations possible)
 
    double E_f = fermi_energy(bands,degen*Nocc,degen);

    for(i=0;i<Norb;i++){
      occ[i].first = bands[i].first;
      occ[i].second = population(bands[i].second,E_f,degen);

    }
  }// pop_opt = 1

  else if(pop_opt==2){  // Fermi distribution (fractional populations possible) 
  // Fractional total number of electrons is possible
 
    double E_f = fermi_energy(bands,Nel,degen);  // otherwise Nel = degen*Nocc

    for(i=0;i<Norb;i++){
      occ[i].first = bands[i].first;
      occ[i].second = population(bands[i].second,E_f,degen);

    }
  }// pop_opt = 2



  if(0){  // for debug, now inactive
    cout<<"In populate:\n";
    cout<<"Occupation numbers:\n";
    for(i=0;i<Norb;i++){
      cout<<"i= "<<i<<" orbital index= "<<occ[i].first<<" occupation_number= "<<occ[i].second<<endl;
    }
  }// if 


}// void populate(int Nocc, int Norb, int degen, double Nel, int pop_opt,vector< pair<int,double> >& bands,vector< pair<int,double> >& occ)


void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P){
// Scales as O(Norb^3)
// P = C * N * C.T(), where N - is diagonal - populations in MO basis

  int a,b,jj,j;

  int Norb = occ.size();
  *P = 0.0;

  int atab = 0;
  for(a=0;a<Norb;a++){

    int btab = 0;
    for(b=0;b<Norb;b++){
      for(jj=0;jj<Norb;jj++){ 

        j = occ[jj].first;       
//          P->M[a*Norb+b] += occ[jj].second*C->M[a*Norb+j]*C->M[b*Norb+j]; // assume coefficients are real  
        P->M[atab+b] += occ[jj].second*C->M[atab+j]*C->M[btab+j]; // assume coefficients are real  


      }// for jj
      btab += Norb;
      
    }// for b

    atab += Norb;
  }// for a

  // For debug, currently inactive
  if(0){
    cout<<"Density matrix:\n";
    cout<<*P<<endl;
    cout<<"tr(density_matrix) = "<<P->tr()<<endl;
  }// restricted

}// void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P)


void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb, MATRIX* Ra, MATRIX* Rb){
// Spin annihilation
// according to Eqs. 2.19 - 2.21 from Pople, Beveridge, Dobosh JCP 47, 2026 (1967)

  int DF = 0;

  if(DF){   cout<<"in annihilate...\n"; }

  int Norb = Pa->num_of_cols;
  MATRIX* Pab;   Pab   = new MATRIX(Norb,Norb);
  MATRIX* Pba;   Pba   = new MATRIX(Norb,Norb);
  MATRIX* Pabab; Pabab = new MATRIX(Norb,Norb);
  MATRIX* Pbaba; Pbaba = new MATRIX(Norb,Norb);

  *Pab = *Pa * *Pb; 
  *Pba = *Pb * *Pa; 
  *Pabab = *Pab * *Pab; 
  *Pbaba = *Pba * *Pba; 

  double Tr_ab = Pab->tr();
  double Tr_abab = Pabab->tr();

  if(DF){
    cout<<"Tr_ab = "<<Tr_ab<<endl;
    cout<<"Tr_abab = "<<Tr_abab<<endl;
  }

  double p = Na; // number of alpha electrons
  double q = Nb; // number of beta electrons
  double s = 0.5*(p-q); 
  double A = q - 2.0*(s+1);                                                                                              
  double M2 = A*A - 2.0*A*Tr_ab + p*q - (p+q)*Tr_ab + 2.0*Tr_ab + 2.0*Tr_ab*Tr_ab - 2.0*Tr_abab;

  if(DF){
    cout<<"p = "<<p<<endl;
    cout<<"q = "<<q<<endl;
    cout<<"s = "<<s<<endl;
    cout<<"A = "<<A<<endl;
    cout<<"M2 = "<<M2<<endl;
  }
  

  // Alpha   
  *Ra =  ( M2 + Tr_ab - q ) * (*Pa)
       + ( p - Tr_ab ) * (*Pb)
       + ( *Pb * *Pab ) 
       + ((p+q) - 4.0*Tr_ab - 3.0 + 2.0*A) * (*Pa * *Pba)
       + ( 2.0*Tr_ab - p + 1.0 - A)*(*Pab + *Pba)
       - 2.0 * (*Pabab + *Pbaba)
       + 4.0 * (*Pabab * *Pa);

  *Ra = (1.0/M2) * (*Ra);

  if(DF){
    cout<<"Original alpha density matrix =\n"<<*Pa<<endl;
    cout<<"Annihilated alpha density matrix =\n"<<*Ra<<endl;
  }


  // Beta
  *Rb =  ( M2 + Tr_ab - p ) * (*Pb)
       + ( q - Tr_ab ) * (*Pa)
       + ( *Pa * *Pba ) 
       + ((p+q) - 4.0*Tr_ab - 3.0 + 2.0*A) * (*Pb * *Pab)
       + ( 2.0*Tr_ab - q + 1.0 - A)*(*Pba + *Pab)
       - 2.0 * (*Pbaba + *Pabab)
       + 4.0 * (*Pbaba * *Pb);

  *Rb = (1.0/M2) * (*Rb);

  if(DF){
    cout<<"Original beta density matrix =\n"<<*Pb<<endl;
    cout<<"Annihilated beta density matrix =\n"<<*Rb<<endl;
  }


  // Clean temporary objects
  delete Pab;
  delete Pba;
  delete Pabab;
  delete Pbaba;

}

void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb){

  int Norb = Pa->num_of_cols;
  MATRIX* Ra;   Ra   = new MATRIX(Norb,Norb);
  MATRIX* Rb;   Rb   = new MATRIX(Norb,Norb);

  annihilate(Na, Nb, Pa, Pb, Ra, Rb);


  *Pa = *Ra;
  *Pb = *Rb;

  
  delete Ra;
  delete Rb;

}


void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt,
               MATRIX* Fao, MATRIX* Sao, MATRIX* C, MATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               MATRIX* P, vector<Timer>& bench_t){
// Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
// In these steps there is no coupling of spin-up and spin-down channels, so they can
// be solved one by one, independently
//exit(0);
  int BM = 1;
    
  // Get electronic structure (wfc and energies) from given Fock matrix
  if(BM){ bench_t[0].start(); }
  if(eigen_method=="generalized"){    ::solve_eigen(Norb, Fao, Sao, E,C);    }// generalized
  else if(eigen_method=="standard"){  
    MATRIX* I; I = new MATRIX(Norb,Norb); *I = 0.0; for(int i=0;i<Norb;i++){ I->M[i*Norb+i] = 1.0; }
    ::solve_eigen(Norb, Fao, I, E,C);       // generalized, but with unit overlap
    delete I;
  }// standard
  if(BM){ bench_t[0].stop(); }

  // Generate and order bands in compressed form from the matrices
  if(BM){ bench_t[1].start(); }
  ::order_bands(Norb, E, bands);
  if(BM){ bench_t[1].stop(); }

  // Populate bands
  if(BM){ bench_t[2].start(); }
  populate(Nocc, Norb, degen, Nel, pop_opt, bands, occ);
  if(BM){ bench_t[2].stop(); }

  // Update density matrix
  if(BM){ bench_t[3].start(); }
  compute_density_matrix(occ, C, P);
  if(BM){ bench_t[3].stop(); }

}//void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt, ....
 


void update_Mull_orb_pop(MATRIX* P, MATRIX* S, vector<double>& Mull_orb_pop_gross,vector<double>& Mull_orb_pop_net){
// Compute Mulliken orbital-resolved populations from the density matrix
// memory for Mull_orb_pop_*  is assumed already allocated
// gross and net populations are computed

  int a,b;
  double tmp_a, tmp_ab;

  int Norb = Mull_orb_pop_gross.size(); // is assumed be the same as for the other array

  for(a=0;a<Norb;a++){
    Mull_orb_pop_gross[a] = 0.0;
    Mull_orb_pop_net[a] = 0.0;
  }

  for(a=0;a<Norb;a++){

    tmp_a = 0.0;
    for(b=0;b<Norb;b++){

      tmp_ab = P->M[a*Norb+b] * S->M[a*Norb+b];
      tmp_a += tmp_ab;

      if(b==a){ Mull_orb_pop_net[a] = tmp_ab; }

    }// for b

    Mull_orb_pop_gross[a] = tmp_a;  

  }// a


/*
  // Alternative definition - test only !!!
  for(a=0;a<Norb;a++){
    for(b=0;b<Norb;b++){

      tmp_ab = P->M[a*Norb+b] * S->M[a*Norb+b];
      Mull_orb_pop_gross[a] += 0.5*tmp_ab;
      Mull_orb_pop_gross[b] += 0.5*tmp_ab;

      if(b==a){ Mull_orb_pop_net[a] = tmp_ab; }

    }// for b
  }// a
*/


}// void update_Mull_orb_pop(MATRIX* P, MATRIX* S, vector<double>& Mull_orb_pop_gross,vector<double>& Mull_orb_pop_net)


void update_Mull_charges(vector<int>& fragment, vector<int>& basis_fo, vector<vector<int> >& at_orbitals,vector<double>& Zeff,
                         vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net,
                         vector<double>& Mull_charges_gross, vector<double>& Mull_charges_net){
// fragment - contains indexes of all nuclei, for which the Mull charges should be updated
// at_orbitals - contains indexes of AO for each of all nuclei in the system  - size is Nnucl
// len(Mull_orb_pop_*) = size of fragment basis = len(basis_fo)
// len(Mull_charges_*) = size of global nuclear array = len(Zeff)


  int Nat_frag = fragment.size(); 

  // Init arrays - only nuclear charges
  for(int a=0;a<Nat_frag;a++){  // O(Nfrag) 
    int n = fragment[a]; // index of a-th nucleus in global array
    Mull_charges_gross[n] = Zeff[n];
    Mull_charges_net[n] = Zeff[n];

    for(int i=0;i<at_orbitals[n].size();i++){  // O(Nnucl)
      int j = at_orbitals[n][i];

      for(int k=0;k<basis_fo.size();k++){      // O(Nbas)
        if(basis_fo[k]==j){  

          Mull_charges_gross[n] -= Mull_orb_pop_gross[k];
          Mull_charges_net[n] -= Mull_orb_pop_net[k];
        }// if
      }// for k
    }// for i
  }// for n

}// void update_Mull_charges(vector<double>& Mull_orb_pop_gross, vector<double>& Mull_orb_pop_net, ...


 
//Deprecate this version!!!
double energy_elec(int Norb,MATRIX* Pao,MATRIX* Hao,MATRIX* Fao){

  cout<<"This version of energy_elec() is going to be deprecated...\n";
  cout<<"If you see this message, then some part of the code you are using should be updated to include newer version of energy_elec() function...\n"; 
  exit(0);

// Compute electronic energy (true for HF-derived methods: HF, CNDO, CNDO/2, INDO)
// this general formula is also true for EHT (F = Hcore, so the energy is simply a weighted sum of the eigenvalues)
  
  int i,ii,j;
  double Eelec = 0.0;

  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

      Eelec += Pao->M[i*Norb+j]*(Hao->M[i*Norb+j] + Fao->M[i*Norb+j]);

    }// for i
  }// for j

  Eelec *= 0.5;

  return Eelec;

}// double energy_elec(int Norb,MATRIX* Pao,MATRIX* Hao,MATRIX* Fao)



double energy_elec(int Norb, 
                   MATRIX* P_alp, MATRIX* P_bet, 
                   MATRIX* Hao_alp, MATRIX* Hao_bet,
                   MATRIX* Fao_alp, MATRIX* Fao_bet,
                   MATRIX* dFao_alp_dP_alp, MATRIX* dFao_alp_dP_bet,
                   MATRIX* dFao_bet_dP_alp, MATRIX* dFao_bet_dP_bet,
                   MATRIX* temp
                  ){
// Compute electronic energy (true for HF-derived methods: HF, CNDO, CNDO/2, INDO)
// this general formula is also true for EHT (F = Hcore, so the energy is simply a weighted sum of the eigenvalues)
// temp - is just a temporary array - preallocate it before calling this function
//        this will give some acceleration if the energy function is called very often

// Scaling - O(N^3)
  
  int i,ii,j;
  double Eelec_alp = 0.0;
  double Eelec_bet = 0.0;
  double Eelec_tot = 0.0;


  *temp = 0.0; 
  if(dFao_alp_dP_alp->max_elt()>1e-10){   *temp += *P_alp * *dFao_alp_dP_alp;    }
  if(dFao_alp_dP_bet->max_elt()>1e-10){   *temp += *P_bet * *dFao_alp_dP_bet;    }
   
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      Eelec_alp += P_alp->M[i*Norb+j]*(Hao_alp->M[i*Norb+j] +( Fao_alp->M[i*Norb+j] + temp->M[i*Norb+j]));
    }// for i
  }// for j
  Eelec_alp *= 0.5;


  *temp = 0.0; 
  if(dFao_bet_dP_alp->max_elt()>1e-10){   *temp += *P_alp * *dFao_bet_dP_alp;    }
  if(dFao_bet_dP_bet->max_elt()>1e-10){   *temp += *P_bet * *dFao_bet_dP_bet;    }
   
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      Eelec_bet += P_bet->M[i*Norb+j]*(Hao_bet->M[i*Norb+j] +( Fao_bet->M[i*Norb+j] + temp->M[i*Norb+j]));
    }// for i
  }// for j
  Eelec_bet *= 0.5;


  Eelec_tot = Eelec_alp + Eelec_bet;


  return Eelec_tot;

}// double energy_elec(...)


void show_bands(int Norb, int Nocc, vector< pair<int,double> >& bands,vector< pair<int,double> >& occ){

  int sz = Norb;//bands.size();
  cout<<"# of bands = "<<sz<<"\n";
  cout<<"# of occupied = "<<Nocc<<endl;
//  cout<<"Fermi energy (eV) = "<<E_f/eV<<endl;

  for(int i=0;i<sz;i++){
    cout<<" band# = "<<i<<" mo indx ="<<bands[i].first
        <<" energy(eV) = "<<bands[i].second/eV
        <<" energy(a.u.) = "<<bands[i].second
        <<" occupation = "<<occ[i].second<<endl;
  }

}// void show_bands(int Norb, int Nocc, vector< pair<int,double> >& bands,vector< pair<int,double> >& occ)
