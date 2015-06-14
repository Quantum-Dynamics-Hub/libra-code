#include "Density_Matrix.h"
#include "Fermi.h"
#include "Bands.h"

#include "../mmath/libmmath.h"
using namespace libmmath;
using namespace libmmath::libmeigen;


namespace libcalculators{


void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P){
// Scales as O(Norb^3)
// P = C * N * C.T(), where N - is diagonal - populations in MO basis
// occ - represents N

  int a,b,jj,j;

  int Norb = occ.size();
  *P = 0.0;

  int atab = 0;
  for(a=0;a<Norb;a++){

    int btab = 0;
    for(b=0;b<Norb;b++){
      for(jj=0;jj<Norb;jj++){ 

        j = occ[jj].first;       
//      P->M[a*Norb+b] += occ[jj].second*C->M[a*Norb+j]*C->M[b*Norb+j]; // this is what we do below
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

}// void compute_density_matrix(...)



void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt, double kT, double etol,
               MATRIX* Fao, MATRIX* Sao, MATRIX* C, MATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               MATRIX* P, int BM, vector<Timer>& bench_t){
// Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix:
// 1) solve  Fao * C  = Sao * C * E
// 2) order bands
// 3) compute P as  P = C * N * C.T(), where N = occ

    
  // Get electronic structure (wfc and energies) from given Fock matrix
  if(BM){ bench_t[0].start(); }
  if(eigen_method=="generalized"){   solve_eigen(Norb, Fao, Sao, E,C);    }// generalized
  else if(eigen_method=="standard"){  
    MATRIX* I; I = new MATRIX(Norb,Norb); *I = 0.0; for(int i=0;i<Norb;i++){ I->M[i*Norb+i] = 1.0; }
    solve_eigen(Norb, Fao, I, E,C);       // generalized, but with unit overlap
    delete I;
  }// standard
  if(BM){ bench_t[0].stop(); }

  // Generate and order bands in compressed form from the matrices
  if(BM){ bench_t[1].start(); }
  order_bands(Norb, E, bands);
  if(BM){ bench_t[1].stop(); }

  // Populate bands
  if(BM){ bench_t[2].start(); }
  populate_bands(Nocc, Norb, degen, Nel, pop_opt, kT, etol, bands, occ);
  if(BM){ bench_t[2].stop(); }

  // Update density matrix
  if(BM){ bench_t[3].start(); }
  compute_density_matrix(occ, C, P);
  if(BM){ bench_t[3].stop(); }

}//void Fock_to_P(...)
 

}//namespace libcalculators
