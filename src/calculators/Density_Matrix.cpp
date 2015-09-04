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

MATRIX compute_density_matrix(boost::python::list occ, MATRIX C){

  int Norb = C.num_of_cols;
  vector< pair<int,double> > int_occ;
  MATRIX P(Norb,Norb);

  convert_1(occ,int_occ);

  compute_density_matrix(int_occ, &C, &P);

  return P;
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
  if(eigen_method=="generalized"){    solve_eigen(Norb, Fao, Sao, E,C);    }// generalized
  else if(eigen_method=="standard"){  
    MATRIX* I; I = new MATRIX(Norb,Norb); *I = 0.0; for(int i=0;i<Norb;i++){ I->M[i*Norb+i] = 1.0; }
    solve_eigen(Norb, Fao, I, E,C);       // generalized, but with unit overlap
    delete I;
  }// standard
  if(BM){ bench_t[0].stop(); }

  // Generate and order bands in compressed form from the matrices
  if(BM){ bench_t[1].start(); }
  order_bands(E, bands);
  if(BM){ bench_t[1].stop(); }

  // Populate bands
  if(BM){ bench_t[2].start(); }
  double kT = 0.025; // 300 K
  double etol = 0.0001; // how accurately determine E_f
  populate_bands(Nel, degen, kT, etol, pop_opt, bands, occ);
  if(BM){ bench_t[2].stop(); }

  // Update density matrix
  if(BM){ bench_t[3].start(); }
  compute_density_matrix(occ, C, P);
  if(BM){ bench_t[3].stop(); }

}//void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt, ....
 


void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                              /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,           /*Outputs*/
               int BM, vector<Timer>& bench_t){                                              /*Benchmarking data*/
// Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix:
// 1) solve  Fao * C  = Sao * C * E
// 2) order bands
// 3) compute P as  P = C * N * C.T(), where N = occ

  int Norb = Fao->num_of_cols;
    
  // Get electronic structure (wfc and energies) from given Fock matrix
  if(BM){ bench_t[0].start(); }
  solve_eigen(Norb, Fao, Sao, E,C); 
  if(BM){ bench_t[0].stop(); }

  // Generate and order bands in compressed form from the matrices
  if(BM){ bench_t[1].start(); }
  order_bands(E, bands);
  if(BM){ bench_t[1].stop(); }

  // Populate bands
  if(BM){ bench_t[2].start(); }
  populate_bands(Nel, degen, kT, etol, pop_opt, bands, occ);
  if(BM){ bench_t[2].stop(); }

  // Update density matrix
  if(BM){ bench_t[3].start(); }
  compute_density_matrix(occ, C, P);
  if(BM){ bench_t[3].stop(); }

}//void Fock_to_P(...)

void Fock_to_P(MATRIX* Fao, MATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               MATRIX* E, MATRIX* C, MATRIX* P,                                                         /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ                       /*Outputs*/
              ){       

  int BM = 0; 
  vector<Timer> bench_t;

  Fock_to_P(Fao, Sao, Nel, degen, kT, etol, pop_opt, E, C, P, bands, occ, BM, bench_t);

}





boost::python::list Fock_to_P(MATRIX Fao, MATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt){ 

  int Norb = Fao.num_of_cols;
  MATRIX E(Norb,Norb);
  MATRIX C(Norb,Norb);
  MATRIX P(Norb,Norb);
  vector< pair<int,double> > bands;
  vector< pair<int,double> > occ;


  Fock_to_P(&Fao, &Sao, Nel, degen, kT, etol, pop_opt, &E, &C, &P, bands, occ);

  boost::python::list res;

  res.append(E);
  res.append(C);
  res.append(P);
  res.append(convert_2(bands));
  res.append(convert_2(occ));

  return res;

}


 

}//namespace libcalculators
