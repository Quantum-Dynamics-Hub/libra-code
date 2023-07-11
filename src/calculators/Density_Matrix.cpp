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
/**
  \file Density_Matrix.cpp
  \brief The file implements functions for density matrix and Fock-to-density calculations
    
*/

#include "Density_Matrix.h"
#include "Fermi.h"
#include "Bands.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libmeigen;


/// libcalculators namespace
namespace libcalculators{

void compute_density_matrix(vector< pair<int,double> >& occ, MATRIX* C, MATRIX* P){
/**
  \brief Straightforward density matrix computation

  Scales as O(Norb^3)
  P = C * N * C.T(), where N - is diagonal - populations in MO basis

  \param[in] occ Occupations of the molecular orbitals (list of MO-index/MO-occupation pairs), represents N
  \param[in] C Pointer the matrix of the MO-LCAO coefficients
  \param[out] P Pointer to the density matrix (in AO basis)

*/


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


void compute_density_matrix(vector< pair<int,double> >& occ, CMATRIX* C, CMATRIX* P){
/**
  \brief Straightforward density matrix computation - complex-valued version

  Scales as O(Norb^3)
  P = C * N * C.H(), where N - is diagonal - populations in MO basis

  \param[in] occ Occupations of the molecular orbitals (list of MO-index/MO-occupation pairs), represents N
  \param[in] C Pointer to the matrix of the MO-LCAO coefficients
  \param[out] P Pointer to the density matrix (in AO basis)

*/


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
        P->M[atab+b] += C->M[atab+j] * occ[jj].second * std::conj(C->M[btab+j]); // assume coefficients are real

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
/**
  \brief Straightforward density matrix computation (Python-friendly)

  Scales as O(Norb^3)
  P = C * N * C.T(), where N - is diagonal - populations in MO basis

  \param[in] occ Occupations of the molecular orbitals (Python list of MO-index/MO-occupation pairs), represents N 
  \param[in] C The matrix of the MO-LCAO coefficients
  Density matrix (in AO basis) is returned as a MATRIX object

*/

  int Norb = C.n_cols;
  vector< pair<int,double> > int_occ;
  MATRIX P(Norb,Norb);

  convert_1(occ,int_occ);

  compute_density_matrix(int_occ, &C, &P);

  return P;
}

CMATRIX compute_density_matrix(boost::python::list occ, CMATRIX C){
/**
  \brief Straightforward density matrix computation (Python-friendly, complex-valued)

  Scales as O(Norb^3)
  P = C * N * C.H(), where N - is diagonal - populations in MO basis

  \param[in] occ Occupations of the molecular orbitals (Python list of MO-index/MO-occupation pairs), represents N 
  \param[in] C The matrix of the MO-LCAO coefficients
  Density matrix (in AO basis) is returned as a MATRIX object

*/

  int Norb = C.n_cols;
  vector< pair<int,double> > int_occ;
  CMATRIX P(Norb,Norb);

  convert_1(occ,int_occ);

  compute_density_matrix(int_occ, &C, &P);

  return P;
}



void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt,
               MATRIX* Fao, MATRIX* Sao, MATRIX* C, MATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               MATRIX* P, vector<Timer>& bench_t){
/**
  \brief Set of instructions to compute density matrix from the Fock matrix

  This is a somewhat older version

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.T(), where N = occ


  \param[in] Norb The number of orbitals = the dimensionality of electronic problem
  \param[in] Nocc The number of occupied (integer occupation) orbitals
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] Nel The number of electrons
  \param[in] eigen_method The flag controlling the assumer orthogonality of basis orbitals:
             eigen_method = "generalized" - non-orthogonal orbitals are assumed, solve a generalized eigenvalue problem
             eigen_method = "standard" - orthogonal orbitals are assumed, solve standard eigenvalue problem
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[in] Fao The pointer to the Fock matrix
  \param[in] Sao The pointer to the AO overlap matrix
  \param[in,out] C The pointer to MO-LCAO matrix
  \param[out] E The pointer to the eigenvalues (of the Fock operator) matrix
  \param[in,out] bands The orbital energies in the vector of pairs format
  \param[in,out] occ The orbital occupancies (MO basis populations) in the vector of pairs format
  \param[in,out] P The pointer to the density matrix 
  \param[in,out] bench_t The benchmarking information: bench_t[0] - eigenvalue solvers, bench_t[1] - ordering bands,
                 bench_t[2] - populate bands, bench_t[3] - density matrix computations
*/

  int BM = 1;
    
  // Get electronic structure (wfc and energies) from given Fock matrix
  if(BM){ bench_t[0].start(); }
  if(eigen_method=="generalized"){   
   solve_eigen(Fao, Sao, E, C, 0);   
  }// generalized
  else if(eigen_method=="standard"){  
    solve_eigen(Fao, E, C, 0);       // generalized, but with unit overlap
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


void Fock_to_P(int Norb,int Nocc, int degen, double Nel, std::string eigen_method, int pop_opt,
               CMATRIX* Fao, CMATRIX* Sao, CMATRIX* C, CMATRIX* E,
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,
               CMATRIX* P, vector<Timer>& bench_t){
/**
  \brief Set of instructions to compute density matrix from the Fock matrix - assuming the complex-valued Fock and MO matrices

  This is a somewhat older version

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.H(), where N = occ


  \param[in] Norb The number of orbitals = the dimensionality of electronic problem
  \param[in] Nocc The number of occupied (integer occupation) orbitals
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] Nel The number of electrons
  \param[in] eigen_method The flag controlling the assumer orthogonality of basis orbitals:
             eigen_method = "generalized" - non-orthogonal orbitals are assumed, solve a generalized eigenvalue problem
             eigen_method = "standard" - orthogonal orbitals are assumed, solve standard eigenvalue problem
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[in] Fao The pointer to the Fock matrix
  \param[in] Sao The pointer to the AO overlap matrix
  \param[in,out] C The pointer to MO-LCAO matrix
  \param[out] E The pointer to the eigenvalues (of the Fock operator) matrix
  \param[in,out] bands The orbital energies in the vector of pairs format
  \param[in,out] occ The orbital occupancies (MO basis populations) in the vector of pairs format
  \param[in,out] P The pointer to the density matrix 
  \param[in,out] bench_t The benchmarking information: bench_t[0] - eigenvalue solvers, bench_t[1] - ordering bands,
                 bench_t[2] - populate bands, bench_t[3] - density matrix computations
*/

  int BM = 1;
    
  // Get electronic structure (wfc and energies) from given Fock matrix
  if(BM){ bench_t[0].start(); }
  if(eigen_method=="generalized"){   
   solve_eigen(Fao, Sao, E, C, 0);   
  }// generalized
  else if(eigen_method=="standard"){  
    solve_eigen(Fao, E, C, 0);       // generalized, but with unit overlap
  }// standard
  if(BM){ bench_t[0].stop(); }

  // Generate and order bands in compressed form from the matrices
  if(BM){ bench_t[1].start(); }
  order_bands(&E->real(), bands);
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
/**
  \brief Set of instructions to compute density matrix from the Fock matrix

  This is a newer version - takes less argiments and makes some inferrences

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.T(), where N = occ

  \param[in] Fao The pointer to the Fock matrix
  \param[in] Sao The pointer to the AO overlap matrix
  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \params[in] kT  Broadening factor for Fermi distribution
  \params[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[out] E The pointer to the eigenvalues (of the Fock operator) matrix
  \param[in,out] C The pointer to MO-LCAO matrix
  \param[in,out] P The pointer to the density matrix 
  \param[in,out] bands The orbital energies in the vector of pairs format
  \param[in,out] occ The orbital occupancies (MO basis populations) in the vector of pairs format
  \param[in] BM Benchmark flag: = 0 - don't do benchmarking, 1 - do it
  \param[in,out] bench_t The benchmarking information: bench_t[0] - eigenvalue solvers, bench_t[1] - ordering bands,
                 bench_t[2] - populate bands, bench_t[3] - density matrix computations
*/
                                                                     

  int Norb = Fao->n_cols;
    
  // Get electronic structure (wfc and energies) from given Fock matrix
  if(BM){ bench_t[0].start(); }
  solve_eigen(Fao, Sao, E,C, 0); 
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


void Fock_to_P(CMATRIX* Fao, CMATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               CMATRIX* E, CMATRIX* C, CMATRIX* P,                                              /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ,           /*Outputs*/
               int BM, vector<Timer>& bench_t){                                              /*Benchmarking data*/
/**
  \brief Set of instructions to compute density matrix from the Fock matrix - complex-valued version

  This is a newer version - takes less argiments and makes some inferrences

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.H(), where N = occ

  \param[in] Fao The pointer to the Fock matrix
  \param[in] Sao The pointer to the AO overlap matrix
  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \params[in] kT  Broadening factor for Fermi distribution
  \params[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[out] E The pointer to the eigenvalues (of the Fock operator) matrix
  \param[in,out] C The pointer to MO-LCAO matrix
  \param[in,out] P The pointer to the density matrix 
  \param[in,out] bands The orbital energies in the vector of pairs format
  \param[in,out] occ The orbital occupancies (MO basis populations) in the vector of pairs format
  \param[in] BM Benchmark flag: = 0 - don't do benchmarking, 1 - do it
  \param[in,out] bench_t The benchmarking information: bench_t[0] - eigenvalue solvers, bench_t[1] - ordering bands,
                 bench_t[2] - populate bands, bench_t[3] - density matrix computations
*/
                                                                     

  int Norb = Fao->n_cols;
    
  // Get electronic structure (wfc and energies) from given Fock matrix
  if(BM){ bench_t[0].start(); }
  solve_eigen(Fao, Sao, E,C, 0); 
  if(BM){ bench_t[0].stop(); }

  // Generate and order bands in compressed form from the matrices
  if(BM){ bench_t[1].start(); }
  order_bands(&E->real(), bands);
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
/**
  \brief Set of instructions to compute density matrix from the Fock matrix

  This is a newer version - even simpler: no benchmarking

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.T(), where N = occ

  \param[in] Fao The pointer to the Fock matrix
  \param[in] Sao The pointer to the AO overlap matrix
  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[out] E The pointer to the eigenvalues (of the Fock operator) matrix
  \param[in,out] C The pointer to MO-LCAO matrix
  \param[in,out] P The pointer to the density matrix 
  \param[in,out] bands The orbital energies in the vector of pairs format
  \param[in,out] occ The orbital occupancies (MO basis populations) in the vector of pairs format
*/


  int BM = 0; 
  vector<Timer> bench_t;

  Fock_to_P(Fao, Sao, Nel, degen, kT, etol, pop_opt, E, C, P, bands, occ, BM, bench_t);

}


void Fock_to_P(CMATRIX* Fao, CMATRIX* Sao, double Nel, double degen, double kT, double etol, int pop_opt, /*Inputs*/
               CMATRIX* E, CMATRIX* C, CMATRIX* P,                                                         /*Outputs*/
               vector< pair<int,double> >& bands, vector< pair<int,double> >& occ                       /*Outputs*/
              ){       
/**
  \brief Set of instructions to compute density matrix from the Fock matrix - complex-valued version

  This is a newer version - even simpler: no benchmarking

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.H(), where N = occ

  \param[in] Fao The pointer to the Fock matrix
  \param[in] Sao The pointer to the AO overlap matrix
  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[out] E The pointer to the eigenvalues (of the Fock operator) matrix
  \param[in,out] C The pointer to MO-LCAO matrix
  \param[in,out] P The pointer to the density matrix 
  \param[in,out] bands The orbital energies in the vector of pairs format
  \param[in,out] occ The orbital occupancies (MO basis populations) in the vector of pairs format
*/


  int BM = 0; 
  vector<Timer> bench_t;

  Fock_to_P(Fao, Sao, Nel, degen, kT, etol, pop_opt, E, C, P, bands, occ, BM, bench_t);

}





boost::python::list Fock_to_P(MATRIX Fao, MATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt){ 
/**
  \brief Set of instructions to compute density matrix from the Fock matrix

  The simplest and Python-friendly version

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.T(), where N = occ

  \param[in] Fao The Fock matrix
  \param[in] Sao The AO overlap matrix
  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  Returns the list of the objects: res[0] = E (eigenvalues matrix), res[1] = C (eigenvectors matrix),
  res[2] = P (density matrix), res[3] = bands (energies, list of lists), res[4] = occ (occupations, list of lists)
*/


  int Norb = Fao.n_cols;
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


boost::python::list Fock_to_P(CMATRIX Fao, CMATRIX Sao, double Nel, double degen, double kT, double etol, int pop_opt){ 
/**
  \brief Set of instructions to compute density matrix from the Fock matrix - complex-valued version

  The simplest and Python-friendly version

  Iterative unit: from a given Hamiltonian (Fock matrix) we obtain density matrix
  In these steps there is no coupling of spin-up and spin-down channels, so they can
  be solved one by one, independently. 
  1) solve  Fao * C  = Sao * C * E
  2) order bands
  3) compute P as  P = C * N * C.H(), where N = occ

  \param[in] Fao The Fock matrix
  \param[in] Sao The AO overlap matrix
  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  Returns the list of the objects: res[0] = E (eigenvalues matrix), res[1] = C (eigenvectors matrix),
  res[2] = P (density matrix), res[3] = bands (energies, list of lists), res[4] = occ (occupations, list of lists)
*/


  int Norb = Fao.n_cols;
  CMATRIX E(Norb,Norb);
  CMATRIX C(Norb,Norb);
  CMATRIX P(Norb,Norb);
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
} /// liblibra

