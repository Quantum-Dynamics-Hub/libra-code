/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Basis_ovlp.cpp
  \brief The file implements function for updating the AO overlap matrix
    
*/

#include "Basis.h"
#include "../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libmeigen;
using namespace libqobjects;


/// libbasis namespace
namespace libbasis{



void update_overlap_matrix(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
                           vector<AO>& basis_ao, MATRIX& Sao){
/**
  \brief Update the oberlap matrix (in AO basis): <AO(i)|AO(j)>
  \param[in] x_period Then number of periodic shells in X direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] y_period Then number of periodic shells in Y direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] z_period Then number of periodic shells in Z direction: 0 - only the central shell, 1 - [-1,0,1], etc.
  \param[in] t1 The periodicity vector along a crystal direction ("X")
  \param[in] t2 The periodicity vector along b crystal direction ("Y")
  \param[in] t3 The periodicity vector along c crystal direction ("Z")
  \param[in] basis_ao The list of all AOs (basis)
  \param[out] Sao The output overlap matrix

  This function can also take periodic images of the system into account
*/


  int i,j,n,I,J;
  VECTOR dIdA,dIdB,TV, Rij,Rij0;
  double dist, dist_min;  
  int opt_x,opt_y,opt_z;

  int Norb = basis_ao.size();

  Sao = 0.0;
 
  for(i=0;i<Norb;i++){
    for(j=i;j<Norb;j++){      

      for(int nx=-x_period;nx<=x_period;nx++){
        for(int ny=-y_period;ny<=y_period;ny++){
          for(int nz=-z_period;nz<=z_period;nz++){

            // This summation corresponds to k = 0 (Gamma-point)    
            TV = nx*t1 + ny*t2 + nz*t3;

            if(i==j){
              AO tmp_ao(basis_ao[i]); 
              tmp_ao.shift_position(TV);
              Sao.M[i*Norb+j] += gaussian_overlap(basis_ao[i],tmp_ao); //,0,dIdA,dIdB,aux,n_aux,TV);
            }
            else{

              basis_ao[j].shift_position(TV);
              Sao.M[i*Norb+j] += gaussian_overlap(basis_ao[i],basis_ao[j]); //,0,dIdA,dIdB,aux,n_aux,TV);
              basis_ao[j].shift_position(-TV);

            }// i != j

          }// for nz
        }// for ny
      }// for nx 

      Sao.M[j*Norb+i] = Sao.M[i*Norb+j];

    }// j
  }// i



/*   !!!!!!!!!!!!!! Don't need the below part untill we get to the gradients !!!!!!!!!!!

  else{

    for(n=0;n<Nnucl;n++){
      *dSao_dx[n] = 0.0;
      *dSao_dy[n] = 0.0;
      *dSao_dz[n] = 0.0;
    }

  
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

        Sao->M[i*Norb+j] = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

        dSao_dx[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.x;
        dSao_dx[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.x;

        dSao_dy[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.y;
        dSao_dy[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.y;

        dSao_dz[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.z;
        dSao_dz[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.z;
*/
/*
      // Tests
      double dx = 0.0001;
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      double s2 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] -= VECTOR(2.0*dx,0.0,0.0);
      double s1 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

      cout<<"Overlap:  der(numer)= "<<0.5*(s2-s1)/dx<<" der(anal)= "<<dIdA.x<<endl;
*/
     

//    }// for j
//  }// for i
//  }// else : method!=eht0

//  for(int n=0;n<Nnucl;n++){
//    cout<<"dSao_dx["<<n<<"]= \n"<<*dSao_dx[n]<<endl;
//  }

//  cout<<"Overlap matrix in AO basis:\n";
//  cout<<*Sao<<endl;
  // So far use my inverter:
//  Sao->Inverse(Sao_inv);

}


void pop_cols(MATRIX& X, MATRIX& x, vector<int>& cols){
// Copies selected columns from X to x 

  int nrows = X.n_rows;
  int sz = cols.size();

  for(int i=0;i<sz;i++){
    for(int n=0; n<nrows; n++){
      x.set(n,i, X.get(n, cols[i]));
    }
  }

}
void pop_cols(CMATRIX& X, CMATRIX& x, vector<int>& cols){
// Copies selected columns from X to x 

  int nrows = X.n_rows;
  int sz = cols.size();

  for(int i=0;i<sz;i++){
    for(int n=0; n<nrows; n++){
      x.set(n,i, X.get(n, cols[i]));
    }
  }

}



void MO_overlap(MATRIX& Smo, vector<AO>& ao_i, vector<AO>& ao_j, MATRIX& Ci, MATRIX& Cj,
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2){
/**
  Finds the keywords and their patterns and extracts the parameters
  \param[in] ao_i, ao_j : atomic orbital basis at different time step.
  \param[in] Ci, Cj : molecular coefficients at different time step.
  \param[in] active_orb_i, active_orb_j The lists of the active orbitals for each MO-AO data set
  This function returns overlap matrix of atomic orbitals with different time step
  like <MO(t)|MO(t+dt)>.

*/

  int Nbas_i = ao_i.size();
  int Nbas_j = ao_j.size();

  // Ci = Nbas_i x Norb_i 
  // Cj = Nbas_j x Norb_j
  int Norb_i = Ci.n_cols;
  int Norb_j = Cj.n_cols;

  int Norb_i_act = active_orb_i.size();  // the number of active orbitals in set i
  int Norb_j_act = active_orb_j.size();  // the number of active orbitals in set j


  // Smo = MATRIX(Norb_i_act, Norb_j_act) <- this is what expected
  if(Smo.n_rows != Norb_i_act){
    cout<<"Error in MO_overlap: The number of rows of the MO overlap matrix"<<Smo.n_rows
        <<" is not consistent with the number of <i| orbitals "<<Norb_i_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }
  if(Smo.n_cols != Norb_j_act){
    cout<<"Error in MO_overlap: The number of cols of the MO overlap matrix"<<Smo.n_cols
        <<" is not consistent with the number of |j> orbitals "<<Norb_j_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }



  // Allocate working memory
  int i,j;
  int is_derivs = 0;
  int is_normalize = 0;
  VECTOR dIdA, dIdB;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

    
  MATRIX* ci; ci = new MATRIX(Nbas_i, Norb_i_act);
  MATRIX* cj; cj = new MATRIX(Nbas_j, Norb_j_act);
  MATRIX* Sao; Sao = new MATRIX(Nbas_i, Nbas_j);

  pop_cols(Ci, *ci, active_orb_i);
  pop_cols(Cj, *cj, active_orb_j);  


  // overlap matrix of S
  for(i=0;i<Nbas_i;i++){  // all orbitals
    for(j=0;j<Nbas_j;j++){

      double dist2 = (ao_i[i].primitives[0].R - ao_j[j].primitives[0].R).length2();       

      if(dist2<max_d2) {  
        double res = gaussian_overlap(ao_i[i], ao_j[j], is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);
        Sao->set(i,j, res); 
      }
      else{  Sao->set(i,j,0.0);  }


    }// for j
  }// for i


  Smo = (*ci).T() * (*Sao) * (*cj);   // ( Norb_i_act x Nbas_i ) x (Nbas_i x Nbas_j) x (Nbas_j x Norb_j_act) = Norb_i_act x Norb_j_act


  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();

  delete Sao;
  delete ci;
  delete cj;


}


void MO_overlap(CMATRIX& Smo, vector<AO>& ao_i, vector<AO>& ao_j, CMATRIX& Ci, CMATRIX& Cj,
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2){
/**
  Finds the keywords and their patterns and extracts the parameters
  \param[in] ao_i, ao_j : atomic orbital basis at different time step.
  \param[in] Ci, Cj : molecular coefficients at different time step.
  \param[in] active_orb_i, active_orb_j The lists of the active orbitals for each MO-AO data set
  This function returns overlap matrix of atomic orbitals with different time step
  like <MO(t)|MO(t+dt)>.

*/

  int Nbas_i = ao_i.size();
  int Nbas_j = ao_j.size();

  // Ci = Nbas_i x Norb_i 
  // Cj = Nbas_j x Norb_j
  int Norb_i = Ci.n_cols;
  int Norb_j = Cj.n_cols;

  int Norb_i_act = active_orb_i.size();  // the number of active orbitals in set i
  int Norb_j_act = active_orb_j.size();  // the number of active orbitals in set j


  // Smo = MATRIX(Norb_i_act, Norb_j_act) <- this is what expected
  if(Smo.n_rows != Norb_i_act){
    cout<<"Error in MO_overlap: The number of rows of the MO overlap matrix"<<Smo.n_rows
        <<" is not consistent with the number of <i| orbitals "<<Norb_i_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }
  if(Smo.n_cols != Norb_j_act){
    cout<<"Error in MO_overlap: The number of cols of the MO overlap matrix"<<Smo.n_cols
        <<" is not consistent with the number of |j> orbitals "<<Norb_j_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }



  // Allocate working memory
  int i,j;
  int is_derivs = 0;
  int is_normalize = 0;
  VECTOR dIdA, dIdB;
  int n_aux = 20;
  vector<double*> auxd(10);
  for(i=0;i<10;i++){ auxd[i] = new double[n_aux]; }

    
  CMATRIX* ci; ci = new CMATRIX(Nbas_i, Norb_i_act);
  CMATRIX* cj; cj = new CMATRIX(Nbas_j, Norb_j_act);
  CMATRIX* Sao; Sao = new CMATRIX(Nbas_i, Nbas_j);

  pop_cols(Ci, *ci, active_orb_i);
  pop_cols(Cj, *cj, active_orb_j);  


  // overlap matrix of S
  for(i=0;i<Nbas_i;i++){  // all orbitals
    for(j=0;j<Nbas_j;j++){

      double dist2 = (ao_i[i].primitives[0].R - ao_j[j].primitives[0].R).length2();       

      if(dist2<max_d2) {  
        double res = gaussian_overlap(ao_i[i], ao_j[j], is_normalize, is_derivs, dIdA, dIdB, auxd, n_aux);
        Sao->set(i,j, res, 0.0); 
      }
      else{  Sao->set(i,j,0.0, 0.0);  }


    }// for j
  }// for i


  Smo = (*ci).H() * (*Sao) * (*cj);   // ( Norb_i_act x Nbas_i ) x (Nbas_i x Nbas_j) x (Nbas_j x Norb_j_act) = Norb_i_act x Norb_j_act


  // Clean working memory
  for(i=0;i<10;i++){ delete [] auxd[i]; }  
  auxd.clear();

  delete Sao;
  delete ci;
  delete cj;


}



void MO_overlap(MATRIX& Smo, MATRIX& Ci, MATRIX& Cj,
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2){
/**
  Finds the keywords and their patterns and extracts the parameters
  \param[in] Ci, Cj : molecular coefficients at different time step.
  \param[in] active_orb_i, active_orb_j The lists of the active orbitals for each MO-AO data set
  This function returns overlap matrix of atomic orbitals with different time step
  like <MO(t)|MO(t+dt)>.

*/

  // Ci = Nbas_i x Norb_i 
  // Cj = Nbas_j x Norb_j
  int Nbas_i = Ci.n_rows;
  int Nbas_j = Cj.n_rows;

  if(Nbas_i != Nbas_j){
    cout<<"Error in MO_overlap: The number of rows of the MO-LCAO matrix Ci "<<Ci.n_rows
        <<" is not equal to the number of rows of the MO-LCAO matrix Cj "<<Cj.n_rows<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }


  int Norb_i = Ci.n_cols;
  int Norb_j = Cj.n_cols;

  int Norb_i_act = active_orb_i.size();  // the number of active orbitals in set i
  int Norb_j_act = active_orb_j.size();  // the number of active orbitals in set j


  // Smo = MATRIX(Norb_i_act, Norb_j_act) <- this is what expected
  if(Smo.n_rows != Norb_i_act){
    cout<<"Error in MO_overlap: The number of rows of the MO overlap matrix"<<Smo.n_rows
        <<" is not consistent with the number of <i| orbitals "<<Norb_i_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }
  if(Smo.n_cols != Norb_j_act){
    cout<<"Error in MO_overlap: The number of cols of the MO overlap matrix"<<Smo.n_cols
        <<" is not consistent with the number of |j> orbitals "<<Norb_j_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }


  MATRIX* ci; ci = new MATRIX(Nbas_i, Norb_i_act);
  MATRIX* cj; cj = new MATRIX(Nbas_j, Norb_j_act);

  pop_cols(Ci, *ci, active_orb_i);
  pop_cols(Cj, *cj, active_orb_j);  


  Smo = (*ci).T() * (*cj);   // ( Norb_i_act x Nbas_i ) x (Nbas_j x Norb_j_act) = Norb_i_act x Norb_j_act


  // Clean working memory
  delete ci;
  delete cj;


}

void MO_overlap(CMATRIX& Smo, CMATRIX& Ci, CMATRIX& Cj,
 vector<int>& active_orb_i, vector<int>& active_orb_j, double max_d2){
/**
  Finds the keywords and their patterns and extracts the parameters
  \param[in] Ci, Cj : molecular coefficients at different time step.
  \param[in] active_orb_i, active_orb_j The lists of the active orbitals for each MO-AO data set
  This function returns overlap matrix of atomic orbitals with different time step
  like <MO(t)|MO(t+dt)>.

*/

  // Ci = Nbas_i x Norb_i 
  // Cj = Nbas_j x Norb_j
  int Nbas_i = Ci.n_rows;
  int Nbas_j = Cj.n_rows;

  if(Nbas_i != Nbas_j){
    cout<<"Error in MO_overlap: The number of rows of the MO-LCAO matrix Ci "<<Ci.n_rows
        <<" is not equal to the number of rows of the MO-LCAO matrix Cj "<<Cj.n_rows<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }


  int Norb_i = Ci.n_cols;
  int Norb_j = Cj.n_cols;

  int Norb_i_act = active_orb_i.size();  // the number of active orbitals in set i
  int Norb_j_act = active_orb_j.size();  // the number of active orbitals in set j


  // Smo = MATRIX(Norb_i_act, Norb_j_act) <- this is what expected
  if(Smo.n_rows != Norb_i_act){
    cout<<"Error in MO_overlap: The number of rows of the MO overlap matrix"<<Smo.n_rows
        <<" is not consistent with the number of <i| orbitals "<<Norb_i_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }
  if(Smo.n_cols != Norb_j_act){
    cout<<"Error in MO_overlap: The number of cols of the MO overlap matrix"<<Smo.n_cols
        <<" is not consistent with the number of |j> orbitals "<<Norb_j_act<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }


  CMATRIX* ci; ci = new CMATRIX(Nbas_i, Norb_i_act);
  CMATRIX* cj; cj = new CMATRIX(Nbas_j, Norb_j_act);

  pop_cols(Ci, *ci, active_orb_i);
  pop_cols(Cj, *cj, active_orb_j);  


  Smo = (*ci).H() * (*cj);   // ( Norb_i_act x Nbas_i ) x (Nbas_j x Norb_j_act) = Norb_i_act x Norb_j_act


  // Clean working memory
  delete ci;
  delete cj;


}


complex<double> SD_overlap(SD& sd_i, SD& sd_j){
/**
  \brief This function computes the overlap of two Slater determinants sd_i and sd_j
  \param[in] sd_i, sd_j : Slater Determinants at possibly different time steps (or at the same one)
  The computed overlap value will be returned

  <SD(t)|SD(t')>. =  det | <psi_i(t) | psi_j(t') > |  

  This is based on the Lowdin formula. See, for instance:
  (1) Ryabinkin, I. G.; Nagesh, J.; Izmaylov, A. F. Fast Numerical Evaluation of Time-Derivative Non-Adiabatic Couplings for Mixed Quantum-Classical Methods. J. Phys. Chem. Lett. 2015, 6, 4200–4203.


*/

  if(sd_i.N_bas != sd_j.N_bas){
    cout<<"Error in SD_overlap: The number of basis functions in which MOs of the SD sd_i are expanded "<<sd_i.N_bas
        <<" is not equal to the number of basis functions in which MOs of the SD sd_j are expanded "<<sd_j.N_bas<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }
  if(sd_i.N!=sd_j.N){
    cout<<"Error in SD_overlap: The number of MOs included in the SD sd_i "<<sd_i.N
        <<" is not equal to the number of MOs included in the SD sd_j"<<sd_j.N<<endl;
    cout<<"Exiting..."; 
    exit(0);
  }

  // Compute the matrix of 1-el MO overlaps
  CMATRIX* Smo; Smo = new CMATRIX(sd_i.N,sd_j.N); // square matrix
  *Smo = sd_i.get().H() * sd_j.get();   // ( sd_i.N x sd_i.N_bas ) x (sd_j.N_bas x sd_j.N) = sd_i.N x sd_j.N

  // Now, apply spin considerations
  for(int i=0;i<sd_i.N;i++){
    for(int j=0;j<sd_j.N;j++){
      if(sd_i.spin[i]!=sd_j.spin[j]){ Smo->set(i,j, 0.0, 0.0);   }
    }// for j
  }// for i

  // Compute the determinant
  double nrm = sd_i.normalization_factor();
  complex<double> res = libmeigen::det(*Smo); 


  // Clean working memory
  delete Smo;

  return res;

}

CMATRIX SD_overlap(vector<SD>& sd_i, vector<SD>& sd_j){
/**
  \brief This function computes the matrix of the SD overlaps from two data sets (e.g. fragments or timesteps)
  \param[in] sd_i, sd_j : Are the lists of SDs belonging to each of the two data sets
  The computed matrix of overlaps value will be returned
*/

  int Ni = sd_i.size();
  int Nj = sd_j.size();

  CMATRIX* SD_ovlp; SD_ovlp = new CMATRIX(Ni,Nj);

  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      SD_ovlp->set(i,j, SD_overlap(sd_i[i], sd_j[j]) ); 
    }// for j
  }// for i

  return *SD_ovlp;

}

void SD_overlap(CMATRIX& SD_ovlp, vector<SD>& sd_i, vector<SD>& sd_j){
/**
  \brief This function computes the matrix of the SD overlaps from two data sets (e.g. fragments or timesteps)
  \param[out] SD_ovlp The matrix storing the results that is to be updated
  \param[in] sd_i, sd_j : Are the lists of SDs belonging to each of the two data sets  
*/

  int Ni = sd_i.size();
  int Nj = sd_j.size();

  if(SD_ovlp.n_rows!=Ni){
    std::cout<<"Error in SD_overlap : the # of rows of the output matrix, SD_ovlp ( "<<SD_ovlp.n_rows
             <<" ) is not equal to the number of Slater Determinants in the first (left) set ( "<<Ni<<" )\n";
    exit(0);
  }
  if(SD_ovlp.n_cols!=Nj){
    std::cout<<"Error in SD_overlap : the # of cols of the output matrix, SD_ovlp ( "<<SD_ovlp.n_cols
             <<" ) is not equal to the number of Slater Determinants in the second (right) set ( "<<Nj<<" )\n";
    exit(0);
  }


  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      SD_ovlp.set(i,j, SD_overlap(sd_i[i], sd_j[j]) ); 
    }// for j
  }// for i


}





}//namespace libbasis
}//namespace liblibra
