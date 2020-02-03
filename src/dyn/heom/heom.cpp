/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file heom.cpp
  \brief The file implements the procedures for the HEOM calculations

  Based on the Fortran 90 code of Amber Jain & Joe Subotnik:  https://github.com/subotnikgroup/HEOM_Amber    
*/

#include "heom.h"
#include "../../math_linalg/liblinalg.h"
#include "libheom.h"

namespace liblibra{
namespace libdyn{
namespace libheom{

using namespace liblinalg;



int compute_nn_tot(int nquant, int KK, int LL){
/**
  \brief Compute the total number of density matrices in the EOM hierarchy

  This number is given by:   (LL+nquant*KK)! / ( LL! * (nquant*KK)! )

  of course, we don't want to evaluate the factorial explicitly

  \param[in] L Level of the HEOM hierarchy
  \param[in] nquant (often denoted M) the number of system-bath coupling terms, in the case
  where each state is coupled to the bath intependently, this number is also equal to the number
  of quantum state, hence the current name
  \param[in] KK the number of the Matsubara terms


*/

  int i;
  double tmp = 1.0;

  for(int i=1; i<=nquant*(KK+1); i++){
    tmp *= (LL+i)/((float)i);
  }

  return (int)tmp;
}


vector<int> allocate_1D(int sz1){

  vector<int> res(sz1, 0);

  return res;

}

vector< vector<int> > allocate_2D(int sz1, int sz2){

  vector< vector<int> > res(sz1, vector<int>(sz2, 0) );

  return res;
}

vector< vector< vector<int> > > allocate_3D(int sz1, int sz2, int sz3){

  vector< vector< vector<int> > > res(sz1, vector< vector<int> >(sz2,  vector<int>(sz3, 0) ) );

  return res;
}





void compute_nn(int nquant, int KK, int LL, 
                vector<int>& map_sum, 
                vector< vector< vector<int> > >& nn){
/**
  \brief ?? 
  \param[in] nquant
  \param[in] KK
  \param[in] LL
  \param[in/out] map_sum  is of size LL+1:  [0, ... , LL]
  \param[in/out] map_nneg, map_nplus, nn  are of size nquant+1 x KK+1 x LL+1:  [0,... nquant] x [0, KK] x [0, ... , LL] each

*/

  // Alternative recursive method
  int n_beg, n_end;
  n_beg = n_end = 0;
  
  for(int n=0; n<=LL; n++){

    compute_nn_sum_L(nquant, KK, n, n_beg, n_end, nn);
    map_sum[n] = n_beg;

    cout<<"n = "<<n<<" n_beg = "<<n_beg<<" n_end = "<<n_end<<endl;
  }
  cout<<" n_end = "<<n_end<<endl;

}


void compute_nn_sum_L(int nquant, int KK, int L, int& n_beg, int& n_end, 
                      vector< vector< vector<int> > >& nn){
/**
  \brief ?? 
  \param[in] L Level of the HEOM hierarchy
  \param[in/out] n_beg ?? at input, state/end point of entries for sum {L-1}; at output, for entries for sum {L}
  \param[out] nn (IMATRIX(nquant, KK+1) ) the vectors of the indices
  \param[in] nquant the number of quantum states
  \param[in] KK the number of the Matsubara terms
  \param[in/out] nn  is of size nquant+1 x KK+1 x LL+1:  [0,... nquant] x [0, KK] x [0, ... , LL] each

*/
  int tot_L, i, j, m, k, a,b;
  int flag_new; 

  vector< vector<int> > nvec(nquant+1, vector<int>(KK+1, 0) );

  if(L==0){

    for(a=0;a<=nquant;a++){  for(b=0;b<=KK;b++){  nn[a][b][1] = 0;   } }
    n_beg=1;
    n_end=1;
  }
  else{

    tot_L=n_end;

    for(i=n_beg; i<=n_end; i++){
      for(m=1; m<=nquant; m++){

        for(k=0; k<=KK; k++){

          for(a=0;a<=nquant;a++){  for(b=0;b<=KK;b++){ nvec[a][b] = nn[a][b][i];   } }          
          nvec[m][k] +=1;


          flag_new=0;
          if(tot_L > n_end){
            for(j=n_end+1; j<=tot_L; j++){

              int same = 1;
              for(a=1;a<=nquant;a++){  for(b=0;b<=KK;b++){  if(nvec[a][b] != nn[a][b][j]){  same = 0; }   } }          
              if(same==1){ flag_new=1; }

            }// for j
          }

          if(flag_new==0){
            tot_L += 1;

            for(a=0;a<=nquant;a++){  for(b=0;b<=KK;b++){ nn[a][b][tot_L] = nn[a][b][i];   } }          
            nn[m][k][tot_L] = nn[m][k][i] + 1; 

          }

        }// for k

      }// for m

    }// for i

    n_beg = n_end+1;
    n_end = tot_L;
    
  }// else


}




vector< vector<int> > index_int2vec(vector< vector< vector<int> > >& nn, int n, int nquant, int KK){
/**
  \brief Computes the correspondence between the many-integers representation of 
  the auxiliary density matrix and their single-integer index

  ** n is input, nvec is output **

  \param[in] nn  The set of all many-integer indices
  \param[in] n  The single-integer index

  
*/
  vector< vector<int> > nvec(nquant+1, vector<int>(KK+1, 0) );

  int a,b;

  for(a=0; a<=nquant; a++){
    for(b=0; b<=KK; b++){
      nvec[a][b] = nn[a][b][n];
    }
  }

  return nvec;

}


int sum2D(vector< vector<int> >& nvec){

  int sz1 = nvec.size();
  int sz2 = nvec[0].size();

  int i,j;
  int summ = 0;

  for(i=0; i<sz1; i++){
    for(j=0; j<sz2; j++){
      summ += nvec[i][j];
    }
  }

  return summ;
}

int index_vec2int(vector< vector< vector<int> > >& nn, vector< vector<int> >& nvec, int LL){
/**
  \brief Computes the correspondence between the many-integers representation of 
  the auxiliary density matrix and their single-integer index

  ** nvec is input, n is output **

  \param[in] nn  The set of all many-integer indices
  \param[in/out] n  The single-integer index
  \param[in/out] nvec  The many-integers index
  \param[in] iflag The direction of the transformation

*/

  int i,j;

  int sz1 = nn.size();
  int sz2 = nn[0].size();
  int nn_tot = nn[0][0].size() - 1;

  int n = 0;

  int l_sum = sum2D(nvec);

  if(l_sum<=LL){

    int found = 0;

    for(int m=1; m<=nn_tot && found==0; m++){ 
    
      int are_equal = 1;
    
      for(i=0; i<sz1 && are_equal; i++){
        for(j=0; j<sz2 && are_equal; j++){
          if(nvec[i][j] != nn[i][j][m]){  are_equal = 0; }
        }
      }
    
      if(are_equal==1){  n = m;  found = 1;   }
    
    }// for m

  }// if l_sum

  return n;
        
}




void compute_map(int nquant, int KK, int LL, 
                 vector< vector< vector<int> > >& nn,
                 vector< vector< vector<int> > >& map_nplus,
                 vector< vector< vector<int> > >& map_nneg){
/**
  \brief Compute the mapping of the overall index of the density matrix to the 
  3 integers that characterize the density matrix. In fact, this function return
  the enumerated list of such 3-integer vectors

  \param[in] L Level of the HEOM hierarchy
  \param[in] nquant (often denoted M) the number of system-bath coupling terms, in the case
  where each state is coupled to the bath intependently, this number is also equal to the number
  of quantum state, hence the current name
  \param[in] KK the number of the Matsubara terms

  \param[in/out] map_nneg, map_nplus, nn  are of size nquant+1 x KK+1 x LL+1:  [0,... nquant] x [0, KK] x [0, ... , LL] each


*/

  int n,m,k;
  int nplus, nneg;
  vector< vector<int> > nvec(nquant+1, vector<int>(KK+1, 0) );
  vector< vector<int> > nvec_plus(nquant+1, vector<int>(KK+1, 0) );
  vector< vector<int> > nvec_neg(nquant+1, vector<int>(KK+1, 0) );


  int nn_tot = compute_nn_tot(nquant, KK, LL);

//  map_nplus = vector< vector< vector<int> > >(nquant+1, vector< vector<int> >(KK+1, vector<int>(nn_tot+1, 0)) );
//  map_nneg = vector< vector< vector<int> > >(nquant+1, vector< vector<int> >(KK+1, vector<int>(nn_tot+1, 0)) );

  for(n=1; n<=nn_tot; n++){

    m=n;
    nvec = index_int2vec(nn, m, nquant, KK);

    for(m=1; m<=nquant; m++){
      for(k=0; k<=KK; k++){

        nvec_plus = nvec;
        nvec_plus[m][k] = nvec_plus[m][k] + 1;

        map_nplus[m][k][n] = index_vec2int(nn, nvec_plus, LL);
       
        if(nvec[m][k] > 0){
          nvec_neg = nvec;
          nvec_neg[m][k] = nvec_neg[m][k] - 1;

          map_nneg[m][k][n] = index_vec2int(nn, nvec_neg, LL);
        }

      }// for k
    }// for m

  }// for n

}


vector<int> filter(vector<CMATRIX>& rho, double tolerance){
/**
   Shi, Chen, Nan, Xu, Yan, JCP 130, 084105 (2009)
*/

  int n;
  int nn_tot = rho.size()-1; 
  vector<int> zero(nn_tot+1, 0);

  for(int n=1; n<=nn_tot; n++){

    if( abs(rho[n].max_elt() ) < tolerance ){
      rho[n] = 0.0;  zero[n] = 1;
    }

  }// for n

  return zero;
}



CMATRIX compute_deriv_n(int n, vector<CMATRIX>& rho, CMATRIX& Ham, 
        double eta, double temperature,
        vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
        vector< vector< vector<int> > >& nn, int KK, vector<int>& zero,
        vector< vector< vector<int> > >& map_nplus, vector< vector< vector<int> > >& map_nneg){
/**

   Compute drho_n/dt = ...     !! Eq. 15, JCP 131, 094502 (2009)

*/

  int nn_tot = rho.size();
  int nquant = rho[0].n_cols;

  int m,k,nplus,nminus;
  CMATRIX drho_n_dt(nquant, nquant);
  CMATRIX commut(nquant, nquant);
  CMATRIX term1(nquant, nquant);
  CMATRIX term2(nquant, nquant);
  CMATRIX mat_tmp(nquant, nquant);
  CMATRIX mat_tmp2(nquant, nquant);

  vector< vector<int> > nvec(nquant, vector<int>(KK+1, 0) );
  vector< vector<int> > nvec_plus(nquant, vector<int>(KK, 0) );
  vector< vector<int> > nvec_neg(nquant, vector<int>(KK, 0) );

  complex<double> iota(0.0, 1.0);
  double kB = boltzmann/hartree;
 

  m=n;
  nvec = index_int2vec(nn, m, nquant, KK);

  drho_n_dt = -iota * ( Ham * rho[n] - rho[n] * Ham ); 


  if(zero[n]==0){  // matrix at n is not filtered out

      for(m=1; m<=nquant; m++){
        for(k=0; k<=KK; k++){
          drho_n_dt -= nvec[m][k] * gamma_matsubara[k] * rho[n];
        }// for k
      } // for m


      complex<double> matsubara_sum(0.0, 0.0);
      for(k=0; k<=KK; k++){
        matsubara_sum += c_matsubara[k]/gamma_matsubara[k];
      }

      double pref = eta * kB * temperature/gamma_matsubara[0] - std::real(matsubara_sum);

      for(m=1; m<=nquant; m++){

        mat_tmp=0.0;
        //mat_tmp.set(m,m, complex<double>(1.0, 0.0));  // matrix |m><m|
        mat_tmp.set(m-1,m-1, complex<double>(1.0, 0.0));  // matrix |m><m|

        commut = mat_tmp * rho[n] - rho[n] * mat_tmp;
        commut = mat_tmp * commut - commut * mat_tmp;

        drho_n_dt -= pref *commut;

      } // for m
  }// if 


  for(m=1; m<=nquant; m++){

      mat_tmp=0.0;
      //mat_tmp.set(m,m, complex<double>(1.0, 0.0));  // matrix |m><m|
      mat_tmp.set(m-1,m-1, complex<double>(1.0, 0.0));  // matrix |m><m|

      mat_tmp2=0.0;

      for(k=0; k<=KK; k++){

        nplus = map_nplus[m][k][n];

        if(nplus> 0 && nplus<=nn_tot){

          if(zero[nplus]==0){
            drho_n_dt -=  iota*( mat_tmp * rho[nplus] - rho[nplus] * mat_tmp) * sqrt((nvec[m][k]+1.0)*abs(c_matsubara[k]));
          }

        }// if 
      }// for k
  } // for m


  for(m=1; m<=nquant; m++){

      mat_tmp=0.0;
      //mat_tmp.set(m,m, complex<double>(1.0, 0.0));  // matrix |m><m|
      mat_tmp.set(m-1,m-1, complex<double>(1.0, 0.0));  // matrix |m><m|

      for(k=0; k<=KK; k++){

        nminus = map_nneg[m][k][n];
        if(nminus>0 && nminus<=nn_tot){

          if(zero[nminus]==0){
            term1 = c_matsubara[k] * mat_tmp * rho[nminus];
            term2 = std::conj(c_matsubara[k]) * rho[nminus] * mat_tmp;

            drho_n_dt -= iota * (term1 - term2) * nvec[m][k] * sqrt(nvec[m][k]/std::real(c_matsubara[k]) );
          }

        }// if 
      }// for k 
  }// for m


  return drho_n_dt;

}



void unpack_rho(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO){
/**
*/

    int i, n;
    int nquant = RHO.n_cols;
    int nn_tot = int(RHO.n_rows / nquant) - 1;

    if(rho_unpacked.size() != (nn_tot+1) ) {   cout<<"ERROR in unpack_rho: \n"; exit(0); }
    

    vector<int> sub_x(nquant, 0);
    vector<int> sub_y(nquant, 0);
    for(i = 0; i<nquant; i++){    sub_y[i] = i;  }

    
    for(n = 0; n<=nn_tot; n++){

      for(i = 0; i<nquant; i++){    sub_x[i] = n*nquant + i;  }
      pop_submatrix(RHO, rho_unpacked[n], sub_x, sub_y );
        
    }

}

void pack_rho(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO){

    int i, n;
    int nquant = RHO.n_cols;
    int nn_tot = int(RHO.n_rows / nquant) - 1;

    if(rho_unpacked.size() != (nn_tot+1) ) {   cout<<"ERROR in pack_rho: \n"; exit(0); }
    

    vector<int> sub_x(nquant, 0);
    vector<int> sub_y(nquant, 0);
    for(i = 0; i<nquant; i++){    sub_y[i] = i;  }

    
    for(n = 0; n<=nn_tot; n++){

      for(i = 0; i<nquant; i++){    sub_x[i] = n*nquant + i;  }
      push_submatrix(RHO, rho_unpacked[n], sub_x, sub_y );
        
    }
}


CMATRIX compute_heom_derivatives(CMATRIX& RHO, bp::dict prms){
    /**

    RHO - CMATRIX((nn_tot+1)*nquant, nquant) - all the density matrices stacked together in a super-column
    prms - Python dict - control parameters
    
    */

    int i, n;
    int nquant = RHO.n_cols;
    int nn_tot = int(RHO.n_rows / nquant) - 1;


    CMATRIX Ham(nquant, nquant);
    double eta = 0.0;
    double temperature = 300.0;
    int KK;

    vector<double> gamma_matsubara;
    vector< complex<double> > c_matsubara;
    intList3 nn, map_nplus, map_nneg;
    intList zero;

    
    std::string key;
    for(i=0;i<len(prms.values());i++){

      key = extract<std::string>(prms.keys()[i]);

      if(key=="Ham"){  Ham = extract<CMATRIX>(prms.values()[i]); }
      if(key=="eta"){  eta = extract<double>(prms.values()[i]); }
      if(key=="temperature"){  temperature = extract<double>(prms.values()[i]); }
      if(key=="KK"){  KK = extract< int >(prms.values()[i]); }

      if(key=="gamma_matsubara"){  gamma_matsubara = extract< doubleList >(prms.values()[i]); }
      if(key=="c_matsubara"){  c_matsubara = extract< complexList >(prms.values()[i]); }

      if(key=="nn"){  nn = extract< intList3 >(prms.values()[i]); }
      if(key=="map_nplus"){  map_nplus = extract< intList3 >(prms.values()[i]); }
      if(key=="map_nneg"){  map_nneg = extract< intList3 >(prms.values()[i]); }
      if(key=="zero"){  zero = extract< intList >(prms.values()[i]); }

    }




    vector<int> sub_x(nquant, 0);
    vector<int> sub_y(nquant, 0);

    for(i = 0; i<nquant; i++){    sub_y[i] = i;  }

    CMATRIX dRHO((nn_tot+1)*nquant, nquant);
    CMATRIX drho(nquant, nquant);
    vector<CMATRIX> rho_unpacked(nn_tot+1, CMATRIX(nquant, nquant));
    vector<CMATRIX> drho_unpacked(nn_tot+1, CMATRIX(nquant, nquant));

    unpack_rho(rho_unpacked, RHO);


    vector< vector<int> > sub_x_all(nn_tot+1, vector<int>(nquant, 0));
    for(n=1; n<=nn_tot; n++){
        for(i = 0; i<nquant; i++){    sub_x_all[n][i] = n*nquant + i;  }
    }


    #pragma omp parallel for
    for(n=1; n<=nn_tot; n++){

        //for(i = 0; i<nquant; i++){    sub_x[i] = n*nquant + i;  }
        
        drho_unpacked[n] = compute_deriv_n(n, rho_unpacked, Ham, eta, temperature, 
          gamma_matsubara, c_matsubara,  nn, KK, zero, map_nplus,  map_nneg);

        //push_submatrix(dRHO, drho, sub_x_all[n], sub_y );

    }

    for(n=1; n<=nn_tot; n++){
        push_submatrix(dRHO, drho_unpacked[n], sub_x_all[n], sub_y );
    }


    return dRHO;

}


void setup_bath(bp::dict prms, vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara){

  int k, i;

  complex<double> one(1.0, 0.0);
  complex<double> iota(0.0, 1.0);

  int KK = 0;
  double eta = 0.0;
  double gamma = 0.0;
  double temperature = 300.0;

  std::string key;
  for(i=0;i<len(prms.values());i++){

    key = extract<std::string>(prms.keys()[i]);

    if(key=="KK"){  KK = extract< int >(prms.values()[i]); }
    if(key=="eta"){  eta = extract<double>(prms.values()[i]); }
    if(key=="gamma"){  gamma = extract<double>(prms.values()[i]); }
    if(key=="temperature"){  temperature = extract<double>(prms.values()[i]); }


  }

  double kB = boltzmann/hartree;
  double kT = kB * temperature;

  c_matsubara = vector< complex<double> >(KK+1, complex<double>(0.0, 0.0) );
  gamma_matsubara = vector<double>(KK+1, 0.0 );

  gamma_matsubara[0] = gamma ;
  c_matsubara[0] = 0.5*eta*gamma * ( 1.0/tan( 0.5 * gamma/kT )*one - iota );

  for(k=1; k<=KK; k++){

    gamma_matsubara[k] = 2.0*k*M_PI*kT;

    double g = gamma_matsubara[k];
    c_matsubara[k] = 2*eta*gamma*kT * g/(g*g - gamma*gamma);
  }


}



}// namespace libheom
}// namespace libdyn
}// liblibra

