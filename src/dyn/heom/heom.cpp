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


vector< vector<int> > gen_next_level(vector<int>& parent){
    /**
    parent = ( n10, n11, ... , n1K    ...     nM0, nM1, ... , nMK)

    This function generates the vectors of integers our of a given vector
    by incrementing every element of the original (parent) vector by 1
 
    The result is the list of the resulting children vectors 

    */

    int d = parent.size();  // shall be M * (K+1)

    vector< vector<int> > children;

    for(int i=0;i<d; i++){
        children.push_back( parent );
        children[i][i] += 1;
    }

    return children;
}


vector< vector<int> > gen_next_level2(vector< vector<int> >& parents){
    /**
    This function is similar to gen_next_level, except it takes a list of
    parent vectors of integers and generates the list of the children for 
    all of the parent inputs.
    */

    vector< vector<int> > next_level;

    int d = parents.size();
    for(int i=0;i<d; i++){
    
        vector< vector<int> > children = gen_next_level(parents[i]);

        int num_children = children.size();

        for(int j=0;j<num_children; j++){
            next_level.push_back(children[j]);
        }
    }

    return next_level;
}


int is_equal(vector<int>& vec1, vector<int>& vec2){

  int res = 1; 
  int sz1 = vec1.size();
  int sz2 = vec2.size();

  if(sz1!=sz2){ res = 0; }
  else{
    for(int i=0; i<sz1; i++){
      if(vec1[i] != vec2[i]){ res = 0; break; }
    }
  }
  return res;
}

int is_included(vector<int>& vec1, vector<vector<int> >& vec){

  int res = 0;
  for(int i=0; i<vec.size(); i++){
    if(is_equal(vec1, vec[i])){ res = 1; break; }
  }
  return res;
}

void gen_hierarchy(int d, int max_tier, int verbosity,
                   vector< vector<int> >& all_vectors, 
                   vector< vector<int> >& vec_plus, 
                   vector< vector<int> >& vec_minus){
    /**

    d  - is the size of the simplex, in HEOM it will be taken as nquant * (KK+1)

    max_tier - is the maximal tier of the vectors to be generated, the tier is
    devined by the sum of the vector elements

    all_vectors - a list of the d-dimensional vectors of ints 

    vec_plus : vec_plus[n][k] - index of the vector for which k-th element is +1 of that of the n-th vector

    vec_minus : vec_minus[n][k] - index of the vector for which k-th element is -1 of that of the n-th vector

    The vec_minus and vec_plus contain -1 for the situation when there is no such a vector included in 
    the current hierarchy structure
 
    */

    int i, j, k;

 
    vector< vector<int> > all_coordinates; //  for each vector - the coordinates are (L, i)
    vector<int> tier_nums; // the number of the nodes for the tier up to a given one

    vector< vector<int> > parents(1, vector<int>(d, 0));

    for(int tier=0; tier<=max_tier; tier++){

        int iparent = 0;
        int num_parents = parents.size();

        for(i=0; i<num_parents; i++){
            if( !is_included(parents[i], all_vectors) ){
                all_vectors.push_back( parents[i] );
                vector<int> coord(2,0);
                coord[0] = tier;
                coord[1] = iparent;
                all_coordinates.push_back( coord );
                iparent += 1;
            }// not yet included
        }// for i

        tier_nums.push_back( all_vectors.size() ) ;

        vector< vector<int> > new_parents = gen_next_level2(parents);
        parents.resize(new_parents.size());
        parents = new_parents;

    }// for tier

    //#==============================================


    int num_nodes = all_vectors.size();
    vec_plus = vector< vector<int> >(num_nodes, vector<int>(d, -1)); 
    vec_minus = vector< vector<int> >(num_nodes, vector<int>(d, -1)); 


    for(i=0; i<num_nodes; i++){
        int L = all_coordinates[i][0];
        int ipos = all_coordinates[i][1];

        for(k=0; k<d; k++){
            vector<int> np = all_vectors[i];
            np[k] += 1;

            int max_range = max_tier;
            if(L<max_tier){
                max_range = tier_nums[L+1];
            }

            for(j = tier_nums[L]; j < max_range; j++){
                if(np == all_vectors[j]){
                    vec_plus[i][k] = j;
                }
            }
        }// for k

        for(k=0; k<d; k++){
            vector<int> nm = all_vectors[i];
            nm[k] -= 1;

            int min_range = 0;
            if(L>=2){
                min_range = tier_nums[L-2];
            }
            for(j=min_range; j<tier_nums[L-1]; j++){
                if(nm == all_vectors[j]){
                    vec_minus[i][k] = j;
                }
            }
        }// for k
    }// for i


    if(verbosity>0){
        cout<<"Number of nodes = "<<num_nodes<<endl;
        cout<<"Tier nums: \n";
        for(i=0; i<tier_nums.size(); i++){
            cout<<"Tier "<<i<<" nums = "<<tier_nums[i]<<endl;
        }
    }

    if(verbosity>1){
        for(i=0; i<num_nodes; i++){

            cout<<"index= "<<i<<" ";
    
            cout<<" vector=["; for(j=0;j<d;j++){ cout<<all_vectors[i][j]<<","; } cout<<"] ";

            cout<<" coordinates=["<<all_coordinates[i][0]<<","<<all_coordinates[i][1]<<"] ";

            cout<<" vec_plus=["; for(j=0;j<d;j++){ cout<<vec_plus[i][j]<<","; } cout<<"] ";

            cout<<" vec_minus=["; for(j=0;j<d;j++){ cout<<vec_minus[i][j]<<","; } cout<<"] ";

            cout<<"\n";
        }
    }

}


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


vector<CMATRIX> initialize_el_phonon_couplings(int nquant){
/**
  Compute the matrices  Q_m = |m><m| - the simplest type of the coupling of electrons to phonons

*/

  vector<CMATRIX> projectors(nquant+1, CMATRIX(nquant, nquant));

  for(int m=1; m<=nquant; m++){

    projectors[m] = 0.0;
    projectors[m].set(m-1,m-1, complex<double>(1.0, 0.0));  // matrix |m><m|
  }

  return projectors;

}


complex<double> compute_matsubara_sum(vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara, int KK){
/**
  Compute the sum of the c/gamma terms over all Matsubara frequencies
*/

  complex<double> matsubara_sum(0.0, 0.0);
  for(int k=0; k<=KK; k++){
    matsubara_sum += c_matsubara[k]/gamma_matsubara[k];
  }

  return matsubara_sum;

}


CMATRIX compute_deriv_n(int n, vector<CMATRIX>& rho, CMATRIX& Ham, vector<CMATRIX>& el_phon_coupl,
        double eta, double temperature,
        vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
        int do_truncate, int do_scale,
        vector< vector< vector<int> > >& nn, int KK, vector<int>& zero,
        vector< vector< vector<int> > >& map_nplus, vector< vector< vector<int> > >& map_nneg        
        ){
/**

   Compute drho_n/dt = ...     !! Eq. 15, JCP 131, 094502 (2009)

   Parameters:
    n : the order of the auxiliary matrix
    rho : density matrices (including the auxiliary ones) of the system
    Ham : system's Hamiltonian 
    projectors : the matrices that describe how each electronic phonon is coupled to various electronic states
       in the simplest picture, when a phonon k is coupled to electronic state m, the matrix projectors[k] will 
       contain 1.0 at the element (m,m) and zeroes everywhere else. You can of course define more general situations
       when one phonon is coupled to many states (and perhaps to their coherences)
    eta : reorganization energy [a.u.]
    temperature : bath temperature [ K ] 
    gamma_matsubara : Matsubara frequencies
    c_matsubara : expansion coefficients of the autocorrelation function in the Matsubara terms
    do_truncate : a flag to control the Ihizaki-Tanimura scheme for truncation -  JPSJ 74 3131, 2005
      0 - don't truncate 
      1 - do truncate 
    do_scale : a flag to control the scaled HEOM approach from JCP 130, 084105 (2009)
      0 - don't use scaling
      1 - do use scaling
    nn : the mapping of the variables
    KK : the number of Matsubara terms
    zero : mask to discard some of the matrices
    map_nplus : mapping of the rho_n+ matrices
    map_minus : mapping of the rho_n- matrices
    
    
   

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

  double scaling_factor = 1.0;
  complex<double> iota(0.0, 1.0);
  double kB = boltzmann/hartree;
 

  m=n;
  nvec = index_int2vec(nn, m, nquant, KK);

  //=============== Liouvillian =====================
  drho_n_dt = -iota * ( Ham * rho[n] - rho[n] * Ham ); 


  if(zero[n]==0){  // matrix at n is not filtered out

    //============= Friction ======================
    for(m=1; m<=nquant; m++){
      for(k=0; k<=KK; k++){
        drho_n_dt -= nvec[m][k] * gamma_matsubara[k] * rho[n];
      }// for k
    } // for m

       
    if(do_truncate==1){
      // Ihizaki-Tanimura scheme for truncation
      // JPSJ 74 3131, 2005
      complex<double> matsubara_sum = compute_matsubara_sum(gamma_matsubara, c_matsubara, KK);
      double pref = eta * kB * temperature/gamma_matsubara[0] - std::real(matsubara_sum);
      
      for(m=1; m<=nquant; m++){         
        commut = el_phon_coupl[m] * rho[n] - rho[n] * el_phon_coupl[m];
        commut = el_phon_coupl[m] * commut - commut * el_phon_coupl[m];
        drho_n_dt -= pref *commut;
      }

    }// if 
  }// zero[n] ==0

  


  for(m=1; m<=nquant; m++){
    for(k=0; k<=KK; k++){
      nplus = map_nplus[m][k][n];

      if(nplus > 0 && nplus <= nn_tot){
        if(zero[nplus]==0){

          scaling_factor = 1.0;          
          if(do_scale==1){
            // Scaled HEOM approach from JCP 130, 084105 (2009)
            scaling_factor = sqrt((nvec[m][k]+1.0)*abs(c_matsubara[k]));
          }

          drho_n_dt -=  iota*( el_phon_coupl[m] * rho[nplus] - rho[nplus] * el_phon_coupl[m]) * scaling_factor;

        }// zero[nplus] == 0
      }// if 

    }// for k
  } // for m


  for(m=1; m<=nquant; m++){
    for(k=0; k<=KK; k++){
      nminus = map_nneg[m][k][n];

      if(nminus > 0 && nminus <= nn_tot){

        if(zero[nminus]==0){

          term1 = c_matsubara[k] * el_phon_coupl[m] * rho[nminus];
          term2 = std::conj(c_matsubara[k]) * rho[nminus] * el_phon_coupl[m];


          scaling_factor = nvec[m][k];
          if(do_scale==1){
            // Scaled HEOM approach from JCP 130, 084105 (2009)
            scaling_factor =  sqrt(nvec[m][k]/abs(c_matsubara[k]));
          } 

           drho_n_dt -= iota * (term1 - term2) * scaling_factor;

        }// if
      }// if 
    }// for k 
  }// for m


  return drho_n_dt;

}




CMATRIX compute_deriv_n_new(int n, vector<CMATRIX>& rho, CMATRIX& Ham, vector<CMATRIX>& el_phon_coupl,
        double eta, double temperature,
        vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
        int do_truncate, int do_scale,
        vector<int>& nonzero,
        vector< vector<int> >& nvectors, 
        vector< vector<int> >& vec_plus, vector< vector<int> >& vec_minus        
        ){
/**

   Compute drho_n/dt = ...     !! Eq. 15, JCP 131, 094502 (2009)

   Parameters:
    n : the index of the ADOs
    rho : density matrices (including the auxiliary ones) of the system
    Ham : system's Hamiltonian 
    projectors : the matrices that describe how each electronic phonon is coupled to various electronic states
       in the simplest picture, when a phonon k is coupled to electronic state m, the matrix projectors[k] will 
       contain 1.0 at the element (m,m) and zeroes everywhere else. You can of course define more general situations
       when one phonon is coupled to many states (and perhaps to their coherences)
    eta : reorganization energy [a.u.]
    temperature : bath temperature [ K ] 
    gamma_matsubara : Matsubara frequencies
    c_matsubara : expansion coefficients of the autocorrelation function in the Matsubara terms
    do_truncate : a flag to control the Ihizaki-Tanimura scheme for truncation -  JPSJ 74 3131, 2005
      0 - don't truncate 
      1 - do truncate 
    do_scale : a flag to control the scaled HEOM approach from JCP 130, 084105 (2009)
      0 - don't use scaling
      1 - do use scaling

    adm_list : the list of the indices of all non-zero auxiliary density matrices (ADMs)
    nvectors : the multi-index vectors
    vec_plus : mapping of the rho_n+ matrices
    vec_minus : mapping of the rho_n- matrices

 
    
    
   

*/

  int nn_tot = rho.size();
  int nquant = rho[0].n_cols;

  int KK = gamma_matsubara.size() - 1; // KK+1 is the number of Matsubara freqs.

  int m,k,d,nplus,nminus;
  CMATRIX drho_n_dt(nquant, nquant);
  CMATRIX commut(nquant, nquant);
  CMATRIX term1(nquant, nquant);
  CMATRIX term2(nquant, nquant);
  CMATRIX mat_tmp(nquant, nquant);
  CMATRIX mat_tmp2(nquant, nquant);

  vector< vector<int> > nvec(nquant, vector<int>(KK+1, 0) );
  vector< vector<int> > nvec_plus(nquant, vector<int>(KK, 0) );
  vector< vector<int> > nvec_neg(nquant, vector<int>(KK, 0) );

  double scaling_factor = 1.0;
  complex<double> iota(0.0, 1.0);
  double kB = boltzmann/hartree;
 

  //=============== Liouvillian =====================
  drho_n_dt = -iota * ( Ham * rho[n] - rho[n] * Ham ); 


  //============= Friction ======================
  complex<double> pref; 

  for(k=0; k<=KK; k++){

    int sum_over_m = 0;
    for(m=0; m<nquant; m++){
      sum_over_m += nvectors[n][m*(KK+1) + k];
    }

    pref +=  sum_over_m * gamma_matsubara[k];
  }// for k

  drho_n_dt -=  pref * rho[n];

     
/*
  if(do_truncate==1){
    // Ihizaki-Tanimura scheme for truncation
    // JPSJ 74 3131, 2005
    complex<double> matsubara_sum = compute_matsubara_sum(gamma_matsubara, c_matsubara, KK);
    double pref = eta * kB * temperature/gamma_matsubara[0] - std::real(matsubara_sum);
    
    for(m=1; m<=nquant; m++){         
      commut = el_phon_coupl[m] * rho[n] - rho[n] * el_phon_coupl[m];
      commut = el_phon_coupl[m] * commut - commut * el_phon_coupl[m];
      drho_n_dt -= pref *commut;
    }

  }// if 
*/

  //===================== rho_n_plus terms ================

  for(m=0; m<nquant; m++){

    if(do_scale==1){

      // Scaled HEOM approach from JCP 130, 084105 (2009)
      // scaling_factor = sqrt((nvec[m][k]+1.0)*abs(c_matsubara[k]));

    }
 
    else{

      CMATRIX sum_rho_over_k(nquant, nquant);

      for(k=0; k<=KK; k++){
        nplus = vec_plus[n][m*(KK+1) + k];

        if(nplus!=-1 && nonzero[nplus]){
          sum_rho_over_k  += rho[nplus];
        }

      }// for k

      drho_n_dt -=  iota*( el_phon_coupl[m] * sum_rho_over_k - sum_rho_over_k * el_phon_coupl[m]);

    }// else - no scaling

  }// for m


  //===================== rho_n_minus terms ================

  for(m=0; m<nquant; m++){

    if(do_scale==1){

      // Scaled HEOM approach from JCP 130, 084105 (2009)
      // scaling_factor =  sqrt(nvec[m][k]/abs(c_matsubara[k]));

    }
 
    else{

      CMATRIX sum1_c_rho_over_k(nquant, nquant);
      CMATRIX sum2_c_rho_over_k(nquant, nquant);

      for(k=0; k<=KK; k++){
        nminus = vec_minus[n][m*(KK+1) + k];

        if(nminus!=-1 && nonzero[nminus]){
          sum1_c_rho_over_k  += (c_matsubara[k] * double(nvectors[n][m*(KK+1) + k])) * rho[nminus];
          sum2_c_rho_over_k  += (std::conj(c_matsubara[k]) * double(nvectors[n][m*(KK+1) + k])) * rho[nminus];
        }

      }// for k

      drho_n_dt -=  iota*( el_phon_coupl[m] * sum1_c_rho_over_k - sum2_c_rho_over_k * el_phon_coupl[m]);

    }// else - no scaling

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

    int n_threads = 1;
    int i, n;
    int nquant = RHO.n_cols;
    int nn_tot = int(RHO.n_rows / nquant) - 1;


    CMATRIX Ham(nquant, nquant);
    double eta = 0.0;
    double temperature = 300.0;
    int KK;

    vector<double> gamma_matsubara;
    vector< complex<double> > c_matsubara;
    int do_truncate = 1;
    int do_scale = 1;
    intList3 nn, map_nplus, map_nneg;
    intList zero;

    vector<CMATRIX> el_phon_couplings; 
    int el_phon_couplings_status = 0; // whether we have read them

    
    std::string key;
    for(i=0;i<len(prms.values());i++){

      key = extract<std::string>(prms.keys()[i]);

      if(key=="Ham"){  Ham = extract<CMATRIX>(prms.values()[i]); }
      if(key=="eta"){  eta = extract<double>(prms.values()[i]); }
      if(key=="temperature"){  temperature = extract<double>(prms.values()[i]); }
      if(key=="KK"){  KK = extract< int >(prms.values()[i]); }

      if(key=="gamma_matsubara"){  gamma_matsubara = extract< doubleList >(prms.values()[i]); }
      if(key=="c_matsubara"){  c_matsubara = extract< complexList >(prms.values()[i]); }
      if(key=="do_truncate"){  do_truncate = extract< int >(prms.values()[i]); }
      if(key=="do_scale"){  do_scale = extract< int >(prms.values()[i]); }

      if(key=="el_phon_couplings"){  
        el_phon_couplings = extract< CMATRIXList >(prms.values()[i]); 
        el_phon_couplings_status = 1;
      }

      if(key=="nn"){  nn = extract< intList3 >(prms.values()[i]); }
      if(key=="map_nplus"){  map_nplus = extract< intList3 >(prms.values()[i]); }
      if(key=="map_nneg"){  map_nneg = extract< intList3 >(prms.values()[i]); }
      if(key=="zero"){  zero = extract< intList >(prms.values()[i]); }

      if(key=="num_threads"){ n_threads = extract<int>(prms.values()[i]); }

    }


    if(el_phon_couplings_status==0){
      el_phon_couplings = initialize_el_phonon_couplings(nquant);
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

    #pragma omp parallel num_threads(n_threads)
    {
        #pragma omp for
        for(n=1; n<=nn_tot; n++){
            drho_unpacked[n] = compute_deriv_n(n, rho_unpacked, Ham, el_phon_couplings, eta, temperature, 
                gamma_matsubara, c_matsubara,  do_truncate, do_scale, nn, KK, zero, map_nplus,  map_nneg);
        }
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

  gamma_matsubara[0] = gamma;
  c_matsubara[0] = 0.5*eta*gamma * ( 1.0/tan( 0.5 * gamma/kT )*one - iota );

  for(k=1; k<=KK; k++){

    gamma_matsubara[k] = 2.0*k*M_PI*kT;

    double g = gamma_matsubara[k];
    c_matsubara[k] = 2*eta*kT * (gamma* g/(g*g - gamma*gamma));
  }


}



}// namespace libheom
}// namespace libdyn
}// liblibra

