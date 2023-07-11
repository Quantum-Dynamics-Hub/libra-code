/*********************************************************************************
* Copyright (C) 2019-2020 Alexey V. Akimov
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
#include "../../util/libutil.h"
#include "../../math_linalg/liblinalg.h"
#include "../../math_specialfunctions/libspecialfunctions.h"
#include "libheom.h"

namespace liblibra{
namespace libdyn{
namespace libheom{

using namespace libutil;
using namespace liblinalg;
using namespace libspecialfunctions;



int compute_nn_tot(int d, int max_tier){
/**
  \brief Compute the total number of density matrices in the EOM hierarchy

  This number is given by:   (LL+nquant*KK)! / ( LL! * (nquant*KK)! )

  of course, we don't want to evaluate the factorial explicitly

  \param[in] L Level of the HEOM hierarchy
  \param[in] nquant (often denoted M) the number of system-bath coupling terms, in the case
  where each state is coupled to the bath intependently, this number is also equal to the number
  of quantum state, hence the current name
  \param[in] KK the number of the Matsubara terms

  d = nquant*(KK+1)  ! The expectation


*/

  int i;
  double tmp = 1.0;

  for(int i=1; i<=d; i++){
    tmp *= (max_tier+i)/((float)i);
  }

  return (int)tmp;
}



vector< vector<int> > gen_next_level(vector< vector<int> >& parents){
    /**
    This function is similar to gen_next_level, except it takes a list of
    parent vectors of integers and generates the list of the children for 
    all of the parent inputs.
    */

    int num_parents = parents.size();
    int num_children = parents[0].size(); // per parent

    vector< vector<int> > next_level;

    for(int i=0; i<num_parents; i++){
        for(int j=0;j<num_children; j++){

            vector<int> child = parents[i];
            child[j] += 1;

            if(! is_included(child, next_level) ){  next_level.push_back(child);  }

        }// for j
    }// for i

    return next_level;
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
            all_vectors.push_back( parents[i] );
            vector<int> coord(2,0);
            coord[0] = tier;
            coord[1] = iparent;
            all_coordinates.push_back( coord );
            iparent += 1;            

        }// for i

        tier_nums.push_back( all_vectors.size() ) ;

        vector< vector<int> > new_parents = gen_next_level(parents);
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



vector<int> filter(vector<CMATRIX>& rho, vector<int>& adm_list, double tolerance, int do_zeroing){
/**
   Gives the list of indices for the density matrices which are assumed to
   be non-zero, also returns the list of flags for whether a given ADM is non-zero

   Shi, Chen, Nan, Xu, Yan, JCP 130, 084105 (2009)
*/

  int nn_tot = rho.size(); 
  vector<int> nonzero(nn_tot, 0);  
  adm_list.clear();

  for(int n=0; n<nn_tot; n++){

    if( abs(rho[n].max_elt() ) > tolerance ){
      nonzero[n] = 1;
      adm_list.push_back(n);
    }
    else{   
      if(do_zeroing==1){
        rho[n] = 0.0;  
      }
    }

  }// for n

  return nonzero;
}





CMATRIX compute_deriv_n(int n, vector<CMATRIX>& rho, CMATRIX& Ham, vector<CMATRIX>& el_phon_coupl,
        double eta, double temperature, 
        vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara,
        int truncation_scheme, complex<double> truncation_prefactor, int do_scale, vector<int>& nonzero,
        vector< vector<int> >& nvectors, vector< vector<int> >& vec_plus, vector< vector<int> >& vec_minus        
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

  CMATRIX sum_rho_over_k(nquant, nquant);
  CMATRIX sum1_c_rho_over_k(nquant, nquant);
  CMATRIX sum2_c_rho_over_k(nquant, nquant);

  double scaling_factor = 1.0;
  complex<double> iota(0.0, 1.0);
  double kB = boltzmann/hartree;
 

  if(nonzero[n]==1){
    //=============== Liouvillian =====================
    drho_n_dt = -iota * ( Ham * rho[n] - rho[n] * Ham ); 


    //============= Friction ======================
    complex<double> pref; 

    for(k=0; k<=KK; k++){
      int sum_over_m = 0;
      for(m=0; m<nquant; m++){  sum_over_m += nvectors[n][m*(KK+1) + k];  }
      pref +=  sum_over_m * gamma_matsubara[k];
    }// for k
    drho_n_dt -=  pref * rho[n];


    //============= Truncation ======================     
    if(truncation_scheme==1 || truncation_scheme==2 || 
       truncation_scheme==3 || truncation_scheme==4){
      // All the following options take the same functional form, except for 
      // the prefactor is computed differently
      // 1 - according to Schulten, with real part of Matsubara terms
      // 2 - according to Schulten, but with full Matsubara terms
      // 3 - according to Shi, with real part of Matsubara terms
      // 4 - according to Shi, but with full Matsubara terms

      // Ihizaki-Tanimura scheme for truncation
      // JPSJ 74 3131, 2005
      //complex<double> pref(eta * kB * temperature/gamma_matsubara[0], 0.0);

      // Taking the real part in the line below is really important! otherwise, the result
      // is quite erratic
      //pref -= std::real( compute_matsubara_sum(gamma_matsubara, c_matsubara, KK) );
      
      for(m=0; m<nquant; m++){         
        commut = el_phon_coupl[m] * rho[n] - rho[n] * el_phon_coupl[m];
        commut = el_phon_coupl[m] * commut - commut * el_phon_coupl[m];
        drho_n_dt -= truncation_prefactor *commut;
      }
    }// if 

  }// not filtered out


  //===================== rho_n_plus terms ================

  for(m=0; m<nquant; m++){

    sum_rho_over_k = 0.0;
    for(k=0; k<=KK; k++){

      nplus = vec_plus[n][m*(KK+1) + k];
      if(nplus!=-1){ 
        if(nonzero[nplus]==1){

          scaling_factor = 1.0;
          if(do_scale==1){ // Scaled HEOM approach from JCP 130, 084105 (2009)
            scaling_factor = sqrt((nvectors[n][m*(KK+1) + k] + 1.0 )*fabs(c_matsubara[k]));
          } 

          sum_rho_over_k  += scaling_factor * rho[nplus];  

        }// not filtered out
      }// a valid entry 
    }// for k
    drho_n_dt -=  iota*( el_phon_coupl[m] * sum_rho_over_k - sum_rho_over_k * el_phon_coupl[m]);

  }// for m


  //===================== rho_n_minus terms ================

  for(m=0; m<nquant; m++){
    sum1_c_rho_over_k = 0.0;  
    sum2_c_rho_over_k = 0.0;

    for(k=0; k<=KK; k++){
      nminus = vec_minus[n][m*(KK+1) + k];

      if(nminus!=-1){
        if(nonzero[nminus]==1){

          scaling_factor = double(nvectors[n][m*(KK+1) + k]);
          if(do_scale==1){ // Scaled HEOM approach from JCP 130, 084105 (2009)

            double denom = fabs(c_matsubara[k]);
            if(denom > 1e-50){
              scaling_factor = sqrt( nvectors[n][m*(KK+1) + k] / denom );
            }
            else{  scaling_factor = 0.0; }

          } 

          sum1_c_rho_over_k  += scaling_factor * c_matsubara[k] * rho[nminus];
          sum2_c_rho_over_k  += scaling_factor * std::conj(c_matsubara[k]) * rho[nminus];
        }// not filtered out
      }// a valid entry

    }// for k
    drho_n_dt -=  iota*( el_phon_coupl[m] * sum1_c_rho_over_k - sum2_c_rho_over_k * el_phon_coupl[m]);

  }// for m


  return drho_n_dt;

}



CMATRIX compute_heom_derivatives(CMATRIX& RHO, bp::dict prms){
    /**

    RHO - CMATRIX((nn_tot)*nquant, nquant) - all the density matrices stacked together in a super-column
    prms - Python dict - control parameters
    
    */

    int n_threads = 1;
    int truncation_scheme = 1;
    int do_scale = 1;


    int i, n;
    int nquant = RHO.n_cols;
    int nn_tot = int(RHO.n_rows / nquant);


    CMATRIX Ham(nquant, nquant);
    double eta = 0.0;
    double gamma = 0.0;
    double temperature = 300.0;
    int KK; KK = 0;

    vector<double> gamma_matsubara;
    vector< complex<double> > c_matsubara;
    intList2 nvec, nvec_plus, nvec_minus;
    intList zero, nonzero, adm_list;

    adm_list = vector<int>(nn_tot, 0);
    for(i=0;i<nn_tot; i++){ adm_list[i] = i; }

    vector<CMATRIX> el_phon_couplings; 
    int el_phon_couplings_status = 0; // whether we have read them

    
    std::string key;
    for(i=0;i<len(prms.values());i++){

      key = extract<std::string>(prms.keys()[i]);

      if(key=="KK"){  KK = extract< int >(prms.values()[i]); }
      if(key=="Ham"){  Ham = extract<CMATRIX>(prms.values()[i]); }
      if(key=="eta"){  eta = extract<double>(prms.values()[i]); }
      if(key=="gamma"){  gamma = extract<double>(prms.values()[i]); }
      if(key=="temperature"){  temperature = extract<double>(prms.values()[i]); }
      if(key=="gamma_matsubara"){  gamma_matsubara = extract< doubleList >(prms.values()[i]); }
      if(key=="c_matsubara"){  c_matsubara = extract< complexList >(prms.values()[i]); }
      if(key=="truncation_scheme"){  truncation_scheme = extract< int >(prms.values()[i]); }
      if(key=="do_scale"){  do_scale = extract< int >(prms.values()[i]); }
      if(key=="el_phon_couplings"){  
        el_phon_couplings = extract< CMATRIXList >(prms.values()[i]); 
        el_phon_couplings_status = 1;
      }
      if(key=="nvec"){  nvec = extract< intList2 >(prms.values()[i]); }
      if(key=="nvec_plus"){  nvec_plus = extract< intList2 >(prms.values()[i]); }
      if(key=="nvec_minus"){ nvec_minus = extract< intList2 >(prms.values()[i]); }
      if(key=="nonzero"){  nonzero = extract< intList >(prms.values()[i]); }
      if(key=="adm_list"){  adm_list = extract< intList >(prms.values()[i]); }

      if(key=="num_threads"){ n_threads = extract<int>(prms.values()[i]); }

    }

    if(el_phon_couplings_status==0){
      el_phon_couplings = initialize_el_phonon_couplings(nquant);
    }

    complex<double> truncation_prefactor(0.0, 0.0);
    complex<double> iota(0.0, 1.0);
    double kB = boltzmann/hartree;

    if(truncation_scheme==1 || truncation_scheme==2){
      // Ihizaki-Tanimura scheme for truncation
      // JPSJ 74 3131, 2005

      truncation_prefactor = complex<double>(eta * kB * temperature/gamma_matsubara[0], 0.0);

      if(truncation_scheme==1){
        // 1 - according to Schulten, with real part of Matsubara terms
        // Taking the real part in the line below is really important! otherwise, the result
        // is quite erratic
        truncation_prefactor -= std::real( compute_matsubara_sum(gamma_matsubara, c_matsubara, KK) );
      }

      if(truncation_scheme==2){
        // 2 - according to Schulten, but with full Matsubara terms
        // But this is how it is supposed to be
        truncation_prefactor -= compute_matsubara_sum(gamma_matsubara, c_matsubara, KK);
      }
    }

    if(truncation_scheme==3 || truncation_scheme==4){
      // 3 - according to Shi, with real part of Matsubara terms
      // 4 - according to Shi, but with full Matsubara terms

      // This is a bit unoptimized approach to call the bath setups here
      int more = 200;
      int KK_ext = KK + more;
      vector<double> gamma_matsubara_ext;
      vector< complex<double> > c_matsubara_ext;
      setup_bath(KK_ext, eta, gamma, temperature, gamma_matsubara_ext, c_matsubara_ext);

      complex<double> sum_KK; sum_KK = compute_matsubara_sum(gamma_matsubara, c_matsubara, KK);
      complex<double> sum_KK_ext; sum_KK_ext = compute_matsubara_sum(gamma_matsubara_ext, c_matsubara_ext, KK+more);

      if(truncation_scheme==3){
        truncation_prefactor = sum_KK_ext - sum_KK;
      }

      if(truncation_scheme==4){
        truncation_prefactor = std::real( sum_KK_ext - sum_KK);
      }


    }


    CMATRIX dRHO(nn_tot*nquant, nquant);
    CMATRIX drho(nquant, nquant);
    vector<CMATRIX> rho_unpacked(nn_tot, CMATRIX(nquant, nquant));
    vector<CMATRIX> drho_unpacked(nn_tot, CMATRIX(nquant, nquant));

    unpack_mtx(rho_unpacked, RHO);

    int nn_nonzero = adm_list.size();

    #pragma omp parallel num_threads(n_threads)
    {
        #pragma omp for
        for(int n1=0; n1<nn_nonzero; n1++){
            drho_unpacked[adm_list[n1]] = compute_deriv_n(adm_list[n1], rho_unpacked, Ham, el_phon_couplings, eta, temperature, 
                gamma_matsubara, c_matsubara, truncation_scheme, truncation_prefactor, do_scale, nonzero, nvec, nvec_plus, nvec_minus);
        }
    }


    pack_mtx(drho_unpacked, dRHO);

    return dRHO;

}




vector<CMATRIX> initialize_el_phonon_couplings(int nquant){
/**
  Compute the matrices  Q_m = |m><m| - the simplest type of the coupling of electrons to phonons

*/

  vector<CMATRIX> projectors(nquant, CMATRIX(nquant, nquant));

  for(int m=0; m<nquant; m++){

    projectors[m] = 0.0;
    projectors[m].set(m,m, complex<double>(1.0, 0.0));  // matrix |m><m|
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


void setup_bath(int KK, double eta, double gamma, double temperature,
                vector<double>& gamma_matsubara, vector< complex<double> >& c_matsubara){
/** 
  KK - (KK+1) is the number of Matsubara frequencies
  eta - reorganization energy of the bath [Ha]
  gamma - the system-bath interaction frequency [a.u.]
  temperature - the bath temperature [K]
*/

  int k, i;

  complex<double> one(1.0, 0.0);
  complex<double> iota(0.0, 1.0);

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







void unpack_mtx(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO){
/**
    RHO = CMATRIX(num_adm * nquant , nquant)
*/

    int i, n;
    int nquant = RHO.n_cols;
    int num_adm = int(RHO.n_rows / nquant);

    if(rho_unpacked.size() != num_adm ) {   cout<<"ERROR in unpack_rho: \n"; exit(0); }
    
    vector<int> sub_x(nquant, 0);
    vector<int> sub_y(nquant, 0);
    for(i = 0; i < nquant; i++){    sub_y[i] = i;  }

    
    for(n = 0; n < num_adm; n++){

      for(i = 0; i<nquant; i++){    sub_x[i] = n*nquant + i;  }
      pop_submatrix(RHO, rho_unpacked[n], sub_x, sub_y );
        
    }

}


void pack_mtx(vector<CMATRIX>& rho_unpacked, CMATRIX& RHO){

    int i, n;
    int nquant = RHO.n_cols;
    int num_adm = int(RHO.n_rows / nquant);

    if(rho_unpacked.size() != num_adm) {   cout<<"ERROR in pack_rho: \n"; exit(0); }
    

    vector<int> sub_x(nquant, 0);
    vector<int> sub_y(nquant, 0);
    for(i = 0; i<nquant; i++){    sub_y[i] = i;  }

    
    for(n = 0; n < num_adm; n++){

      for(i = 0; i < nquant; i++){    sub_x[i] = n*nquant + i;  }
      push_submatrix(RHO, rho_unpacked[n], sub_x, sub_y );
        
    }
}


void scale_rho(vector<CMATRIX>& rho, vector<CMATRIX>& rho_scaled, bp::dict prms){


    vector< complex<double> > c_matsubara;
    intList2 nvec;

    int i, n,m,k;
    std::string key;

    for(i=0;i<len(prms.values());i++){

      key = extract<std::string>(prms.keys()[i]);

      if(key=="nvec"){  nvec = extract< intList2 >(prms.values()[i]); }
      if(key=="c_matsubara"){  c_matsubara = extract< complexList >(prms.values()[i]); }

    }

    int nn_tot = rho.size();
    int nquant = rho[0].n_cols;
    int KK = c_matsubara.size() - 1; 


    for(n=0; n<nn_tot; n++){

      double scaling_factor = 1.0;

      for(m=0; m<nquant; m++){
        for(k=0; k<=KK; k++){

          int nn = nvec[n][m*(KK+1) + k];
          scaling_factor *= FACTORIAL( nn );
          scaling_factor *= FAST_POW( fabs(c_matsubara[k]), nn);          
        }// for k
      }// for m

      rho_scaled[n] = rho[n] * sqrt(1.0/scaling_factor);
   
    }// for n

}


void inv_scale_rho(vector<CMATRIX>& rho, vector<CMATRIX>& rho_scaled, bp::dict prms){


    vector< complex<double> > c_matsubara;
    intList2 nvec;

    int i, n,m,k;
    std::string key;

    for(i=0;i<len(prms.values());i++){

      key = extract<std::string>(prms.keys()[i]);

      if(key=="nvec"){  nvec = extract< intList2 >(prms.values()[i]); }
      if(key=="c_matsubara"){  c_matsubara = extract< complexList >(prms.values()[i]); }

    }

    int nn_tot = rho.size();
    int nquant = rho[0].n_cols;
    int KK = c_matsubara.size() - 1; 


    for(n=0; n<nn_tot; n++){

      double scaling_factor = 1.0;

      for(m=0; m<nquant; m++){
        for(k=0; k<=KK; k++){

          int nn = nvec[n][m*(KK+1) + k];
          scaling_factor *= FACTORIAL( nn );
          scaling_factor *= FAST_POW( fabs(c_matsubara[k]), nn);          
        }// for k
      }// for m

      rho[n] = rho_scaled[n] * sqrt(scaling_factor);
   
    }// for n

}





}// namespace libheom
}// namespace libdyn
}// liblibra

