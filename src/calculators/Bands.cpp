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
  \file Bands.cpp
  \brief The file implements functions for ordering, converting, and printing bands (energies and occupations) information
    
*/

#include "Bands.h"
#include "Fermi.h"
#include "../math_specialfunctions/libspecialfunctions.h"
#include "../Units.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;


/// libcalculators namespace
namespace libcalculators{


void convert_1(boost::python::list bands,  vector< pair<int,double> >& int_bands){
/**
  \brief Auxiliary converter function

  This function converts Python list into vector of pairs, so we can conveniently 
  conver input in Python-friendly functions, to use internal functions

  \param[in] bands The Python input - it is expected to be the list of 2-element lists
  \param[out] int_bands The C++ output - the vector of pairs
*/

  int Norb = len(bands);

  for(int i=0;i<Norb;i++){
    // Extract data from external list of lists
    boost::python::list bnd = boost::python::extract<boost::python::list>(bands[i]);

    int indx = boost::python::extract<int>(bnd[0]);
    double en = boost::python::extract<double>(bnd[1]);

    // Create internal representation
    std::pair<int,double> x;
    x = make_pair(indx, en);

    int_bands.push_back(x);

  }// for i

}

boost::python::list convert_2( vector< pair<int,double> >& bands){
/**
  \brief Auxiliary converter function

  This function converts C++ vector of pairs into Python list

  \param[in] bands The C++ input - the vector of pairs
  The function returns a list of 2-element lists
*/

  boost::python::list res;

  for(int i=0;i<bands.size();i++){
    boost::python::list bnd; 
 
    bnd.append(bands[i].first);
    bnd.append(bands[i].second);

    res.append(bnd);
  }

  return res;

}



void order_bands(MATRIX* E, vector< pair<int,double> >& bands){
/**
  \brief Ordering of bands 

  This function takes the diagonal elements of the E matrix, orders them and packs into a vector of pairs. Each 
  pair is contains the index of the state (in the original, petentially disordered, structure E) and the corresponding
  value of the matrix element.

  \param[in] E The pointer to the diagonal matrix with energies
  \param[in,out] bands The packed structure containing ordered eigenvalues and their original ordering indices
*/


  int Norb = E->n_cols;

  std::vector< pair<int,double> > in;
  std::pair<int,double> x;
  int i;

  for(i=0;i<Norb;i++){
    x = make_pair(i, E->M[i*Norb+i]);
    in.push_back(x);
  }

  if(0){ // For debug
    cout<<"Unordered bands:\n";
    for(i=0;i<Norb;i++){  cout<<"i= "<<i<<" orb_indx = "<<in[i].first<<" E[i]= "<<in[i].second<<endl; }
  }

  libspecialfunctions::merge_sort(in,bands);  

  if(0){ // For debug
    cout<<"Ordered bands:\n";
    for(i=0;i<Norb;i++){  cout<<"i= "<<i<<" orb_indx = "<<in[i].first<<" E[i]= "<<in[i].second<<endl; }
  }

  
}// order_bands(...)


boost::python::list order_bands(MATRIX E){
/**
  \brief Ordering of bands (Python-friendly version)

  This function takes the diagonal elements of the E matrix, orders them and packs into a vector of pairs. Each 
  pair is contains the index of the state (in the original, petentially disordered, structure E) and the corresponding
  value of the matrix element.

  \param[in] E The diagonal matrix with energies
  The function returns the packed structure (represented as a list of 2-element lists) containing ordered eigenvalues 
  and their original ordering indices
*/


  vector< pair<int,double> > bands;
  order_bands(&E, bands);

  return convert_2(bands);

}


void populate_bands(double Nel, double degen, double kT, double etol, int pop_opt,
         vector< pair<int,double> >& bands,vector< pair<int,double> >& occ){
/**
  \brief Compute populations of bands

  The function computes populations of bands, depending on the chosen poulation scheme and the corresponding parameters.

  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[in] bands The packed structure containing ordered eigenvalues and their original ordering indices
  \param[in,out] occ The packed structure containing populations of the energy levels packed in bands

*/

  int Norb = bands.size();
  int i;
  occ = bands;    // maybe not the best way, but it allocates memory

  if(pop_opt==0){  // integer populations - populate first Nocc bands

    double sum = 0.0;
    for(i=0;i<Norb;i++){
      occ[i].first = bands[i].first;

      if(sum<Nel){ occ[i].second = degen; sum+= degen; }
      else{ occ[i].second = 0.0; }

    }// for i
  }// pop_opt = 0

  else if(pop_opt==1){  // Fermi distribution (fractional populations possible)

    vector<double> e(Norb,0.0);  for(i=0;i<Norb;i++){ e[i] = bands[i].second; }
 
    double E_f = fermi_energy(e, Nel, degen, kT, etol);

    for(i=0;i<Norb;i++){
      occ[i].first = bands[i].first;
      occ[i].second = fermi_population(e[i],E_f,degen,kT);

    }// for i
  }// pop_opt = 1


  if(0){  // for debug, now inactive
    cout<<"In populate:\n";
    cout<<"Occupation numbers:\n";
    for(i=0;i<Norb;i++){
      cout<<"i= "<<i<<" orbital index= "<<occ[i].first<<" occupation_number= "<<occ[i].second<<endl;
    }
  }// if 


}// populate_bands(...)


boost::python::list populate_bands(double Nel, double degen, double kT, double etol, int pop_opt,
         boost::python::list bands){
/**
  \brief Compute populations of bands (Python-friendly version)

  The function computes populations of bands, depending on the chosen poulation scheme and the corresponding parameters.

  \param[in] Nel The number of electrons
  \param[in] degen Dengeneracy of orbitals (the maximal number of electrons that can occupay one orbital)
  \param[in] kT  Broadening factor for Fermi distribution
  \param[in] etol Tolerance level (stop when 0.5*|e_f(old) - e_f(new)|<tol)
  \param[in] pop_opt The flag controlling the population scheme
             pop_opt = 0 - integer occupation numbers will be used (good in many standard cases)
             pop_opt = 1 - fractional occupations will be possible (can help in difficult cases)
  \param[in] bands The packed structure containing ordered eigenvalues and their original ordering indices
  The function returns the packed structure (list of 2-element lists) containing populations of the energy levels packed in bands

*/


  // General bands:
  int Norb = len(bands);
  vector< pair<int,double> > int_bands;

  // Extract data from external list of lists and create internal representation
  convert_1(bands, int_bands);

  // Call c++ version
  vector< pair<int,double> > occ;
  populate_bands(Nel, degen, kT, etol, pop_opt, int_bands, occ);


  // Now create list of lists
  boost::python::list res;

  for(int i=0;i<Norb;i++){
    boost::python::list oc; 
 
    oc.append(occ[i].first);
    oc.append(occ[i].second);

    res.append(oc);
  }

  return res;
}



void show_bands(int Norb, int Nocc, vector< pair<int,double> >& bands,vector< pair<int,double> >& occ){
/**
  \brief Formatter printout of the bands

  The function will print energy levels and the corresponding populations

  \param[in] Norb The number of orbitals to show (starting with index 0) - usually all
  \param[in] Nocc The number of occupied orbitals - we just print this number. This parameter is not used in any other way.
  \param[in] bands The packed structure containing ordered eigenvalues and their original ordering indices
  \param[in] occ The packed structure containing populations of the energy levels packed in bands

*/

// Show only Norb bands

  int sz = Norb;
  cout<<"# of bands = "<<sz<<"\n";
  cout<<"# of occupied = "<<Nocc<<endl;

  for(int i=0;i<sz;i++){
    cout<<" band# = "<<i<<" mo indx ="<<bands[i].first
        <<" energy(eV) = "<<bands[i].second/eV
        <<" energy(a.u.) = "<<bands[i].second
        <<" occupation = "<<occ[i].second<<endl;
  }

}// show_bands(...)

}// namespace libcalculators

}// liblibra


