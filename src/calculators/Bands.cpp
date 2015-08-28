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

#include "Bands.h"
#include "Fermi.h"

namespace libcalculators{


void convert_1(boost::python::list bands,  vector< pair<int,double> >& int_bands){
// Converst list to vector< pair<int, double> >
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

  // Now create list of lists
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

  int Norb = E->num_of_cols;

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

  libmmath::libspecialfunctions::merge_sort(in,bands);  

  if(0){ // For debug
    cout<<"Ordered bands:\n";
    for(i=0;i<Norb;i++){  cout<<"i= "<<i<<" orb_indx = "<<in[i].first<<" E[i]= "<<in[i].second<<endl; }
  }

  
}// order_bands(...)


boost::python::list order_bands(MATRIX E){

  vector< pair<int,double> > bands;
  order_bands(&E, bands);

  return convert_2(bands);

}


void populate_bands(double Nel, double degen, double kT, double etol, int pop_opt,
         vector< pair<int,double> >& bands,vector< pair<int,double> >& occ){

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
