#include "Bands.h"
#include "Fermi.h"

namespace libcalculators{


void order_bands(int Norb, MATRIX* E, vector< pair<int,double> >& bands){

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


void populate_bands(int Nocc, int Norb, int degen, double Nel, int pop_opt, double kT, double etol,
         vector< pair<int,double> >& bands,vector< pair<int,double> >& occ){

  int i;
  occ = bands;    // maybe not the best way, but it allocates memory

  if(pop_opt==0){  // integer populations - populate first Nocc bands

    for(i=0;i<Norb;i++){
      occ[i].first = bands[i].first;

      if(i<Nocc){ occ[i].second = (double)degen; }
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
