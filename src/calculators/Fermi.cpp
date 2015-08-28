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

#include "Fermi.h"

namespace libcalculators{

double fermi_population(double e,double ef,double degen, double kT){
// 
// Compute Fermi-based occupation of energy level with energy e
// ef - is the assumed Fermi energy
//  double kT = 0.00095; // kT in a.u. for T = 300 K

  double argg = ((e-ef)/kT);
  double pop;

  if(argg>50.0){ pop = 0.0;  }
  else if(argg<-50.0){   pop = degen;   }
  else{ pop = degen/(1.0 + std::exp(argg));  }  

  return pop;

}// fermi_population



double fermi_integral(std::vector<double>& bnds, double ef, double degen, double kT){
// Compute integral(sum):
//
//   sum [degen / (1 + exp(e-ef))]   
//    i
//
// where N - is a number of electrons (valence)
//  double kT = 0.00095; // kT in a.u. for T = 300 K

  int Norb = bnds.size();

  double sum = 0.0;
  for(int i=0;i<Norb;i++){

    double argg = ((bnds[i] - ef)/kT);

    if(argg>50.0){   }
    else if(argg<-50.0){   sum += degen;   }
    else{  sum += degen/(1.0 + std::exp(argg));   }  

  }// for i

  return sum;

}// double fermi_integral

double fermi_integral(boost::python::list bnds,double ef,double degen, double kT){

  int sz = boost::python::len(bnds);
  vector<double> int_bnds(sz,0.0);
  for(int i=0;i<sz;i++){ 
    int_bnds[i] = boost::python::extract<double>(bnds[i]);
  }

  return fermi_integral(int_bnds, ef, degen, kT);

}



double fermi_energy(std::vector<double>& bnds,double Nel,double degen, double kT, double etol){
// Computes Fermi energy by solving equation
// fermi_integral( ...  ef ... )  = Nel
//
// Using bisection method

  double ef_l,ef_m,ef_r,i_l,i_m,i_r;
  double err = 2.0*etol;
  int Norb = bnds.size();

  ef_l = bnds[0] - 10.0;
  ef_r = bnds[Norb-1] + 10.0;

  i_l = fermi_integral(bnds,ef_l,degen,kT) - Nel;
  i_r = fermi_integral(bnds,ef_r,degen,kT) - Nel;

  do{
    ef_m = 0.5*(ef_l + ef_r);
    i_m = fermi_integral(bnds,ef_m,degen,kT) - Nel;

    if(0){ // For debug
      cout<<"ef_l= "<<ef_l<<" ef_m= "<<ef_m<<" ef_r= "<<ef_r<<endl;
      cout<<"i_l= "<<i_l<<" i_m= "<<i_m<<" i_r= "<<i_r<<endl;
    }

    int var;
    if(i_m*i_r<=0.0 && i_m*i_l>=0.0){ var = 1; }
    else if(i_m*i_r>=0.0 && i_m*i_l<=0.0){ var = 2; }
    else{ cout<<"Error in fermi_energy\n"; exit(0); }

    switch(var){
      case 1: {i_l = i_m; ef_l = ef_m; break; }
      case 2: {i_r = i_m; ef_r = ef_m; break; }
      default: break;
    }

    err = 0.5*(ef_r-ef_l); 
//    err = 0.5*(i_r - i_l);
    
  }while(err>etol);
 
  return 0.5*(ef_l+ef_r);

}// double fermi_energy

double fermi_energy(boost::python::list bnds,double Nel,double degen, double kT, double etol){

  int sz = boost::python::len(bnds);
  vector<double> int_bnds(sz,0.0);
  for(int i=0;i<sz;i++){ 
    int_bnds[i] = boost::python::extract<double>(bnds[i]);
  }

  return fermi_energy(int_bnds, Nel, degen, kT, etol);

}


}// namespace libcalculators

