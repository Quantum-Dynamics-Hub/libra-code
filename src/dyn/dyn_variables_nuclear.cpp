/*********************************************************************************
* Copyright (C) 2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_variables_nuclear.cpp
  \brief The file implements the methods for nuclear dyn vars
*/

#include "dyn_variables.h"
#include "../util/libutil.h"
#include "../converters/libconverters.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libutil;
using namespace libconverters;


/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;



void dyn_variables::init_nuclear_dyn_var(bp::dict _params, Random& rnd){
/**
    """
    Args:
        q ( list of doubles ): the mean values of coordinates for all DOFs [ units: a.u.]
        p ( list of doubles ): the mean values of momenta for all DOFs [ units: a.u. ]
        mass ( list of doubles ): masses of all nuclear DOFs [ units: a.u. ]

        params ( dictionary ): control parameters
 
            * **params["init_type"]** ( int ): the type of sampling of nuclear DOFs
     
                - 0 : initialize ```ntraj``` identical copies of coordinates and momenta
 
                - 1 : keep all coordinates identical, but sample momenta from the normal 
                    distribution with a given width in each dimension

                - 2 : keep all momenta identical, but sample coordinates from the normal 
                    distribution with a given width in each dimension

                - 3 : sample both coordinates and momenta from the normal 
                    distributions with given widths in each dimension

            * **params["force_constant"]** ( list of double ): force constants involved in the Harmonic
                oscillator model: U = (1/2) * k * x^2, and omega = sqrt( k / m )
                These parameters define the Harmonic oscillator ground state wavefunctions from which
                the coordinates and momenta are sampled. [ units: a.u. = Ha/Bohr^2, default: [0.001] ]
                The length should be consistent with the length of Q, P and M arguments
                The ground state of the corresponding HO is given by:
                 psi_0 = (alp/pi)^(1/4)  exp(-1/2 * alp * x^2), with alp = m * omega / hbar
                The corresponding probability density is distributed as:
                Tully uses psi(x) ~ exp(-(x/sigma)^2) so 1/sigma^2 = alp/2 => alp = 2/sigma^2 = m sqrt (k/m) = sqrt( k * m )

                To make sigma = 20 / p0 => k * m = 4 /sigma^4 = 4/ (20/p0)^4 = 4* p^4 / 20^4 =  0.000025 * p^4 
                So: k =  0.000025 * p^4  / m
                

            * **params["ntraj"]** ( int ): the number of trajectories - the parameter defines the
                number of columns the output matrices will have [ default: 1 ]



        rnd ( Random ): random numbers generator object


    Returns:
        q, p, iM:  where:

            * q ( MATRIX(ndof, ntraj) ) : coordinates for all trajectories
            * p ( MATRIX(ndof, ntraj) ) : momenta for all trajectories
            * iM ( MATRIX(ndof, 1) ) : inverse masses of all DOFs (same across the trajectories)

    """
*/


  //# Read the parameters
  bp::list critical_params; 
  bp::dict default_params;
  bp::dict params(_params);

  default_params["init_type"] = 0;
  default_params["q"] = vector<double>(ndof, 0.0);
  default_params["p"] = vector<double>(ndof, 0.0);
  default_params["mass"] = vector<double>(ndof, 1.0);
  default_params["force_constant"] = vector<double>(ndof, 0.001);

  check_input(params, default_params, critical_params);


  int init_type;
  vector<double> _Q;
  vector<double> _P;
  vector<double> _M;
  vector<double> force_constant;

  int idof;

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);

    ///================= Computing Hamiltonian-related properties ====================
    if(key=="init_type") {  init_type = bp::extract<int>(params.values()[i]); }
    else if(key=="q") {  _Q = liblibra::libconverters::Py2Cpp<double>( bp::extract< bp::list >(params.values()[i]) ); }
    else if(key=="p") {  _P = liblibra::libconverters::Py2Cpp<double>( bp::extract< bp::list >(params.values()[i]) ); }
    else if(key=="mass") {  _M = liblibra::libconverters::Py2Cpp<double>( bp::extract< bp::list >(params.values()[i]) ); }
    else if(key=="force_constant") {  force_constant = liblibra::libconverters::Py2Cpp<double>( bp::extract< bp::list >(params.values()[i]) ); }
  }




  if( !(init_type==0 || init_type==1 || init_type==2 || init_type==3) ){
    cout<<"WARNINIG in init_nuclear_dyn_var: \
           the init_type = "<<init_type<<" is not known\
           Allowed values are: [0, 1, 2, 3]\n";
  }

  if(_Q.size() != _P.size()){
    cout<<"ERROR in init_nuclear_dyn_var: \
          the length of input Q is = "<<_Q.size()<<", \
          the length of input P is = "<<_P.size()<<", but they should be equal to each other\n";
    exit(0);
  }

  if(_Q.size() != _M.size()){
    cout<<"ERROR in init_nuclear_dyn_var: \
          the length of input Q is = "<<_Q.size()<<", \
          the length of input M is = "<<_M.size()<<", but they should be equal to each other\n";
    exit(0);
  }

  if(_Q.size() != force_constant.size()){
    cout<<"ERROR in init_nuclear_dyn_var: \
          the length of input Q is = "<<_Q.size()<<", \
          the length of input force_constant is = "<<force_constant.size()<<", but they should be equal to each other\n";
    exit(0);
  }


  /// At this point, it is safe to define ndof:
  for(idof=0; idof<ndof; idof++){
      iM->set(idof, 0, 1.0/_M[idof]);
  }


  // Mean values
  MATRIX mean_q(ndof, 1);
  MATRIX mean_p(ndof, 1);

  for(idof=0; idof<ndof; idof++){
    mean_q.set(idof, 0, _Q[idof]);
    mean_p.set(idof, 0, _P[idof]);
  }

  // Deviations
  MATRIX sigma_q(ndof, 1);
  MATRIX sigma_p(ndof, 1);

  if(init_type==0){
    for(idof=0; idof<ndof; idof++){
      sigma_q.set(idof, 0, 0.0);
      sigma_p.set(idof, 0, 0.0);
    }
  }
  else if(init_type==1){
    for(idof=0; idof<ndof; idof++){
      sigma_q.set(idof, 0, 0.0);
      sigma_p.set(idof, 0, sqrt( 0.5*sqrt((force_constant[idof]*_M[idof])) ));
    }
  }
  else if(init_type==2){
    for(idof=0; idof<ndof; idof++){
      sigma_q.set(idof, 0, sqrt( 0.5*sqrt(1.0/(force_constant[idof]*_M[idof])) ));
      sigma_p.set(idof, 0, 0.0);
    }
  }
  else if(init_type==3){
    for(idof=0; idof<ndof; idof++){
      sigma_q.set(idof, 0, sqrt( 0.5*sqrt(1.0/(force_constant[idof]*_M[idof])) ));
      sigma_p.set(idof, 0, sqrt( 0.5*sqrt((force_constant[idof]*_M[idof])) ));
    }
  }


  // Now sample the values
  liblibra::libspecialfunctions::sample(q, mean_q, sigma_q, rnd);
  liblibra::libspecialfunctions::sample(p, mean_p, sigma_p, rnd);

}


double dyn_variables::compute_average_kinetic_energy(){
  double res = 0.0;

  for(int itraj = 0; itraj < ntraj; itraj++){
    for(int idof = 0; idof < ndof; idof++){
      res += p->get(idof, itraj) * p->get(idof, itraj) * iM->get(idof, 0);
    }    
  }
  return 0.5*res/ float(ntraj);
}

double dyn_variables::compute_average_kinetic_energy(vector<int>& which_dofs){
  double res = 0.0;

  for(int itraj = 0; itraj < ntraj; itraj++){
    for(auto idof: which_dofs){
      res += p->get(idof, itraj) * p->get(idof, itraj) * iM->get(idof, 0);
    }    
  }
  return 0.5*res / float(which_dofs.size() );
}


double dyn_variables::compute_kinetic_energy(int itraj){
  double res = 0.0;

  for(int idof = 0; idof < ndof; idof++){
    res += p->get(idof, itraj) * p->get(idof, itraj) * iM->get(idof, 0);
  }
  res *= 0.5;

  return res;
}


double dyn_variables::compute_kinetic_energy(int itraj, vector<int>& which_dofs){
  double res = 0.0;

  for(auto idof : which_dofs){
    res += p->get(idof, itraj) * p->get(idof, itraj) * iM->get(idof, 0);
  }
  res *= 0.5;

  return res;
}



vector<double> dyn_variables::compute_kinetic_energies(){
  vector<double> res(ntraj, 0.0);

  for(int itraj = 0; itraj < ntraj; itraj++){
    for(int idof = 0; idof < ndof; idof++){
      res[itraj] += p->get(idof, itraj) * p->get(idof, itraj) * iM->get(idof, 0);
    }    
    res[itraj] *= 0.5;
  }

  return res;

}

vector<double> dyn_variables::compute_kinetic_energies(vector<int>& which_dofs){
  vector<double> res(ntraj, 0.0);

  for(int itraj = 0; itraj < ntraj; itraj++){
    for(auto idof : which_dofs){
      res[itraj] += p->get(idof, itraj) * p->get(idof, itraj) * iM->get(idof, 0);
    }    
    res[itraj] *= 0.5;
  }

  return res;

}





}// namespace libdyn
}// liblibra

