/*********************************************************************************
* Copyright (C) 2017-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Interactions.h
  \brief The file describes functions and classes for molecular-mechanical Hamiltonian calculations 
*/

#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "../../cell/libcell.h"
#include "../../pot/libpot.h"
#include "../../chemobjects/libchemobjects.h"
#include "../../forcefield/libforcefield.h"


/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{

using namespace libchemobjects;
using namespace libchemobjects::libchemsys;
using namespace libpot;
using namespace libcell;
using namespace libforcefield;

using namespace boost::python;



map<std::string, double> dict2map(boost::python::dict d);

template<class X>
vector<X> value_subset(vector<X>& r_in, vector<int>& indxs){
/**
  This is an auxiliary function that sets the values in r_out to the values of  
  a sub-set of pointers in r_in. The subset is determined by the vector of indices "indxs"
  
*/
    vector<X> r_out;

    int sz_in  = r_in.size();
    int sz_sub = indxs.size();


    for(int i=0;i<sz_sub;i++){

        if(indxs[i]>=sz_in){ 
            cout<<"Error: Sub-set indxed "<<indxs[i]<<" is beyond the size of the set"<<sz_in<<endl;    
            exit(0);
        }// if

      r_out.push_back( r_in[indxs[i]] );
     
    }// for i

    return r_out;

}// value_subset


template<class X>
void address_subset(vector<X*>& r_out, vector<X*>& r_in, vector<int>& indxs){
/**
  This is an auxiliary function that sets the pointers in r_out to the 
  a sub-set of pointers in r_in. The subset is determined by the vector of indices "indxs"
  
*/
    int sz_out = r_out.size();
    int sz_in  = r_in.size();

    if(sz_out!=indxs.size()){
        cout<<"The size of the output vector "<<sz_out<<" is not equal to the size to the subset "<<indxs.size()<<endl;
        cout<<"Are you sure you have allocated memory for the r_out?\n";
        exit(0);
    }

    for(int i=0;i<sz_out;i++){

        if(indxs[i]>=sz_in){ 
            cout<<"Error: Sub-set indxed "<<indxs[i]<<" is beyond the size of the set"<<sz_in<<endl;    
            exit(0);
         }// if

        r_out[i] = r_in[indxs[i]];
     
    }// for i

}// create_subset

template<class X>
void address_subset(vector<X*>& r_out, vector<X>& r_in, vector<int>& indxs){
/**
  This is an auxiliary function that sets the pointers in r_out to the addresses of a 
  a sub-set of VECTORS in r_in. The subset is determined by the vector of indices "indxs"
  
*/
    int sz_out = r_out.size();
    int sz_in  = r_in.size();

    if(sz_out!=indxs.size()){
        cout<<"The size of the output vector "<<sz_out<<" is not equal to the size to the subset "<<indxs.size()<<endl;
        cout<<"Are you sure you have allocated memory for the r_out?\n";
        exit(0);
    }

    for(int i=0;i<sz_out;i++){

        if(indxs[i]>=sz_in){ 
            cout<<"Error: Sub-set indxed "<<indxs[i]<<" is beyond the size of the set"<<sz_in<<endl;    
            exit(0);
         }// if

        r_out[i] = &r_in[indxs[i]];
     
    }// for i

}// create_subset




//vector<VECTOR> value_subset(vector<VECTOR>& r_in, vector<int>& indxs);
//void address_subset(vector<VECTOR*>& r_out, vector<VECTOR*>& r_in, vector<int>& indxs);
//void address_subset(vector<VECTOR*>& r_out, vector<VECTOR>& r_in, vector<int>& indxs);


class Interaction_N_Body{

public:

    // Data
    int Nbody;           ///< order of the N-body interactions

    vector<VECTOR*> r;   ///< list of the pointers to the atomic coordinates
    vector<VECTOR*> t;   ///< list of the pointers to the vectors representing periodic translations of the actoms
    vector<VECTOR*> f;   ///< list of the pointers to the forces acting on atoms
    vector<double*> q;   ///< list of the pointers to the atomic charges
    MATRIX* Hess;        ///< Hessian
    vector<int> Hess_stenc; ///<stencil for the Hessian due to given interaction

    int is_active;        ///< Flag showing if this interaction is active
    int int_type;         ///< Type of interaction
    int functional;       ///< Type of potential to use for given interaction type

/**
  In the derived classes:

                               int_type           functional

   Interaction_N_Body             1               Ewald_3D      (0)
                                                   vdw_LJ       (1)
                                                   vdw_LJ1      (2)
                                                  LJ_Coulomb    (3)
                                  
   Interaction_2_Body             2

      Bond_Interaction            20               Harmonic     (0)
                                                    Quartic     (1)
                                                     Morse      (2)

      VdW_Interaction             21                   LJ       (0)
                                                  Buffered14_7  (1)
                                                     Morse      (2)

      Elec_Interaction            22                 Coulomb    (0)


   Interaction_3_Body             3

      Angle_Interaction           30                Harmonic    (0)
                                                    Fourier     (1)
                                                Fourier_General (2)
                                                Fourier_Special (3)
                                                  Harmonic_Cos  (4)
                                            Harmonic_Cos_General(5)
                                                      Cubic     (6)

   Interaction_4_Body             4

      Dihedral_Interaction        40               General0     (0)
                                                   General1     (1)
                                                   General2     (2)
                                                   General3     (3)
                                                   Fourier0     (4)
                                                   Fourier1     (5)

      OOP_Interaction             41               Fourier      (0)
                                                   Wilson       (1)
                                                   Harmonic     (2)

*/

    double energy;        ///< Energy for the given Hamiltonian and the status flag
    MATRIX3x3 stress_at;  ///< atomic stress tensor and status



    // Constructor & Destructor
    Interaction_N_Body();
    Interaction_N_Body(int Nbody_);

    // Methods
    void set_coords(VECTOR& r_, int indx_);
    void set_coords(vector<VECTOR*>& r_, vector<int>& indxs_);
    void set_coords(vector<VECTOR>& r_, vector<int>& indxs_);

    void set_transl(VECTOR& t_, int indx_);
    void set_transl(vector<VECTOR*>& t_, vector<int>& indxs_);
    void set_transl(vector<VECTOR>& t_, vector<int>& indxs_);

    void set_forces(VECTOR& f_, int indx_);
    void set_forces(vector<VECTOR*>& f_, vector<int>& indxs_);
    void set_forces(vector<VECTOR>& f_, vector<int>& indxs_);

    void set_charges(double& q_, int indx_);
    void set_charges(vector<double*>& q_, vector<int>& indxs_);
    void set_charges(vector<double>& q_, vector<int>& indxs_);

    void set_hessian(MATRIX& Hess_, vector<int>& Hess_stenc_);


    void set_functional(std::string f);
  
};









}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra


#endif // INTERACTIONS_H
