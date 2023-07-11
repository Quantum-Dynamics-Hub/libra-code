/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Hamiltonian_MM.h
  \brief The file describes functions and classes for molecular-mechanical Hamiltonian calculations 
*/

#ifndef HAMILTONIAN_MM_H
#define HAMILTONIAN_MM_H

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





class Hamiltonian_MM{
/**
  This class represents the classical interaction: bond (2-body), angle (3-body), etc... up to many-body interactions like
  Ewald summation in periodic box
*/


  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Hamiltonian_MM&); // Copies the content which is defined

//----------------- Types of interaction supported ---------------
  struct bond_interaction{
    int id1,id2;
    VECTOR* r1;  VECTOR* r2;
    VECTOR* g1;  VECTOR* g2;
    VECTOR* m1;  VECTOR* m2;
    VECTOR* f1;  VECTOR* f2;
    VECTOR* t1;  VECTOR* t2;  // periodic translations of each atoms
    double K,D,r0,alpha;
  };
  struct angle_interaction{
    int id1,id2,id3;
    VECTOR* r1;  VECTOR* r2; VECTOR* r3;
    VECTOR* g1;  VECTOR* g2; VECTOR* g3;
    VECTOR* m1;  VECTOR* m2; VECTOR* m3;
    VECTOR* f1;  VECTOR* f2; VECTOR* f3;
    VECTOR* t1;  VECTOR* t2; VECTOR* t3; // periodic translations of each atoms
    double k_theta, theta_0, cos_theta_0, C0,C1,C2;
    int coordination;
  }; 
  struct dihedral_interaction{
    int id1,id2,id3,id4;
    VECTOR* r1;  VECTOR* r2; VECTOR* r3; VECTOR* r4;
    VECTOR* f1;  VECTOR* f2; VECTOR* f3; VECTOR* f4;
    VECTOR* g1;  VECTOR* g2; VECTOR* g3; VECTOR* g4;
    VECTOR* m1;  VECTOR* m2; VECTOR* m3; VECTOR* m4;
    double Vphi,phi0,Vphi1,Vphi2,Vphi3;
    int opt,n;
  };
  struct oop_interaction{
    int id1,id2,id3,id4;
    VECTOR* r1;  VECTOR* r2; VECTOR* r3; VECTOR* r4;
    VECTOR* g1;  VECTOR* g2; VECTOR* g3; VECTOR* g4;
    VECTOR* m1;  VECTOR* m2; VECTOR* m3; VECTOR* m4;
    VECTOR* f1;  VECTOR* f2; VECTOR* f3; VECTOR* f4;    
    double K,C0,C1,C2,xi_0;
    int opt;
  };
  struct vdw_interaction{
    int id1,id2;
    VECTOR* r1;  VECTOR* r2;
    VECTOR* g1;  VECTOR* g2;
    VECTOR* m1;  VECTOR* m2;
    VECTOR* f1;  VECTOR* f2;
    double sigma,epsilon,D,r0,alpha,scale;
    int is_cutoff;
    double R_on,R_off;
    double R_on2,R_off2;
    vector<triple> images;  int is_images;
    triple central_translation; int is_central_translation;
    int time; // time since last recalculation of this pair
    VECTOR r1_old; VECTOR r2_old;   
    double* displr1_2;
    double* displr2_2; 
    double* displT_2;
  };
  struct elec_interaction{
    int id1,id2;
    VECTOR* r1;  VECTOR* r2;
    VECTOR* g1;  VECTOR* g2;
    VECTOR* m1;  VECTOR* m2;
    VECTOR* f1;  VECTOR* f2;
    double q1,q2,J,xi1,xi2,eps,delta,scale;   
    int is_cutoff;
    double R_on,R_off;
    double R_on2,R_off2;
  };
  struct gay_berne_interaction{
    int id1,id2;
    VECTOR* r1;  VECTOR* r2; VECTOR* u1; VECTOR* u2;
    VECTOR* f1;  VECTOR* f2; VECTOR* t1; VECTOR* t2;
    double di,dj,li,lj,e0,rat,dw,mu,nu;
  };
  struct mb_interaction{ // many(multi)-body interaction
    int sz,nexcl;
    int* id; int* excl1; int* excl2;
    double* scale;
    VECTOR** r;
    VECTOR** g;
    VECTOR** m; 
    VECTOR** f;
    double** q;
    double** epsilon;
    double** sigma;
    double** displr_2;
    double*  displT_2;
    int is_cutoff;
    double R_on,R_off;
    double R_on2,R_off2;
    double elec_etha;
    vector< vector<triple> > images;  int is_images;
    vector<triple> central_translation; int is_central_translation;
    vector< vector<quartet> > at_neib;
    vector< vector<excl_scale> > excl_scales; 
    int time; // time since last recalculation of this pair

  };
//-------------------------------------------------------------------

  int is_active;                                  // Flag showing if this interaction is active
  int int_type;              int is_int_type;     // Type of interaction
  int functional;            int is_functional;   // Type of potential to use for given interaction type
  int respa_type;            int is_respa_type;   // Type of the potential in RESPA methods (fast = 0, medium = 1, slow = 2)

  MATRIX3x3* Box;  
  MATRIX3x3 Box_old;
  int kx,ky,kz;          // translation vectors of first center
//  int kx2,ky2,kz2;       // translation vectors of second center
  bond_interaction*      data_bond; 
  angle_interaction*     data_angle;
  dihedral_interaction*  data_dihedral;
  oop_interaction*       data_oop;
  vdw_interaction*       data_vdw;
  elec_interaction*      data_elec;
  gay_berne_interaction* data_gay_berne;
  mb_interaction*        data_mb;

public:

  double energy;          int is_energy;     ///< Energy for the given Hamiltonian and the status flag
  MATRIX3x3 hessian;      int is_hessian;    ///< Hessian for the given Hamiltonian and the status flag
  MATRIX3x3 stress_at;    int is_stress_at;  ///< atomic stress tensor and status
  MATRIX3x3 stress_fr;    int is_stress_fr;  ///< fragmental stress tensor and status
  MATRIX3x3 stress_ml;    int is_stress_ml;  ///< molecular stress tensor and status

  //----------- Basic class operations ---------------------------
  // Defined in Hamiltonian_MM.cpp
  Hamiltonian_MM();                   ///< constructor
  Hamiltonian_MM(const Hamiltonian_MM&); ///< copy-constructor
 ~Hamiltonian_MM();                   ///< destructor
  Hamiltonian_MM& operator=(const Hamiltonian_MM&); ///< assignment operator
  friend int operator == (const Hamiltonian_MM& i1, const Hamiltonian_MM& i2);

  void show_info();

  //---------------------------------------------------------------
  int get_type()  { return int_type; }   ///< Returns the type of interaction for this Hamiltonian (bonds, angles, etc.)
  int get_status(){ return is_active; }  ///< Returns the "active" status of the interaction (Hamiltonian) - if it is deactivated, 
                                         ///< the interaction is not computed

  //------- Interface : Defined in Hamiltonian_MM_methods1.cpp -------------
  //General manipulations
  void activate()  { is_active = 1; }    ///<  Makes this interaction active
  void deactivate(){ is_active = 0; }    ///< Makes this interaction inactive
  void set_pbc(MATRIX3x3*,int,int,int);
  int is_origin();
  void set_respa_type(int int_type_,int respa_type_){ 
    if(int_type_==int_type && respa_type>=0){ respa_type = respa_type_; is_respa_type = 1; }
  }
  int get_respa_type(){ return respa_type; }   ///< Returns the RESPA type
  // 2, 3, 4 - atomic interactions
  void set_interaction_type_and_functional(std::string t,std::string f);
  void set_2a_interaction(std::string t,std::string f,
                          int,int,
                          VECTOR& r1,VECTOR& r2,
                          VECTOR& g1,VECTOR& g2,
                          VECTOR& m1,VECTOR& m2,
                          VECTOR& f1,VECTOR& f2,
                          double& dr1_2, double& dr2_2, double& dT_2,
                          map<std::string,double> params);

  void set_3a_interaction(std::string t,std::string f,
                          int,int,int,
                          VECTOR& r1,VECTOR& r2,VECTOR& r3,
                          VECTOR& g1,VECTOR& g2,VECTOR& g3,
                          VECTOR& m1,VECTOR& m2,VECTOR& m3,
                          VECTOR& f1,VECTOR& f2,VECTOR& f3,
                          map<std::string,double> params);

  void set_4a_interaction(std::string t,std::string f,
                          int,int,int,int,
                          VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,
                          VECTOR& g1,VECTOR& g2,VECTOR& g3,VECTOR& g4,
                          VECTOR& m1,VECTOR& m2,VECTOR& m3,VECTOR& m4,
                          VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,
                          map<std::string,double> params);

  void set_mb_interaction(std::string t,std::string f,
                          int, int*, VECTOR**, VECTOR**,VECTOR**,
                          VECTOR**,
                          double**,double**, double**,int nexcl, int* excl1, int* excl2, double* scale,
                          double**, double*,
                          vector< vector<excl_scale> >& excl_scales,
                          map<std::string,double> params);

  // 2 - fragment interactions
  void set_2f_interaction(std::string t,std::string f,
                          int,int,
                          VECTOR& r1,VECTOR& r2,VECTOR& u1,VECTOR& u2,
                          VECTOR& f1,VECTOR& f2,VECTOR& t1,VECTOR& t2,
                          map<std::string,double> params);

  double calculate(int&);
  double calculate(int,int&);
  
};



class listHamiltonian_MM{
/**
  This class represents a collection of classical interactions in a molecular (or solid-state) system.
  This is essentially a classical-mechanical Hamiltonian (although it is called a listHamiltonian - this is only for the
  uniformity with QM Hamiltonian and also for better flexibility - this way, the actual Hamiltonian can be a multi-resolution
  Hamiltonian, so each component can be tackeld by different MM level of theory)
*/


public:

    listHamiltonian_MM(){ ;; }


    vector<Hamiltonian_MM> interactions;  ///< The list of classical interaction (individual, primitive Hamiltonians)
    vector<int>     active_interactions;  ///< The list with the indices showing which of these interactions are actually active
                                          ///< Note that only "active" interactions are computed when "compute" is applied

    std::string stress_opt;int is_stress_opt;  ///< The option that selects which type of stress to compute: "at", "fr", "ml"
    MATRIX3x3 stress_at;   int is_stress_at;   ///< Total atomic stress and the status flag
    MATRIX3x3 stress_fr;   int is_stress_fr;   ///< Total fragment stress and the status flag
    MATRIX3x3 stress_ml;   int is_stress_ml;   ///< Total molecular stress and the status flag
    MATRIX3x3 hessian;     int is_hessian;     ///< Total hessian and the status flag

    // RESPA auxiliary variables
    vector<VECTOR> respa_f_fast,respa_f_medium;  ///< RESPA forces: fast and medium components
    vector<VECTOR> respa_t_fast,respa_t_medium;  ///< RESPA torques: fast and medium components
    MATRIX3x3 respa_s_fast,respa_s_medium;     
    double respa_E_fast,respa_E_medium;          ///< RESPA energies: fast and medium components


  //----------- Defined in Hamiltonian_MM_methods2.cpp ------------------
  // Interaction related functions:
  int is_new_interaction(Hamiltonian_MM&);
  void show_interactions_statistics();
  void show_interactions();

  int set_atom_types(System& syst,vector<int>& lst,ForceField& ff);
  int set_fragment_types(System& syst, vector<int>& lst,ForceField& ff);

  bool is_active(Atom&,Atom&);
  bool is_active(Atom&,Atom&,Atom&);
  bool is_active(Atom&,Atom&,Atom&,Atom&);

  void set_atom_interactions_for_atoms(System& syst,string int_type,vector<Atom>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff,int verb);
  void set_group_interactions_for_atoms(System& syst,string int_type,vector<Group>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff);

  void set_interactions_for_atoms(System& syst, boost::python::list,boost::python::list,ForceField&,int verb, int assign_rings);
  void set_interactions_for_fragments(System& syst, boost::python::list,boost::python::list,ForceField&);


  void apply_pbc_to_interactions(System& syst, int int_type,int nx,int ny,int nz);
  void set_respa_types(std::string inter_type,std::string respa_type);



};




}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra


#endif // HAMILTONIAN_MM_H
