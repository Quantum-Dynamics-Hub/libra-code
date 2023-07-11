/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file AO.h
  \brief The file describes: a) the AO class that represents atomic orbtals as a linear combination of 
  Gaussian primitives; b) related functions
    
*/

#ifndef AO_H
#define AO_H

#include "PrimitiveG.h"

#include "../math_linalg/liblinalg.h"
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;


namespace libqobjects{


class AO{
/**
  Objects of this class are atomic orbitals
  AO = linear combination of primitive gaussians (contraction)
*/

public:

  // Members
  std::string element;              int is_element;       ///< name of chemical element on which this AO is localized (if so)
  std::string ao_shell;             int is_ao_shell;      ///< name of the AO shell e.g. 1s, 2p, 3p, 3d, etc.
  std::string ao_shell_type;        int is_ao_shell_type; ///< type of the AO shell e.g. s, p, d,  etc.
  std::string ao_name;              int is_ao_name;       ///< name of the AO e.g. 1s, 2px, 3dxy,  etc.  int at_indx;

  // Angular momentum is defined for AOs, but each AO may have many Gaussians with different
  // alpha constants
  int x_exp;                        int is_x_exp;         ///< angular momentum quantum number (x projection) - if all primitives have it same
  int y_exp;                        int is_y_exp;         ///< angular momentum quantum number (y projection) - if all primitives have it same
  int z_exp;                        int is_z_exp;         ///< angular momentum quantum number (z projection) - if all primitives have it same
  int expansion_size;               int is_expansion_size;///< The number of Primitive Gaussians in the contraction
  vector<PrimitiveG> primitives;    int is_primitives;    ///< Primitive Gaussians that constitute this AO
  vector<double> coefficients;      int is_coefficients;  ///< Contraction coefficients: these coefficients correspond to normalized primitive Gaussians



  //----------------- Function members --------------------
  AO();  ///   c-tor
  AO(const AO&); /// cc-tor
  AO& operator=(const AO&);

  // Construction and print  
  void init();
  void clear();
  void add_primitive(double c,PrimitiveG g);
  void show_info();

  // Computations
  double compute(VECTOR&);
  double norm2();
  double norm1();
  double normalization_factor();
  void normalize();


  // Transformations
  void shift_position(VECTOR);
  void set_position(VECTOR);
  void shift_position_const_ref(const VECTOR&);
  void set_position_const_ref(const VECTOR&);




  friend int operator == (const AO& g1, const AO& g2){
   /*
    int res = ((g1.x_exp==g2.x_exp) && (g1.y_exp==g2.y_exp) && (g1.z_exp==g2.z_exp)
              && (g1.expansion_size==g2.expansion_size) && (g1.ao_name==g2.ao_name)
              && (g1.ao_shell == g2.ao_shell) && (g1.ao_shell_type == g2.ao_shell_type) 
              && (g1.element == g2.element) 
              );
    if(res){ 
      for(int i=0;i<g1.expansion_size;i++){  res *= (g1.coefficients[i]==g2.coefficients[i]); }
    }
    if(res){ 
      for(int i=0;i<g1.expansion_size;i++){  res *= (g1.primitives[i]==g2.primitives[i]); }
    }
    return res;
   */
   return (&g1 == &g2);
  }


};


//====================== Overlaps ================================
// Versions with references
double gaussian_overlap
( AO& AOa, AO& AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
);
double gaussian_overlap( AO& AOa, AO& AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB );
boost::python::list gaussian_overlap( AO& AOa, AO& AOb,int is_normalize, int is_derivs);
double gaussian_overlap(AO& AOa, AO& AOb,int is_normalize);
double gaussian_overlap(AO& AOa, AO& AOb);


// Versions with pointers - only for C++
double gaussian_overlap
( AO* AOa, AO* AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
);
double gaussian_overlap( AO* AOa, AO* AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB );
double gaussian_overlap(AO* AOa, AO* AOb,int is_normalize);
double gaussian_overlap(AO* AOa, AO* AOb);


//====================== Moments ================================
// Versions with references
double gaussian_moment
( AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB, vector<double*>& auxd,int n_aux
);
double gaussian_moment( AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdG, VECTOR& dIdB );
boost::python::list gaussian_moment( AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize, int is_derivs);
double gaussian_moment(AO& AOa, PrimitiveG& G, AO& AOb,int is_normalize);
double gaussian_moment(AO& AOa, PrimitiveG& G, AO& AOb);


// Versions with pointers - only for C++
double gaussian_moment
( AO* AOa, PrimitiveG& G, AO* AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB, vector<double*>& auxd,int n_aux
);
double gaussian_moment( AO* AOa, PrimitiveG& G, AO* AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdR, VECTOR& dIdB );
double gaussian_moment(AO* AOa, PrimitiveG& G, AO* AOb,int is_normalize);
double gaussian_moment(AO* AOa, PrimitiveG& G, AO* AOb);


//====================== Pseudopotentials ================================
// Versions with references
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb, int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux      );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb,int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB  );
boost::python::list pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                                AO& AOa, AO& AOb, int is_normalize, int is_derivs );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO& AOa, AO& AOb,int is_normalize );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R, AO& AOa, AO& AOb );


// Versions with pointers - only for C++
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb, int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB,
                   vector<double*>& auxd,int n_aux      );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb,int is_normalize, 
                   int is_derivs, VECTOR& dIdR, VECTOR& dIdA, VECTOR& dIdB  );
boost::python::list pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                                AO* AOa, AO* AOb, int is_normalize, int is_derivs );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R,
                   AO* AOa, AO* AOb,int is_normalize );
double pseudopot02(double C0, double C2, double alp, const VECTOR& R, AO* AOa, AO* AOb );



//====================== Multipoles ================================
// Versions with references
VECTOR transition_dipole_moment
( AO& AOa, AO& AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
);
VECTOR transition_dipole_moment
( AO& AOa, AO& AOb,int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
);
boost::python::list transition_dipole_moment( AO& AOa, AO& AOb, int is_normalize,int is_derivs);
VECTOR transition_dipole_moment( AO& AOa, AO& AOb, int is_normalize);
VECTOR transition_dipole_moment( AO& AOa, AO& AOb);


// Versions with pointers - only for C++
VECTOR transition_dipole_moment
( AO* AOa, AO* AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
);
VECTOR transition_dipole_moment
( AO* AOa, AO* AOb,int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
);
boost::python::list transition_dipole_moment( AO* AOa, AO* AOb, int is_normalize,int is_derivs);
VECTOR transition_dipole_moment( AO* AOa, AO* AOb, int is_normalize);
VECTOR transition_dipole_moment( AO* AOa, AO* AOb);



//====================== Derivative coupling integrals ================================
// Versions with references
VECTOR derivative_coupling_integral
( AO& AOa, AO& AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
);
VECTOR derivative_coupling_integral
( AO& AOa, AO& AOb,int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
);
boost::python::list derivative_coupling_integral( AO& AOa, AO& AOb, int is_normalize,int is_derivs);
VECTOR derivative_coupling_integral( AO& AOa, AO& AOb, int is_normalize);
VECTOR derivative_coupling_integral( AO& AOa, AO& AOb);


// Versions with pointers - only for C++
VECTOR derivative_coupling_integral
( AO* AOa, AO* AOb,
  int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB,
  vector<double*>& auxd,int n_aux
);
VECTOR derivative_coupling_integral
( AO* AOa, AO* AOb,int is_normalize,int is_derivs, MATRIX3x3& dMdA, MATRIX3x3& dMdB
);
boost::python::list derivative_coupling_integral( AO* AOa, AO* AOb, int is_normalize,int is_derivs);
VECTOR derivative_coupling_integral( AO* AOa, AO* AOb, int is_normalize);
VECTOR derivative_coupling_integral( AO* AOa, AO* AOb);



//====================== Kinetic integral ================================
// Versions with references
double kinetic_integral
( AO& AOa, AO& AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
);
double kinetic_integral( AO& AOa, AO& AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB );
boost::python::list kinetic_integral( AO& AOa, AO& AOb,int is_normalize, int is_derivs);
double kinetic_integral(AO& AOa, AO& AOb,int is_normalize);
double kinetic_integral(AO& AOa, AO& AOb);


// Versions with pointers - only for C++
double kinetic_integral
( AO* AOa, AO* AOb,int is_normalize, int is_derivs,
  VECTOR& dIdA, VECTOR& dIdB, vector<double*>& auxd,int n_aux
);
double kinetic_integral( AO* AOa, AO* AOb,int is_normalize, int is_derivs, VECTOR& dIdA, VECTOR& dIdB );
double kinetic_integral(AO* AOa, AO* AOb,int is_normalize);
double kinetic_integral(AO* AOa, AO* AOb);


//====================== Nuclear Attraction Integral ================================
// Versions with references
double nuclear_attraction_integral
( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);
double nuclear_attraction_integral
( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC
);
boost::python::list nuclear_attraction_integral( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize, int is_derivs );
double nuclear_attraction_integral( AO& AOa, AO& AOb, VECTOR& Rc, int is_normalize );
double nuclear_attraction_integral( AO& AOa, AO& AOb, VECTOR& Rc );


// Versions with pointers - only for C++
double nuclear_attraction_integral
( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);
double nuclear_attraction_integral
( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC
);
boost::python::list nuclear_attraction_integral( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize, int is_derivs );
double nuclear_attraction_integral( AO* AOa, AO* AOb, VECTOR& Rc, int is_normalize );
double nuclear_attraction_integral( AO* AOa, AO* AOb, VECTOR& Rc );


//====================== Electron Repulsion Integral ================================
// Versions with references
double electron_repulsion_integral
( AO& AOa, AO& AOb, AO& AOc, AO& AOd,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);
double electron_repulsion_integral
( AO& AOa, AO& AOb, AO& AOc, AO& AOd,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
);
boost::python::list electron_repulsion_integral( AO& AOa, AO& AOb, AO& AOc, AO& AOd, int is_normalize, int is_derivs);
double electron_repulsion_integral( AO& AOa, AO& AOb, AO& AOc, AO& AOd,int is_normalize);
double electron_repulsion_integral( AO& AOa, AO& AOb, AO& AOc, AO& AOd);


// Versions with pointers - only for C++
double electron_repulsion_integral
( AO* AOa, AO* AOb, AO* AOc, AO* AOd,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD,
  vector<double*>& aux,int n_aux,vector<VECTOR*>& auxv,int n_auxv
);
double electron_repulsion_integral
( AO* AOa, AO* AOb, AO* AOc, AO* AOd,
  int is_normalize, 
  int is_derivs,  VECTOR& DA,VECTOR& DB,VECTOR& DC,VECTOR& DD
);
boost::python::list electron_repulsion_integral( AO* AOa, AO* AOb, AO* AOc, AO* AOd, int is_normalize, int is_derivs);
double electron_repulsion_integral( AO* AOa, AO* AOb, AO* AOc, AO* AOd,int is_normalize);
double electron_repulsion_integral( AO* AOa, AO* AOb, AO* AOc, AO* AOd);



typedef std::vector<AO> AOList; ///< This is the data type for representing vector of AO objects



}// namespace libqobjects
}// namespace liblibra



#endif // AO_H
