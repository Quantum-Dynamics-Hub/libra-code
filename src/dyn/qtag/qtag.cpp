/*********************************************************************************
* Copyright (C) 2022 Matthew Dutra, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file qtag.cpp
  \brief The file implements various auxiliary functions for QTAG method
    
*/

#include "qtag.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

/// libqtag namespace
namespace libqtag{


double qtag_momentum(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff){
/**
 Returns the momentum *mom* calculated for each basis function according to p=Im(grad(psi)/psi). 
       Also returns the corresponding real component *r*, which can be used in updates of the basis phase parameter *s*.

*/


  int i,j, idof;
  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  complex<double> dzt;
  double dq, nrm;

  CMATRIX mom(ndof,ntraj); // real = mom, imag = r
  CMATRIX dz(ndof,1);
  complex<double> one(0.0, 1.0);

  for(i=0; i<ntraj; i++){
    complex<double> z(0.0,0.0);
    dz *= 0.0;

    for(j=0; j<ntraj; j++){
      complex<double> nrm(1.0, 0.0);

      for(idof=0; idof<ndof; idof++){
          dq = q.get(idof, i) - q.get(idof, j);
          double argg = (p.get(idof, j) * dq + s.get(idof, j) );
          complex<double> tmp( cos(argg), sin(argg) );
          nrm *= exp( -0.5 * alp.get(idof, j) * dq*dq ) * tmp * pow( (alp.get(idof, i)/M_PI), 0.25);
     
      }// for idof

      z += Coeff.get(j) * nrm;

      for(idof=0; idof<ndof; idof++){
        dq = q.get(idof, i) - q.get(idof, j);
        dzt = -(alp.get(idof, j) * dq - one * p.get(idof, j)) * Coeff.get(j) * nrm;
        dz.add(idof, 0, dzt);
      }// for idof
      
    }// for j

    for(idof=0; idof<ndof; idof++){
      mom.set(i, idof,  dz.get(idof)/z );
    }// for idof

  }// for i

}




}// namespace libqtag
}// namespace libdyn
}// liblibra


