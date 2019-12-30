/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Wfcgrid2_ColbertMiller.cpp
  \brief The file implements the TD-SE integrators based on the Colbert-Miller DVR
  method     
*/

#include "Wfcgrid2.h"
#include "../../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{

using namespace libmeigen;


int points_on_same_line(int idof, vector<int>& pt1, vector<int>& pt2){
/**

 This function checks if the two multidimensional points are on
 the same line, which is the axis along the specified direction idof, 

 Returns:
 0 (false) - if the projections along any but `idof` axis are distinct

 1 (true) - if the projections along all but `idof` axes are the same

*/

  int dof;
  int res = 1;

  int ndof = pt1.size();

  for(dof=0; dof<ndof && res==1; dof++){

    if(dof!=idof){

      if(pt1[dof] != pt2[dof]){  

        res = 0;

      }
    }// exclude a specified dof axiss 

  }// for all dofs

  return res;
}


vector<CMATRIX> Wfcgrid2::T_PSI(vector<CMATRIX>& inp_psi, vector<int>& bc_type, vector<double>& mass, complex<double> scaling){
/**
  \brief Application of the kinetic energy operator to the wavefunction, according to Colbert-Miller formula
  
  T_psi = scaling * T * inp_psi

  bc_type - the boundary condition along each dof:

  0 - finite (the most general case) - (a, b) - use rmin and rmax for that dof, in this case
  1 - (-inf, +inf)
  2 - (0, +inf)

*/

  int idof, ipt1, ipt2, i, j, istate;
  vector<int> pt1(ndof, 0);
  vector<int> pt2(ndof, 0);  


  vector< double > A(ndof, 0.0);

  ///========= Common pre-factors: different for each BC type ==========
  for(idof=0; idof<ndof; idof++){  

    if(bc_type[idof]==0){
      /// ==========(rmin, rmax) ==========
      double L = rmax[idof] - rmin[idof];
      A[idof] = 0.25 * M_PI * M_PI / (mass[idof] * L * L);

    }

    else if(bc_type[idof]==1 || bc_type[idof]==2){
      /// ==========(-inf, +inf) ==========
      A[idof] = 0.5 / (mass[idof] * dr[idof]*dr[idof]);

    }

  }


  vector<CMATRIX> res(Npts, CMATRIX(nstates, 1) );
  CMATRIX T(nstates, nstates);

  double T_ij = 0.0;

  for(ipt1=0; ipt1<Npts; ipt1++){

    pt1 = gmap[ipt1];

    for(ipt2=0; ipt2<Npts; ipt2++){

      pt2 = gmap[ipt2];

      for(idof=0; idof<ndof; idof++){

        if( points_on_same_line(idof, pt1, pt2)==1 ){

          i = pt1[idof];
          j = pt2[idof];


          if(bc_type[idof]==0){
          /// ==========(rmin, rmax) ==========

            if(i==j){  

              if(i>0){
                double si = sin(M_PI*i/npts[idof]);
                T_ij = (2.0 * npts[idof] * npts[idof] + 1.0)/3.0;
                T_ij -= 1.0/(si*si);
              }
              else{  T_ij = 0.0; }

            }
            else{   
              if(i>0 && j>0){

                double si = sin(0.5*M_PI*(i-j)/npts[idof]);
                T_ij = 1.0/(si*si);

                si = sin(0.5*M_PI*(i+j)/npts[idof]);
                T_ij -= 1.0/(si*si);

                T_ij *= pow(-1.0, i-j);
              }
              else{ T_ij = 0.0; }
            }

            T_ij *= A[idof];

          }// bc_type == 0


          else if(bc_type[idof]==1){
          /// ==========(-inf, +inf) ==========          

            if(i==j){  T_ij = M_PI * M_PI/3.0;  }
            else{   T_ij = 2.0/pow((i-j),2);    }

            T_ij *= pow(-1.0, i-j);
            T_ij *= A[idof];

          }// bc_type == 1  

          else if(bc_type[idof]==2){
          /// ==========(0, +inf) ==========

            if(i==j){  
              if(i>0){
                T_ij = M_PI * M_PI/3.0 - 0.5/double(i*i); 
              }
            }
            else{   T_ij = 2.0/double((i-j)*(i-j)) - 2.0/double((i+j)*(i+j)) ;    }

            T_ij *= pow(-1.0, i-j);
            T_ij *= A[idof];

          }// bc_type == 2



          for(istate=0; istate<nstates; istate++){
            T.set(istate, istate, complex<double>(T_ij, 0.0));
          }
        
          res[ipt1] += T * inp_psi[ipt2];

        }

      }// for idof 

    }// for ipt2
  }// for ipt1


  // Scale
  for(ipt1=0; ipt1<Npts; ipt1++){
    res[ipt1] *= scaling;
  }


  return res;

}


vector<CMATRIX> Wfcgrid2::T_PSI_adi(vector<int>& bc_type, vector<double>& mass, complex<double> scaling){

  return T_PSI(PSI_adi, bc_type, mass, scaling);

}

vector<CMATRIX> Wfcgrid2::T_PSI_dia(vector<int>& bc_type, vector<double>& mass, complex<double> scaling){

  return T_PSI(PSI_dia, bc_type, mass, scaling);

}





}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

