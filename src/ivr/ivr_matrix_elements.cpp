/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file ivr_matrix_elements.cpp
  \brief These are the C++ re-implementation of the Fortran codes from the Ananth group.
         The original codes can be found here:
         https://github.com/AnanthGroup/SC-IVR-Code-Package    

  According to original documentation:

! This file contains subroutines that 
! compute coherent state overlaps matrix elements
! for any multidimensional  SC-IVR:
!                                                                         
! Coherent state overlap                                           
!                                                                         
! Operator B  matrix elements:                                   
!   - Forward-Backward: B = B(q,p)
!   - Forward-Forward : B = B(q,p)
!   - Husimi-IVR      : B = B(q,p)
!   - LSC-IVR         : B = B(q,p)


*/

#include "ivr.h"


/// liblibra namespace
namespace liblibra{


namespace libivr{

complex<double> CS_overlap(MATRIX& q, MATRIX& p, MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& invWidth0){
/**
  \brief Multi-dimensional  Coherent State (CS) overlap:  <p0,q0|pIn,qIn>

  \param[in] q - coordinates (Ndof x 1 matrix)
  \param[in] p - momenta (Ndof x 1 matrix)
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] invWidth0 - a Ndof x Ndof matrix with the inverse width parameters of coherent states (CS) 
      
  Return value: a complex number representing an overlap of 2 CSs

*/

  if(qIn.n_rows!=pIn.n_rows){
    cout<<"Error in CS_overlap: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }

  if(q.n_rows!=p.n_rows){
    cout<<"Error in CS_overlap: the input matrices q and p have different # of rows\n";
    exit(0);
  }

  
  int Ndof = qIn.n_rows;


  // To avoid the full-scale vector-matrix multiplication:
  MATRIX Wq(Ndof, 1), iWp(Ndof, 1);

  for(int i=0;i<Ndof;i++){
    Wq.M[i] = Width0.get(i,i) * (q.get(i) - qIn.get(i));
    iWp.M[i] = invWidth0.get(i,i) * (p.get(i) - pIn.get(i));
  }

  double pr1 = ((q - qIn).T() * Wq).M[0];
  double pr2 = ((p - pIn).T() * iWp).M[0];
  double ar = 0.5 * ((p + pIn).T() * (q - qIn)).M[0];

  return exp(-0.25*(pr1+pr2)) * complex<double>(cos(ar), sin(ar)); 

}


complex<double> mat_elt_FB_B(MATRIX& q, MATRIX& p, int opt, int lab){
/**
  \brief Operator B matrix element:  B = B(q,p) Forward-Backward

  \param[in] q - coordinates (Ndof x 1 matrix)
  \param[in] p - momenta (Ndof x 1 matrix)
  \param[in] opt - type of the operator to compute
             opt = 0:  B(q,p) = q
             opt = 1:  B(q,p) = p
  \param[in] lab - is the label indicating the obsvered mode of the system

  Here we use operators B = B(q) = q or B = B(p) = p
  Other position based operators should be changed to B = B(q)
  Example: If B = q^2-q, write Bq = q**2-q



*/
  double res = 0.0;
  if(opt==0){  res = q.get(lab); }
  else if(opt==1){  res = p.get(lab); }

  return complex<double>(res,0.0); 

}



complex<double> mat_elt_FF_B(MATRIX& q, MATRIX& p, MATRIX& qp, MATRIX& pp, MATRIX& WidthT, MATRIX& invWidthT, int opt, int lab){
/**
  \brief Operator B matrix element:  B = B(q,p) Forward-Forward

  \param[in] q - coordinates (Ndof x 1 matrix)
  \param[in] p - momenta (Ndof x 1 matrix)
  \param[in] opt - type of the operator to compute
             opt = 0:  B(q,p) = q
             opt = 1:  B(q,p) = p
  \param[in] lab - is the label indicating the obsvered mode of the system

  Here we use operators B = B(q) = q or B = B(p) = p
  Another expression must be derived for another operator of type B(q),
  see documentation for details.
  q = qt, p = pt, qp = qt', pp = pt'

*/

  if(q.n_rows!=p.n_rows){
    cout<<"Error in mat_elt_FF_B: the input matrices q and p have different # of rows\n";
    exit(0);
  }

  if(qp.n_rows!=pp.n_rows){
    cout<<"Error in mat_elt_FF_B: the input matrices qp and pp have different # of rows\n";
    exit(0);
  }

  
  int Ndof = q.n_rows;


  // To avoid the full-scale vector-matrix multiplication:
  MATRIX Wq(Ndof, 1), iWp(Ndof, 1);

  for(int i=0;i<Ndof;i++){
    Wq.M[i] = WidthT.get(i,i) * (q.get(i) - qp.get(i));
    iWp.M[i] = invWidthT.get(i,i) * (p.get(i) - pp.get(i));
  }

  double pr1 = ((q - qp).T() * Wq).M[0];
  double pr2 = ((p - pp).T() * iWp).M[0];
  double ar = 0.5 * ((p + pp).T() * (q - qp)).M[0];
  complex<double> ovlp = exp(-0.25*(pr1+pr2)) * complex<double>(cos(ar), sin(ar)); 
  complex<double> pr3;

  if(opt==0){  // B = B(q) = q
    pr3 = complex<double>(0.5 * (q.get(lab) + qp.get(lab)) , invWidthT.get(lab,lab)*(p.get(lab)-pp.get(lab)) );
  }
  else if(opt==1){  // B = B(p) = p
    pr3 = complex<double>(0.5 * (p.get(lab) + pp.get(lab)) , -invWidthT.get(lab,lab)*(q.get(lab)-qp.get(lab)) );
  }

  return pr3 * ovlp;

}


complex<double> mat_elt_HUS_B(MATRIX& q, MATRIX& p, int opt, int lab){
/**
  \brief Operator B matrix element:  B = B(q,p) Husimi - same as Forward-Backward

  \param[in] q - coordinates (Ndof x 1 matrix)
  \param[in] p - momenta (Ndof x 1 matrix)
  \param[in] opt - type of the operator to compute
             opt = 0:  B(q,p) = q
             opt = 1:  B(q,p) = p
  \param[in] lab - is the label indicating the obsvered mode of the system

  Here we use operators B = B(q) = q or B = B(p) = p
  Another expression must be derived with another operator of type B=B(q,p) 
  See documentation for more details

*/

  double res = 0.0;
  if(opt==0){  res = q.get(lab); }
  else if(opt==1){  res = p.get(lab); }

  return complex<double>(res,0.0); 

}

complex<double> mat_elt_LSC_B(MATRIX& q, MATRIX& p, int opt, int lab){
/**
  \brief Operator B matrix element:  B = B(q,p) LSC - same as Forward-Backward

  \param[in] q - coordinates (Ndof x 1 matrix)
  \param[in] p - momenta (Ndof x 1 matrix)
  \param[in] opt - type of the operator to compute
             opt = 0:  B(q,p) = q
             opt = 1:  B(q,p) = p
  \param[in] lab - is the label indicating the obsvered mode of the system

  Here we use operators B = B(q) = q or B = B(p) = p
  Another expression must be derived with another operator of type B=B(q,p) 
  See documentation for more details

*/

  double res = 0.0;
  if(opt==0){  res = q.get(lab); }
  else if(opt==1){  res = p.get(lab); }

  return complex<double>(res,0.0); 

}




}/// namespace libivr

}/// namespace liblibra



