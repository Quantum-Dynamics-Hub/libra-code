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
#if defined(USING_PCH)
#include "../pch.h"
#else
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#endif 

#include "mEigen.h"

/// liblibra namespace
namespace liblibra{

using namespace Eigen;
using namespace std;
using namespace liblinalg;

/// libmeigen namespace
namespace libmeigen{

double det(MATRIX& A){
/** Compute the determinat of a real-valued square matrix 
    Wrapper of Eigen3 MatrixXd .determinant method
*/
 
  if(A.n_cols!=A.n_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.n_cols;
  int i,j;

  MatrixXd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  return a.determinant();


}

double FullPivLU_det(MATRIX& A){
/** Compute the determinat of a real-valued square matrix  
    Alternative version, based on LU decompostion
    Wrapper of the Eigen::FullPivLU<MatrixXd> class
*/

  if(A.n_cols!=A.n_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.n_cols;
  int i,j;

  MatrixXd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  Eigen::FullPivLU<MatrixXd> lu(a);

  return lu.determinant();


}



complex<double> det(CMATRIX& A){
/** Compute the determinat of a complex-valued square matrix 
    Wrapper of Eigen3 Matrixcd .determinant method
*/


  if(A.n_cols!=A.n_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i

  return a.determinant();

}

complex<double> FullPivLU_det(CMATRIX& A){
/** Compute the determinat of a complex-valued square matrix  
    Alternative version, based on LU decompostion
    Wrapper of the Eigen::FullPivLU<MatrixXcd> class
*/

  if(A.n_cols!=A.n_rows){
    std::cout<<"Error in det(MATRIX): Can not compute a determinant of non-square matrix\n";
    exit(0);
  }

  int N = A.n_cols;
  int i,j;

  MatrixXcd a(N,N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      a(i,j) = A.M[i*N+j];
    }// for j
  }// for i


  Eigen::FullPivLU<MatrixXcd> lu(a);

  return lu.determinant();


}



}// namespace libmeigen
}// namespace liblibra
