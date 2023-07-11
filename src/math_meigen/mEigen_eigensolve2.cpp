/*********************************************************************************
* Copyright (C) 2015-2019 Alexey V. Akimov
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




void solve_eigen(CMATRIX* H, CMATRIX* E, CMATRIX* C, int symm){
/** Solve H * C = C * E
 i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
 C[i][j] - is the weight of j-th AO in i-th MO

 ACHTUNG:    EigenSolver<MatrixXd> solution(A); - messes up the ordering of eigenvalues
 GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B); does not!, so 
 here we actually specialize the generalized eigensolver.
*/

  int N_bas = H->n_cols;
  CMATRIX* S; S = new CMATRIX(N_bas,N_bas); S->diag(N_bas, 1.0);

  solve_eigen(H, S, E, C, symm);  

  delete S;
}


void solve_eigen(CMATRIX& H, CMATRIX& E, CMATRIX& C, int symm){
/** Solve H * C = C * E */
  solve_eigen(&H, &E, &C, symm);  
}





void solve_eigen(MATRIX* H, CMATRIX* E, CMATRIX* C, int symm){
/** Solve H * C = C * E
 i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
 C[i][j] - is the weight of j-th AO in i-th MO

 ACHTUNG:    EigenSolver<MatrixXd> solution(A); - messes up the ordering of eigenvalues
 GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B); does not!, so 
 here we actually specialize the generalized eigensolver.
*/

  int N_bas = H->n_cols;
  MATRIX* S; S = new MATRIX(N_bas,N_bas); S->diag(N_bas, 1.0);

  solve_eigen(H, S, E, C, symm);  

  delete S;
}


void solve_eigen(MATRIX& H, CMATRIX& E, CMATRIX& C, int symm){
/** Solve H * C = C * E */
  solve_eigen(&H, &E, &C, symm);  
}





void solve_eigen(MATRIX* H, MATRIX* E, MATRIX* C, int symm){
/** Solve H * C = C * E
 i-th column of C contains i-th MO (coefficients of expansion in terms of AOs)
 C[i][j] - is the weight of j-th AO in i-th MO

 ACHTUNG:    EigenSolver<MatrixXd> solution(A); - messes up the ordering of eigenvalues
 GeneralizedSelfAdjointEigenSolver<MatrixXd> solution(A,B); does not!, so 
 here we actually specialize the generalized eigensolver.
*/

  int N_bas = H->n_cols;
  MATRIX* S; S = new MATRIX(N_bas,N_bas); S->diag(N_bas, 1.0);

  solve_eigen(H, S, E, C, symm);  

  delete S;
}


void solve_eigen(MATRIX& H, MATRIX& E, MATRIX& C, int symm){
/** Solve H * C = C * E */
  solve_eigen(&H, &E, &C, symm);  
}




}// namespace libmeigen
}// namespace liblibra
