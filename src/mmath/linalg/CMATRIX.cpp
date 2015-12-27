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

#include "CMATRIX.h"
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <iostream>

using namespace std;


namespace libmmath{
namespace liblinalg{


void CMATRIX::set(int indx,double value1, double value2){

    M[indx] = complex<double>(value1,value2);
}

void CMATRIX::set(int indx, complex<double> value1){

    M[indx] = value1;
}


complex<double> CMATRIX::get(int indx){

    return M[indx];
}

void CMATRIX::set(int row,int col,double value1,double value2){

   M[row*n_cols+col] = complex<double>(value1,value2);
}

void CMATRIX::set(int row,int col, complex<double> value1){

   M[row*n_cols+col] = value1;
}


complex<double> CMATRIX::get(int row,int col){

  return M[row*n_cols+col];
}


CMATRIX::CMATRIX(vector<vector<double> >& re_part,vector<vector<double> >& im_part){
/*****************************************************************
  Constructor: creates a CMATRIX from 2 2D-arrays - real and imaginary parts
*****************************************************************/
  if(re_part.size()!=im_part.size()){ 
    cout<<"Error in CMATRIX constructor: y-dimensions(num of rows) of the real and imaginary arrays are not equal\n"; exit(0); 
  }
  if(re_part[0].size()!=im_part[0].size()){
    cout<<"Error in CMATRIX constructor: x-dimensions(num of cols) of the real and imaginary arrays are not equal\n"; exit(0);
  }
  n_rows = re_part.size();
  n_cols = re_part[0].size();
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>(re_part[i][j],im_part[i][j]); n++;
    }
  }
 
}

CMATRIX::CMATRIX(MATRIX& re_part){

  n_rows = re_part.num_of_rows;
  n_cols = re_part.num_of_cols;
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>(re_part.get(i,j), 0.0); n++;
    }
  }


}

CMATRIX::CMATRIX(MATRIX& re_part,MATRIX& im_part){

  if(re_part.num_of_rows!=im_part.num_of_rows){ 
    cout<<"Error in CMATRIX constructor: num_of_rows of the real and imaginary matrices are not equal\n"; exit(0); 
  }
  if(re_part.num_of_cols!=im_part.num_of_cols){
    cout<<"Error in CMATRIX constructor: num_of_cols of the real and imaginary arrays are not equal\n"; exit(0);
  }
  n_rows = re_part.num_of_rows;
  n_cols = re_part.num_of_cols;
  n_elts = n_rows * n_cols;

  M = new complex<double>[n_elts];
  int n = 0;
  for(int i=0;i<n_rows;i++){ 
    for(int j=0;j<n_cols;j++){
      M[n] = complex<double>(re_part.get(i,j),im_part.get(i,j)); n++;
    }
  }

}




CMATRIX::CMATRIX(const CMATRIX& obj){
  n_rows = obj.n_rows;
  n_cols = obj.n_cols;
  n_elts = obj.n_elts;

  M = new complex<double>[n_elts];
  for(int i=0;i<n_elts;i++){ M[i] = obj.M[i];  }

//  cout<<"Copy constructor\n";
}

CMATRIX CMATRIX::operator-(){
  CMATRIX tmp(n_rows,n_cols);
  for(int i=0;i<n_elts;i++){ tmp.M[i]=-M[i]; }
  return tmp;
}

CMATRIX CMATRIX::operator*(const CMATRIX& ob){
// (n_rows x ob.n_cols) = (n_rows x n_cols) * (ob.n_rows * ob.n_cols)
// n_cols must be equal to ob.n_rows
  if(n_cols!=ob.n_rows){ 
    std::cout<<"Matrix multiplication error: Dimensions of operands must match\n";
    std::cout<<"You try to muplitpy CMATRIX "<<n_rows<<" by "<<n_cols<<" and the CMATRIX "
             <<ob.n_rows<<" by "<<ob.n_cols<<"\n";
    std::cout<<"Exiting...\n";
    exit(0);  
  }
  else{
    int n=ob.n_cols;
    int kn; // k*n
    int rn; // row*n
    int rncols; // row*n_cols
   
    CMATRIX Temp(n_rows,n);

    rn = 0;
    rncols = 0;
    for(int row=0;row<n_rows;row++){
      for(int col=0;col<n;col++){
        complex<double> d(0.0,0.0);
        /* Standard formulation
        for(int k=0;k<n_cols;k++){ d+=M[row*n_cols+k]*ob.M[k*n+col];  }
        Temp.M[row*n+col] = d;
        */

        // More efficient formulation
        kn = 0;
        for(int k=0;k<n_cols;k++){ d+=M[rncols+k]*ob.M[kn+col]; kn += n; }
        Temp.M[rn+col] = d;
      }
      rn += n;
      rncols += n_cols;
    }
    return Temp;
  }
}


void CMATRIX::dot(const CMATRIX& ob1,const CMATRIX& ob2){
// Direct product of two matrices - element-wise multiplication
// Dimensions of ob1 and ob2 must be equal

  if(ob1.n_cols!=ob2.n_cols){
    std::cout<<"Error in direct CMATRIX multiplication: Number of columns of multiplying matrices must be equal\n";
    std::cout<<"Exiting...\n";
    exit(0);     
  }

  if(ob1.n_rows!=ob2.n_rows){
    std::cout<<"Error in direct CMATRIX multiplication: Number of rows of multiplying matrices must be equal\n";
    std::cout<<"Exiting...\n";
    exit(0);     
  }

  if(n_cols!=ob2.n_cols){
    std::cout<<"Error in direct CMATRIX multiplication: Number of columns of the multiplying CMATRIX and target CMATRIX must be equal\n";
    std::cout<<"Exiting...\n";
    exit(0);     
  }

  if(n_rows!=ob2.n_rows){
    std::cout<<"Error in direct CMATRIX multiplication: Number of rows of the multiplying CMATRIX and target CMATRIX must be equal\n";
    std::cout<<"Exiting...\n";
    exit(0);     
  }


  for(int row=0;row<n_rows;row++){
    for(int col=0;col<n_cols;col++){
      M[row*n_cols+col] = ob1.M[row*n_cols+col] * ob2.M[row*n_cols+col];
    }
  }
}


CMATRIX CMATRIX::operator+(const CMATRIX& ob){ 
  CMATRIX Temp(n_rows,n_cols);
  for(int i=0;i<n_elts;i++) {Temp.M[i]=M[i]+ob.M[i];}
  return Temp;
}

CMATRIX CMATRIX::operator-(const CMATRIX& ob){
  CMATRIX Temp(n_rows,n_cols);
  for(int i=0;i<n_elts;i++) {Temp.M[i]=M[i]-ob.M[i];}
  return Temp;
}

void CMATRIX::operator+=(const CMATRIX& ob){
  for(int i=0;i<n_elts;i++) {M[i]+=ob.M[i];}
}

void CMATRIX::operator-=(const CMATRIX& ob){
  for(int i=0;i<n_elts;i++) {M[i]-=ob.M[i];}
}

void CMATRIX::operator*=(const double& f){
  for(int i=0;i<n_elts;i++) {M[i]*=f;}
}

void CMATRIX::operator*=(const complex<double>& f){
  for(int i=0;i<n_elts;i++) {M[i]*=f;}
}


void CMATRIX::operator*=(const CMATRIX& ob){
// (n_rows x ob.n_cols) = (n_rows x n_cols) * (ob.n_rows * ob.n_cols)
// n_cols must be equal to ob.n_rows
  if(n_cols!=ob.n_rows){ 
    std::cout<<"Matrix multiplication error: Dimensions of operands must match\n";
    std::cout<<"You try to muplitpy CMATRIX "<<n_rows<<" by "<<n_cols<<" and the CMATRIX "
             <<ob.n_rows<<" by "<<ob.n_cols<<"\n";
    std::cout<<"Exiting...\n";
    exit(0);  
  }
  else{
    int n=ob.n_cols;
    // Counters
    int rncols; // row*n_cols
    int kn;     // k*n
    int rn;     // row*n
    complex<double> *TM;
    TM = new complex<double>[n_rows*n];
    //CMATRIX Temp(n_rows,n);

    rn = 0;
    rncols = 0;
    for(int row=0;row<n_rows;row++){
      for(int col=0;col<n;col++){
        complex<double> d(0.0,0.0);

        kn = 0;
        for(int k=0;k<n_cols;k++){ d+=M[rncols+k]*ob.M[kn+col]; kn += n; }
        TM[rn+col] = d;
      }
      rn += n;
      rncols += n_cols;

    }

    for(int i=0;i<n_elts;i++){ M[i] = TM[i]; }

    delete [] TM;
  }
  
}


CMATRIX CMATRIX::operator/(double num){ 
  CMATRIX m(n_rows,n_cols);
  for(int i=0;i<n_elts;i++){  m.M[i] = M[i]/num;  }
  return m;
}

CMATRIX CMATRIX::operator/(complex<double> num){
  CMATRIX m(n_rows,n_cols);
  for(int i=0;i<n_elts;i++){  m.M[i] = M[i]/num;  }
  return m;
}

CMATRIX& CMATRIX::operator=(const CMATRIX& ob){

  //if(this == &ob){ return *this; }
  //else{

  n_rows = ob.n_rows;
  n_cols = ob.n_cols;
  n_elts = ob.n_elts;

  // Lets comment below 2 lines - it may be not completely correct to do it this way
  // but assignment means there is a left-hand side already allocated - there is no
  // reason to delete storage and reallocate it again - just copy data staight ahead
  // This gives speedup of factor 2

  //delete [] M;  // This is very important and not very obvious - without this line - there will be a memory leak
  //M = new complex<double>[n_elts];
  
  //for(int i=0;i<n_elts;i++){ M[i] = ob.M[i];  }

  memcpy(M,ob.M,sizeof(complex<double>)*n_elts);  // this is slightly more efficient version than above
  
  return *this;


  //}

}

CMATRIX CMATRIX::operator=(double num){
  for(int i=0;i<n_elts;i++){ M[i] = num;  }
  return *this;
}

CMATRIX CMATRIX::operator=(complex<double> num){
  for(int i=0;i<n_elts;i++){ M[i] = num;  }
  return *this;
}

CMATRIX operator*(const double& f,  const CMATRIX& m1){
  CMATRIX m(m1.n_rows,m1.n_cols);
  for(int i=0;i<m1.n_elts;i++){  m.M[i]=m1.M[i]*f;  }
  return m;
}

CMATRIX operator*(const CMATRIX &m1, const double  &f){
  CMATRIX m(m1.n_rows,m1.n_cols);
  for(int i=0;i<m1.n_elts;i++){  m.M[i]=m1.M[i]*f;  }
  return m;
}

CMATRIX operator*(const complex<double>& f,  const CMATRIX& m1){
  CMATRIX m(m1.n_rows,m1.n_cols);
  for(int i=0;i<m1.n_elts;i++){  m.M[i]=m1.M[i]*f;  }
  return m;
}

CMATRIX operator*(const CMATRIX &m1, const complex<double>  &f){
  CMATRIX m(m1.n_rows,m1.n_cols);
  for(int i=0;i<m1.n_elts;i++){  m.M[i]=m1.M[i]*f;  }
  return m;
}



int operator ==(const CMATRIX& m1,const CMATRIX& m2){
        int res=1;
        for(int i=0;i<m1.n_elts;i++) {if(m1.M[i]!=m2.M[i]) {res=0;} else;}
        return res;
}


ostream& operator<<(ostream &strm,CMATRIX ob){
  strm.setf(ios::showpoint);
  for(int i=0;i<ob.n_rows;i++){
    for(int j=0;j<ob.n_cols;j++){
      strm.precision(8);
      strm.width(10);
      strm<<left<<ob.M[i*ob.n_cols+j]<<"  ";
    }
    strm<<endl;
  }
  return strm;
}
istream& operator>>(istream& strm,CMATRIX &ob){
//     Do not defined for general case       !!!      
  return strm;
}

void CMATRIX::show_matrix(){

  std::cout.setf(ios::showpoint);
  for(int i=0;i<n_rows;i++){
    for(int j=0;j<n_cols;j++){
      std::cout.precision(8);
      std::cout.width(10);
      std::cout<<left<<M[i*n_cols+j]<<"  ";
    }
    std::cout<<std::endl;
  }

}




CMATRIX CMATRIX::conj(){
  CMATRIX m(n_rows,n_cols);
  for(int i=0;i<m.n_elts;i++){ m.M[i] = std::conj(M[i]); }
  return m;
}

CMATRIX CMATRIX::T(){
  CMATRIX m(n_cols,n_rows);
  for(int i=0;i<n_rows;i++){
    for(int j=0;j<n_cols;j++){
      m.M[j*n_rows+i] = M[i*n_cols+j];
    }
  }
  return m;
}

CMATRIX CMATRIX::H(){
  CMATRIX m(n_cols,n_rows);
  for(int i=0;i<n_rows;i++){
    for(int j=0;j<n_cols;j++){
      m.M[j*n_rows+i] = std::conj(M[i*n_cols+j]);
    }
  }
  return m;
}

void CMATRIX::load_identity(){
  for(int i=0;i<n_elts;i++){ M[i] = complex<double>(0.0,0.0); }
  for(i=0;i<n_rows;i++){ M[i*n_cols+i] = complex<double>(1.0,0.0); }
}


CMATRIX CMATRIX::col(int i){
// takes given column and makes it n x 1 CMATRIX
  CMATRIX tmp(n_rows,1);
  for(int j=0;j<n_rows;j++){ tmp.M[j] = M[j*n_cols+i]; }
  return tmp;
}

CMATRIX CMATRIX::row(int i){
// takes given row and makes it 1 x n CMATRIX
  CMATRIX tmp(1,n_cols);
  for(int j=0;j<n_cols;j++){ tmp.M[j] = M[i*n_cols+j]; }
  return tmp;
}



void CMATRIX::max_nondiagonal(int& row,int& col){
  double maxeps = norm(M[1]); row = 0; col = 1;
  double eps;
  for(int r=0;r<n_rows;r++){
    for(int c=r+1;c<n_cols;c++){
      eps = norm(M[r*n_cols+c]);
      if(eps>=maxeps){ row = r; col = c; maxeps = eps; }
    }
  }  
}
/*
void CMATRIX::inverse(CMATRIX& inv,double EPS,int max_num_iter,int is_cycle,int alg){
// Based on: EVECT * EVAL = M * EVECT =>  M * (EVECT* EVAL^-1 * EVECT^-1)
// But EVECT^-1 = EVECT^T => M^-1 = (EVECT* EVAL^-1 * EVECT^T)
  if(n_rows==n_cols){
    CMATRIX eval(n_rows,n_cols);
    CMATRIX evec(n_rows,n_cols);
    CMATRIX einv(n_rows,n_cols);

    eigen(eval,evec,EPS,max_num_iter,is_cycle,alg);
    for(int i=0;i<n_rows;i++){ einv.M[i*n_cols+i] = 1.0/eval.M[i*n_cols+i]; }

    inv = evec*einv*(evec.T());
  }  
  else{ cout<<"Warning: in CMATRIX::inverse - CMATRIX is not square\n"; }
}
*/


void CMATRIX::eigen0(CMATRIX& EVAL, CMATRIX& EVECT,double EPS,int max_num_iter,int is_cycle,int alg) {
// Description: Jacobi Eigenvalue Solver - only for complex hermitian CMATRIX!
// EVECT * EVAL  =  M * EVECT
// V = P^T
// EVAL = V_M V_{M-1} ... V_0 * M * V_0^T * V_1^T ... V_M^T = Q^T * M * Q
// EVECT = Q = V_0^T * V_1^T ... V_M^T
// http://coderov.net/vma/140-eigenvalues/862-direct-method-of-rotation.html  <- this is strange, so use
// http://en.wikipedia.org/wiki/Jacobi_method_for_complex_Hermitian_matrices
// Note: Wikipedia source contains an error: CMATRIX element for m=q and n=p should be changed from
// exp(-i*teta1)*cos(teta2) to -i*exp(-i*teta1)*cos(teta2) !!!
// Also in formula for tan(phi2) I assumed that the real part of H_{p,q} is used!

// For Shur decomposition see: http://ndickson.wordpress.com/2011/07/13/jacobi-eigenvalue-algorithm-schur-decomposition-and-wikipedia/
// My derivations are:
// V = |  c  -s* |
//     |  s   c* |
// // Case 1:                    // Case 2: (this is also the case when Mii = Mjj)
// c = 1 + sqrt(1 + a*b)         c = sqrt(b)
// -s* = a                       s = conj(sqrt(a))
// a = 2.0*Mij/(Mii-Mjj)         a = Mij/norm
// b = 2.0*Mji/(Mii-Mjj)         b = Mji/norm, norm =sqrt(|Mij|^2+|Mji|^2)

// New parameters: 
// is_cycle:
//       0 - will be using max non-diagonal element
//       1 - cyclic order will be used
//
// alg: - choose algorithm
//       0 - Shur decomposition
//       1 - Jacobi rotations


  int n = n_rows; // = n_cols
  int row, col, i, j, k, num_iter;
  double val,phi,eps;

  CMATRIX V(n,n);
  CMATRIX temp(n,n);

  for(i=0;i<n_elts;i++){     temp.M[i] = M[i];   }
  EVECT.load_identity();
  num_iter = 0;

/*
  // Define convergence criteria.
  k=0; eps = 0.0;
  for(i=0;i<temp.n_rows;i++){
    for(j=0;j<temp.n_cols;j++){
      if(i!=j) {eps+=norm(temp.M[k]); }
       k++;
    }
  }
*/
  
  row = 0; col = 1;

  do{
  //while(eps>EPS){
    num_iter++;

    //cout<<"num_iter = "<<num_iter<<"  eps = "<<eps<<endl;
    if(!is_cycle){  temp.max_nondiagonal(row,col); }   
    
    if(alg==0){
    // Shur rotation   
      complex<double> c,s;
      double L;
      // Case 1
      if((norm(temp.M[row*n_cols+col])<1e+28*norm(temp.M[row*n_cols+row]-temp.M[col*n_cols+col]) ) &&
         (norm(temp.M[row*n_cols+row]-temp.M[col*n_cols+col])>0.0)
        )
      {
        complex<double> a = 2.0*temp.M[row*n_cols+col]/(temp.M[row*n_cols+row]-temp.M[col*n_cols+col]);
        complex<double> b = 2.0*temp.M[col*n_cols+row]/(temp.M[row*n_cols+row]-temp.M[col*n_cols+col]);      
        c = 1.0 + sqrt(1.0 + a*b);
        s = std::conj(-a);
      }
      // Case 2
      else{
        complex<double> a = temp.M[row*n_cols+col];
        complex<double> b = temp.M[col*n_cols+row];
        double nrm = sqrt(norm(a)+norm(b));

        if(nrm>0.0){
          c = sqrt(b/nrm);
          s = std::conj(sqrt(a/nrm));
        }else{
          c = 1.0;
          s = 0.0;
        }

      }
      L = sqrt(norm(c) + norm(s));
      c = c/L;  s = s/L;

      V.load_identity();
      V.M[row*n_cols + row] = c;   V.M[row*n_cols + col] = std::conj(-s);
      V.M[col*n_cols + row] = s;   V.M[col*n_cols + col] = std::conj(c);
   }// Shur decomposition

   else if(alg==1){
     // Jacobi rotation
     double phi1, phi2;
     phi1 = atan2(    temp.M[row*n_cols+col].imag(), temp.M[row*n_cols+col].real());
     phi2 = atan2(2.0*temp.M[row*n_cols+col].real(),(temp.M[row*n_cols+row].real()-temp.M[col*n_cols+col].real()));
     double tet1,tet2;
     tet1 = 0.25*(M_PI - 2.0*phi1);
     tet2 = 0.5*phi2;

     double s1,s2,c1,c2;
     s1 = sin(tet1); s2 = sin(tet2);
     c1 = cos(tet1); c2 = cos(tet2);
 
     V.load_identity();
     V.M[row*n_cols + row] = complex<double>(-s1*s2,-c1*s2);   V.M[row*n_cols + col] = complex<double>(s1*c2,-c1*c2);
     V.M[col*n_cols + row] = complex<double>(-s1*c2,-c1*c2);   V.M[col*n_cols + col] = complex<double>(-s1*s2,c1*s2);

   }// Jacobi rotation



    //temp = V*temp*V.H();
    temp = V*temp*V.H();
 
    EVECT = EVECT*V.H();
    
    k = 0; eps = 0.0;
    for(i=0;i<temp.n_rows;i++){
      for(j=0;j<temp.n_cols;j++){
        if(i!=j) {eps+=norm(temp.M[k]); }
        k++;
      }// for j
    }// for i

    //cout<<"num_iter = "<<num_iter<<" eps = " <<eps<<endl;

    

    if(is_cycle){
      if(row<(n-2) && col<(n-1)) { col++; }
      else if(row<(n-2) && col==(n-1)) { row++; col = row + 1;}
      else if(row==(n-2) && col==(n-1)){ row = 0; col = 1; }
    }


//  }// while eps>EPS
  }while(eps>EPS && num_iter<max_num_iter);

  if(eps>EPS){
    cout<<"Error: In void CMATRIX::eigen(CMATRIX& EVAL, CMATRIX& EVECT,double EPS,int max_num_iter,int is_cycle,int alg)\n";
    cout<<"Number of iterations num_iter = "<<num_iter<<" exceeded maximal number of iterations max_num_iter = "<<max_num_iter<<endl;
    cout<<"Precision achieved eps = "<<eps<<" is lower then requested accuracy EPS = "<<EPS<<endl;
    cout<<"Convergense failed. Exiting...\n";
    exit(0);
  }



  EVAL = temp;
}

void CMATRIX::QR(CMATRIX& w,CMATRIX& R){
/****************************************************************************
 Very helpful resource:
 http://www.math.umn.edu/~olver/aims_/qr.pdf

 QR decomposition is basically the Gram-Schmidt orthogonaliztion procedue:

 M = Q * R, where Q - is a set of orthonormal vectors and R are the weights
*****************************************************************************/
  int row, col, i, j, k;
  double nrm; // norm
  complex<double> dot; // dot product
  int n = n_rows; // = n_cols

  for(i=0;i<n_elts;i++){     w.M[i] = M[i];   }


  for(i=0;i<n;i++){


    if(i>0){
      // w_k = w_k - (w_k,w_i)*w_i  k = i, i+1, ... n
      for(k=i;k<n;k++){

        dot = complex<double>(0.0,0.0);        
        for(j=0;j<n;j++){ dot = dot + (std::conj(w.M[j*n+k])*w.M[j*n+(i-1)]); }// for j

        dot = std::conj(dot); // This is very tricky part!!!  - arises in case of complex matrixes 

        for(j=0;j<n;j++){ w.M[j*n+k] = w.M[j*n+k] - dot*w.M[j*n+(i-1)];}


      }// for k
    }// i > 0

    
    // Simply normalize i-th column-vector
    nrm = 0.0;
    for(j=0;j<n;j++){ nrm += (std::conj(w.M[j*n+i])*w.M[j*n+i]).real(); }// for j
    nrm = sqrt(nrm);
    for(j=0;j<n;j++){ w.M[j*n+i] /= nrm; }


  }// for i


  // Now for R-CMATRIX
  R = complex<double>(0.0,0.0);
  for(i=0;i<n;i++){
    for(j=i;j<n;j++){      
      for(k=0;k<n;k++){
        // R[i][j] = w_j * u_i, note w - is actually original CMATRIX M, while u is what is now w.
        // Note: For complex (this) case the actual definition of the R[i][j] coefficients is:
        // R[j][i] = (w_j^*  x  u_i)^* = w_j * u_i^*, where ^* - denotes complex conjugation
        R.M[i*n+j] += (M[k*n+j]) * std::conj(w.M[k*n+i]);
      }// for k
    }// for j
  }// for i


}


void CMATRIX::QR1(CMATRIX& w,CMATRIX& R){
/****************************************************************************
 Very helpful resource:
 http://www.math.umn.edu/~olver/aims_/qr.pdf

 QR decomposition is basically the Gram-Schmidt orthogonaliztion procedue:

 M = Q * R, where Q - is a set of orthonormal vectors and R are the weights

 This version is designed for tridiagonal matrices
*****************************************************************************/
  int row, col, i, j, k;
  double nrm; // norm
  complex<double> dot; // dot product
  int n = n_rows; // = n_cols

  for(i=0;i<n_elts;i++){     w.M[i] = M[i];   }


  for(i=0;i<n;i++){


    if(i>0){
      // w_k = w_k - (w_k,w_i)*w_i  k = i, i+1
      for(k=i;k<=min((n-1),(i+1));k++){

        dot = complex<double>(0.0,0.0);
        // k = i, i+1 - two or 1 term in dot product computations
        for(j=0;j<=min((i+2),(n-1));j++){ dot += (std::conj(w.M[j*n+k])*w.M[j*n+(i-1)]); }        

        dot = std::conj(dot); // This is very tricky part!!!  - arises in case of complex matrixes 


        // k = i or i+1
        for(j=0;j<=min(i,(n-1));j++) { w.M[j*n+k] = w.M[j*n+k] - dot*w.M[j*n+(i-1)];}


      }// for k
    }// i > 0

    
    // Simply normalize i-th column-vector
    nrm = 0.0;
    for(j=0;j<=min(n-1,(i+1));j++){ nrm += (std::conj(w.M[j*n+i])*w.M[j*n+i]).real(); }// for j
    nrm = sqrt(nrm);
    for(j=0;j<=min(n-1,(i+1));j++){ w.M[j*n+i] /= nrm; }


  }// for i


  // Now for R-CMATRIX
  R = complex<double>(0.0,0.0);
  for(i=0;i<n;i++){
    for(j=i;j<=min(n-1,i+2);j++){      
      for(k=0;k<n;k++){
        // R[i][j] = w_j * u_i, note w - is actually original CMATRIX M, while u is what is now w.
        // Note: For complex (this) case the actual definition of the R[i][j] coefficients is:
        // R[j][i] = (w_j^*  x  u_i)^* = w_j * u_i^*, where ^* - denotes complex conjugation
        //if(j-i==0 || j-i==1 || j-i==2){
          R.M[i*n+j] += (M[k*n+j]) * std::conj(w.M[k*n+i]);
        //}
      }// for k
    }// for j
  }// for i


}


void CMATRIX::bin_dump(std::string filename){

  std::ofstream f(filename.c_str(), ios::out|ios::binary);

  if(f.is_open()){
    f.seekp(0);
    f.write((char*)M, sizeof(complex<double>)*n_elts);
    f.close();    
  }
  else{  cout<<"File "<<filename<<" cann't be open\n"; }
}
 
void CMATRIX::bin_load(std::string filename){

  std::ifstream f(filename.c_str(), ios::in|ios::binary);

  if(f.is_open()){
    f.seekg(0);
    f.read((char *)M, sizeof(complex<double>)*n_elts);
    f.close();   
  }
  else{  cout<<"File "<<filename<<" cann't be open\n"; }


}



void qr(double EPS,int n,CMATRIX& eval,vector<double>& Eval){
// --------- Recursive QR iterations ------------
// eval - is the input tridiagonal CMATRIX
// n - is a size of the problem

  CMATRIX Q(n,n);
  CMATRIX R(n,n);
  int iter = 0;
  int stop = 0;


  do{
       
    
    eval.QR1(Q,R);

    eval = 0.0;

    // The following stepas are basically the efficient way to do:
    // eval = R * Q
    // Fill out the main diagonal
    for(int i=0;i<n;i++){ 
      if(i==n-1){  eval.M[i*n+i]  = R.M[i*n+i]*Q.M[i*n+i];   }          // only 1 term here
      else{        eval.M[i*n+i]  = R.M[i*n+i]*Q.M[i*n+i] + R.M[i*n+(i+1)]*Q.M[(i+1)*n+i];}  // in fact just only 2 terms here
      // Wilkinson shift:
      //eval.M[i*n+i] += mu;
 
    }
    // Fill out upper diagonal
    for(i=0;i<n-1;i++){ 
      // j = i+1
      if(i==n-2){  eval.M[i*n+(i+1)]  = R.M[i*n+i]*Q.M[i*n+(i+1)] + 
                                        R.M[i*n+(i+1)]*Q.M[(i+1)*n+(i+1)];}   // only 2 terms here
      else{        eval.M[i*n+(i+1)]  = R.M[i*n+i]*Q.M[i*n+(i+1)] + 
                                        R.M[i*n+(i+1)]*Q.M[(i+1)*n+(i+1)] +
                                        R.M[i*n+(i+2)]*Q.M[(i+2)*n+(i+1)];}  // all 3 terms

      // The lower diagonal - is by hermitian symmetry:
      eval.M[(i+1)*n+i] = complex<double>(eval.M[i*n+(i+1)].real(), -eval.M[i*n+(i+1)].imag());
    }


    // m has a tridiagonal form, so judge convergence by the elements in
    // the closest off-diagonal
    stop = 0;
    for(i=0;i<(n-1);i++){  
      if(  (fabs(eval.M[i*n+(i+1)].real())<EPS) && (fabs(eval.M[i*n+(i+1)].imag())<EPS) ){

         // Element (i, j=i+1) is  "zero"
         int sz_up = i+1;
         int sz_dn = n-i-1;

         if(sz_up==1 && sz_dn==1){ // Here we just finished 2x2 CMATRIX - done
           Eval[0] = eval.M[0].real();
           Eval[1] = eval.M[3].real();
         }
         else if(sz_up==1 && sz_dn>1){
           CMATRIX dn(sz_dn,sz_dn); dn = 0.0;
           // copy diagonal elements
           for(int j=i+1;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=i+1;j<(n-1);j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=i+2;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))-1] = eval.M[j*n+j-1]; }

           vector<double> Eval_tmp(sz_dn,0.0);
           qr(EPS,sz_dn,dn,Eval_tmp);

           Eval[0] = eval.M[0].real();
           for(j=1;j<n;j++){ Eval[j] = Eval_tmp[j-1]; } 
           Eval_tmp.clear();

         }

         else if(sz_up>1 && sz_dn==1){
           CMATRIX up(sz_up,sz_up); up = 0.0;
           // copy diagonal elements
           for(int j=0;j<(n-1);j++){ up.M[j*sz_up + j] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=0;j<(n-2);j++){ up.M[j*sz_up + j+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=1;j<(n-1);j++){ up.M[j*sz_up + j-1] = eval.M[j*n+j-1]; }

           vector<double> Eval_tmp(sz_up,0.0);
           qr(EPS,sz_up,up,Eval_tmp);


           for(j=0;j<(n-1);j++){ Eval[j] = Eval_tmp[j]; } 
           Eval[n-1] = eval.M[(n-1)*n+(n-1)].real();
           Eval_tmp.clear();

         }

         else{
         // General case - both matrices are at least 2x2
           CMATRIX dn(sz_dn,sz_dn); dn = 0.0;
           // copy diagonal elements
           for(int j=i+1;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=i+1;j<(n-1);j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=i+2;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))-1] = eval.M[j*n+j-1]; }

           CMATRIX up(sz_up,sz_up); up = 0.0;
           // copy diagonal elements
           for(j=0;j<(i+1);j++){ up.M[j*sz_up + j] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=0;j<i;j++){ up.M[j*sz_up + j+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=1;j<(i+1);j++){ up.M[j*sz_up + j-1] = eval.M[j*n+j-1]; }


           vector<double> Eval_tmp_up(sz_up,0.0);
           qr(EPS,sz_up,up,Eval_tmp_up);

           vector<double> Eval_tmp_dn(sz_dn,0.0);
           qr(EPS,sz_dn,dn,Eval_tmp_dn);

           for(j=0;j<sz_up;j++){ Eval[j] = Eval_tmp_up[j]; } 
           for(j=i+1;j<n;j++){ Eval[j] = Eval_tmp_dn[j-(i+1)]; } 

           Eval_tmp_up.clear();
           Eval_tmp_dn.clear();
 
         }



         stop = 1;
      }// if

    }// for i

    iter++;
  }while(!stop);


}


void qr(double EPS,int n,CMATRIX& eval,vector<double>& Eval,CMATRIX& Evec){
// Overloading qr function - to keep track of the transformation
// CMATRIX
// --------- Recursive QR iterations ------------
// eval - is the input tridiagonal CMATRIX
// n - is a size of the problem
  double mu,d,a1,b2;
  int i1,i2;
  CMATRIX Q(n,n);
  CMATRIX R(n,n);
  CMATRIX Q_tmp(n,n);  Q_tmp.load_identity();
  CMATRIX I(n,n); I.load_identity();
  int iter = 0;
  int stop = 0;

  Evec.load_identity();

  do{
       
    // Perhaps they mean - minimal diagonal value
    a1 = eval.M[0].real(); i1 = 0;
    for(int i=1;i<n;i++){  d = eval.M[i*n+i].real(); if(d<a1){ a1 = d; i1 = i;} }

/*
    // Wilkinson shifts
    b    = eval.M[(n-2)*n+(n-1)];

    delta = 0.5*(an_1 - an);
    if(delta<0){ sgn = -1.0; }
    else if(delta>0){ sgn = 1.0; }
    else if(delta==0.0){ sgn = 0.0; }

    b2    = norm(b);
    mu    = an - sgn*b2/(fabs(delta) + sqrt(delta*delta + b2));

*/
    if(n==2){ mu = 0.0; }
    else{ mu = a1; }


    (eval-mu*I).QR1(Q,R);

    Evec *= Q;

    eval = 0.0;

    // The following steps are basically the efficient way to do:
    // eval = R * Q
    // Fill out the main diagonal
    for(i=0;i<n;i++){ 
      if(i==n-1){  eval.M[i*n+i]  = R.M[i*n+i]*Q.M[i*n+i];   }          // only 1 term here
      else{        eval.M[i*n+i]  = R.M[i*n+i]*Q.M[i*n+i] + R.M[i*n+(i+1)]*Q.M[(i+1)*n+i];}  // in fact just only 2 terms here
      // Shift:
      eval.M[i*n+i] += mu;
 
    }
    // Fill out upper diagonal
    for(i=0;i<n-1;i++){ 
      // j = i+1
      if(i==n-2){  eval.M[i*n+(i+1)]  = R.M[i*n+i]*Q.M[i*n+(i+1)] + 
                                        R.M[i*n+(i+1)]*Q.M[(i+1)*n+(i+1)];}   // only 2 terms here
      else{        eval.M[i*n+(i+1)]  = R.M[i*n+i]*Q.M[i*n+(i+1)] + 
                                        R.M[i*n+(i+1)]*Q.M[(i+1)*n+(i+1)] +
                                        R.M[i*n+(i+2)]*Q.M[(i+2)*n+(i+1)];}  // all 3 terms

      // The lower diagonal - is by hermitian symmetry:
      eval.M[(i+1)*n+i] = complex<double>(eval.M[i*n+(i+1)].real(), -eval.M[i*n+(i+1)].imag());
    }


    // m has a tridiagonal form, so judge convergence by the elements in
    // the closest off-diagonal
    stop = 0;
    for(i=0;i<(n-1);i++){  
      if(  (fabs(eval.M[i*n+(i+1)].real())<EPS) && (fabs(eval.M[i*n+(i+1)].imag())<EPS) ){

         // Element (i, j=i+1) is  "zero"
         int sz_up = i+1;
         int sz_dn = n-i-1;

         if(sz_up==1 && sz_dn==1){ // Here we just finished 2x2 CMATRIX - done
           Eval[0] = eval.M[0].real();
           Eval[1] = eval.M[3].real();
         }
         else if(sz_up==1 && sz_dn>1){
           CMATRIX dn(sz_dn,sz_dn); dn = 0.0;
           // copy diagonal elements
           for(int j=i+1;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=i+1;j<(n-1);j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=i+2;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))-1] = eval.M[j*n+j-1]; }

           vector<double> Eval_tmp(sz_dn,0.0);
           CMATRIX Q_dn(sz_dn,sz_dn);
           qr(EPS,sz_dn,dn,Eval_tmp,Q_dn);

           for(j=i+1;j<n;j++){
             for(int k=i+1;k<n;k++){
               Q_tmp.M[j*n+k] = Q_dn.M[(j-(i+1))*sz_dn + (k-(i+1))];
             }
           }

           Eval[0] = eval.M[0].real();
           for(j=1;j<n;j++){ Eval[j] = Eval_tmp[j-1]; } 
           Eval_tmp.clear();

         }

         else if(sz_up>1 && sz_dn==1){
           CMATRIX up(sz_up,sz_up); up = 0.0;
           // copy diagonal elements
           for(int j=0;j<(n-1);j++){ up.M[j*sz_up + j] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=0;j<(n-2);j++){ up.M[j*sz_up + j+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=1;j<(n-1);j++){ up.M[j*sz_up + j-1] = eval.M[j*n+j-1]; }

           vector<double> Eval_tmp(sz_up,0.0);
           CMATRIX Q_up(sz_up,sz_up);
           qr(EPS,sz_up,up,Eval_tmp,Q_up);

           for(j=0;j<(i+1);j++){
             for(int k=0;k<(i+1);k++){
               Q_tmp.M[j*n+k] = Q_up.M[j*sz_up+k];
             }
           }

           for(j=0;j<(n-1);j++){ Eval[j] = Eval_tmp[j]; } 
           Eval[n-1] = eval.M[(n-1)*n+(n-1)].real();
           Eval_tmp.clear();

         }

         else{
         // General case - both matrices are at least 2x2
           CMATRIX dn(sz_dn,sz_dn); dn = 0.0;
           // copy diagonal elements
           for(int j=i+1;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=i+1;j<(n-1);j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=i+2;j<n;j++){ dn.M[(j-(i+1))*sz_dn + (j-(i+1))-1] = eval.M[j*n+j-1]; }

           CMATRIX up(sz_up,sz_up); up = 0.0;
           // copy diagonal elements
           for(j=0;j<(i+1);j++){ up.M[j*sz_up + j] = eval.M[j*n+j]; }
           // upper off-diagonal elements
           for(j=0;j<i;j++){ up.M[j*sz_up + j+1] = eval.M[j*n+j+1]; }
           // lower off-diagonal elements
           for(j=1;j<(i+1);j++){ up.M[j*sz_up + j-1] = eval.M[j*n+j-1]; }

           vector<double> Eval_tmp_up(sz_up,0.0);
           CMATRIX Q_up(sz_up,sz_up);
           qr(EPS,sz_up,up,Eval_tmp_up,Q_up);

           vector<double> Eval_tmp_dn(sz_dn,0.0);
           CMATRIX Q_dn(sz_dn,sz_dn);
           qr(EPS,sz_dn,dn,Eval_tmp_dn,Q_dn);


           for(j=0;j<(i+1);j++){
             for(int k=0;k<(i+1);k++){
               Q_tmp.M[j*n+k] = Q_up.M[j*sz_up+k];
             }
           }

           for(j=i+1;j<n;j++){
             for(int k=i+1;k<n;k++){
               Q_tmp.M[j*n+k] = Q_dn.M[(j-(i+1))*sz_dn + (k-(i+1))];
             }
           }


           for(j=0;j<sz_up;j++){ Eval[j] = Eval_tmp_up[j]; } 
           for(j=i+1;j<n;j++){ Eval[j] = Eval_tmp_dn[j-(i+1)]; } 

           Eval_tmp_up.clear();
           Eval_tmp_dn.clear();
 
         }

         Evec *= Q_tmp;

         stop = 1;
      }// if

    }// for i

    iter++;


  }while(!stop);

}

void CMATRIX::eigen(double EPS,CMATRIX& EVAL,CMATRIX& EVECT,int opt){
// This is just a convenient interface
// opt - is option - choose the method

  vector<double> Eval(n_rows,0.0);
  EVAL = 0.0;

  if(n_rows==2){ opt = 1; } // some tweak

  if(opt==1){   eigen0(EVAL,EVECT,EPS,10000,0,0);  }  // this works up to ~ n =50 
  else if(opt==2){  eigen2(EPS,Eval,EVECT);
    for(int i=0;i<n_rows;i++){ EVAL.M[i*n_cols+i] = Eval[i]; } // this is fastest version, able to work up to ~n = 250
  }
  else if(opt==3){  eigen3(EPS,Eval,EVECT);
    for(int i=0;i<n_rows;i++){ EVAL.M[i*n_cols+i] = Eval[i]; } // this works up to ~n = 200, but is slow
  }
  

}

void CMATRIX::eigen1(double EPS,vector<double>& Eval){
//-------------------------------------------------------------
// We do the job in reductionist way - once one of the elements on 
// the off-diagonal is smaller than EPS - we split the CMATRIX into
// 2 blocks and then deal with each other independently - deflation
//-------------------------------------------------------------

  int n = n_rows; // = n_cols

  CMATRIX Q(n,n);
  CMATRIX R(n,n);
  CMATRIX eval(n,n);


  tridiagonalize(eval);

  // To avoid error propagation we just set the off-tridiagonal elements to 0.0
  int k = 0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(abs(i-j)>1){ eval.M[k] = 0.0; } k++;
    }
  }

  qr(EPS,n,eval,Eval);

}


void CMATRIX::eigen2(double EPS,vector<double>& Eval,CMATRIX& Evec){
//-------------------------------------------------------------
// This is practically the same version as eigen1, but we also 
// keep track of the transformation matrixes - so to compute all
// eigenvectors
// The relation is:
// this * Evec = Evec * Eval  or (because Evec.H() * Evec = I)
// Eval = Evec.H() * this * Evec
//-------------------------------------------------------------

  int n = n_rows; // = n_cols

  Evec.load_identity();
  CMATRIX Q(n,n);
  CMATRIX T(n,n);

  // M =  H * T * H,  H^H = H^-1 = H, here H = Evec - to save space
  tridiagonalize(T,Evec);

  // To avoid error propagation we just set the off-tridiagonal elements to 0.0
  int k = 0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(abs(i-j)>1){ T.M[k] = 0.0; } k++;
    }
  }

  qr(EPS,n,T,Eval,Q);//    T = Q.H()*Eval*Q
  Evec *= Q;

}

void CMATRIX::eigen3(double EPS,vector<double>& Eval,CMATRIX& Evec){
//-------------------------------------------------------------
// This is practically the same version as eigen1, but we also 
// keep track of the transformation matrixes - so to compute all
// eigenvectors
// The relation is:
// this * Evec = Evec * Eval  or (because Evec.H() * Evec = I)
// Eval = Evec.H() * this * Evec
//-------------------------------------------------------------

  int n = n_rows; // = n_cols

  CMATRIX X(n,1);
  complex<double> gs( (1.0/(sqrt(2.0)*n)), (1.0/(sqrt(2.0)*n)) );

  CMATRIX m(n,n); m = *this;

  m.eigen1(EPS,Eval); // compute all eigenvalues


  for(int i=0;i<n;i++){

    // Instead of : A = m-Eval[i]*I;
    for(int j=0;j<n;j++){ m.M[j*n+j] -= Eval[i]; }
    // Initial guess
    for(j=0;j<n;j++){ X.M[j] = gs; }

    solve_linsys1(m,X,EPS,10000,1.4); // ~1.4 is optimum

    // Restore original m
    for(j=0;j<n;j++){ m.M[j*n+j] += Eval[i]; }

    // Compute norm of the solution vector
    double nrm = 0.0;
    for(j=0;j<n;j++){ nrm += norm(X.M[j]); }
    nrm = sqrt(1.0/nrm);
 
    // Normalize solution vector
    for(j=0;j<n;j++){ Evec.M[j*n+i] = nrm*X.M[j]; }

  }//for i

}



void CMATRIX::tridiagonalize(CMATRIX& T){
/****************************************************************************
  Here we will transfrom the CMATRIX to tridiagonal form (T) using the Householder
  transformations. We assume our CMATRIX is Hermitian or symmetric, otherwise
  the algorithm may not work

  General identity behind Householder:
  Y = P * X
  then P = I - 2*W*W^T
  W = (Y - X)/|Y-X|
  http://math.fullerton.edu/mathews/n2003/HouseholderMod.html

  This algorithm scales as O(N^3)
*****************************************************************************/
  int i,j,k,n;
  double nrm;
  complex<double> alp;

  n = n_rows; // = n_cols


  CMATRIX w(n,1); 
  CMATRIX v(n,1);

  for(i=0;i<n_elts;i++){ T.M[i] = M[i]; }

  for(i=0;i<(n-2);i++){

    w = 0.0;
    nrm = 0.0; for(j=i+1;j<n;j++){ nrm += (std::conj(T.M[j*n+i])*T.M[j*n+i]).real(); }  nrm = sqrt(nrm);

    if(abs(T.M[(i+1)*n+i])==0.0){ alp = complex<double>(1.0,0.0); }
    else{ alp = (T.M[(i+1)*n+i] / abs(T.M[(i+1)*n+i])) ; }
    
    w.M[i+1] = T.M[(i+1)*n+i] - alp*nrm; 

    nrm = (std::conj(w.M[i+1])*w.M[i+1]).real(); // norm of new vector x-y
    for(j=i+2;j<n;j++){ w.M[j] = T.M[j*n+i]; nrm += (std::conj(w.M[j])*w.M[j]).real(); }
    nrm = sqrt(nrm);

    // Normalize new vector (w)
    for(j=i+1;j<n;j++){ w.M[j] = w.M[j] / nrm; }
    
    // The following commented lines are only for mathematical and historical reason
    // Finally, projector
    // P = iden - 2.0 * w * (w.H());      
    // Do the transformation of the CMATRIX
    // This is basic algorithm
    //T = P * T * P; // P = P.H = P^-1   // <-- this is most costly (both memory and time) place!

    // Let's optimize it:
    v = T * w;    
    T =  (T - 2.0*w*(v.H()) - 2.0*v*(w.H()) + 4.0*(w.H()*v).M[0]*w*(w.H()) );
    

  }// for i - transformation index

}


void CMATRIX::tridiagonalize(CMATRIX& T,CMATRIX& H){
/****************************************************************************
  Here we will transfrom the CMATRIX to tridiagonal form (T) using the Householder
  transformations. We assume our CMATRIX is Hermitian or symmetric, otherwise
  the algorithm may not work

  The full transformation is kept in CMATRIX H

  General identity behind Householder:
  Y = P * X
  then P = I - 2*W*W^T
  W = (Y - X)/|Y-X|
  http://math.fullerton.edu/mathews/n2003/HouseholderMod.html

  This algorithm scales as O(N^3)
*****************************************************************************/
  int i,j,k,n;
  double nrm;
  complex<double> alp;

  n = n_rows; // = n_cols

  CMATRIX I(n,n);  I.load_identity();
  CMATRIX tmp1(n,n);
  CMATRIX tmp2(n,n);
  CMATRIX P(n,n);
  H.load_identity();
  CMATRIX w(n,1); 
  CMATRIX v(n,1);

  for(i=0;i<n_elts;i++){ T.M[i] = M[i]; }

  for(i=0;i<(n-2);i++){

    w = 0.0;
    nrm = 0.0; for(j=i+1;j<n;j++){ nrm += (std::conj(T.M[j*n+i])*T.M[j*n+i]).real(); }  nrm = sqrt(nrm);

    if(abs(T.M[(i+1)*n+i])==0.0){ alp = complex<double>(1.0,0.0); }
    else{ alp = (T.M[(i+1)*n+i] / abs(T.M[(i+1)*n+i])) ; }
    
    w.M[i+1] = T.M[(i+1)*n+i] - alp*nrm; 

    nrm = (std::conj(w.M[i+1])*w.M[i+1]).real(); // norm of new vector x-y
    for(j=i+2;j<n;j++){ w.M[j] = T.M[j*n+i]; nrm += (std::conj(w.M[j])*w.M[j]).real(); }
    nrm = sqrt(nrm);

    // Normalize new vector (w)
    for(j=i+1;j<n;j++){ w.M[j] = w.M[j] / nrm; }
    
    // The following commented lines are only for mathematical and historical reason
    // Finally, projector
    // P = iden - 2.0 * w * (w.H());      
    // Do the transformation of the CMATRIX
    // This is basic algorithm
    //T = P * T * P; // P = P.H = P^-1   // <-- this is most costly (both memory and time) place!

    // Let's optimize it:
    v = T * w;    
    P = w*(w.H());

//  The following 2 lines are what we actually doing:
    T =  (T - 2.0*w*(v.H()) - 2.0*v*(w.H()) + 4.0*(w.H()*v).M[0]*P );
    P *= -2.0;
    P += I;
    H = H * P;  

  }// for i - transformation index

}



CMATRIX exp(CMATRIX& m1,complex<double> scl,double eps){
/****************************************************************************
  Computes  exp(m1*scl)
  Works only for Hermitian m1: m1.H() = m1
*****************************************************************************/
  if(m1.n_rows!=m1.n_cols){ cout<<"Error in exp: Can not exponentiate non-square CMATRIX\n"; exit(0); }
  int n = m1.n_rows;
  CMATRIX evec(n,n),eval(n,n);//,inv_evec(n,n);

  m1.eigen(eps,eval,evec,2);
  for(int i=0;i<n;i++){  eval.M[i*n+i] = exp(eval.M[i*n+i].real()*scl); }
  //evec.direct_inverse(eps,inv_evec);  inv_evec = evec.H()
  return (evec*eval*evec.H());
}
 

CMATRIX sin(CMATRIX& m1,complex<double> scl, double eps){
/****************************************************************************
  Computes sin(m1*scl)
  Works only for Hermitian m1: m1.H() = m1
*****************************************************************************/
  if(m1.n_rows!=m1.n_cols){ cout<<"Error in exp: Can not exponentiate non-square CMATRIX\n"; exit(0); }
  int n = m1.n_rows;
  CMATRIX evec(n,n),eval(n,n);//,inv_evec(n,n);

  m1.eigen(eps,eval,evec,2);
  for(int i=0;i<n;i++){  eval.M[i*n+i] = sin(eval.M[i*n+i].real()*scl); }
//  evec.direct_inverse(eps,inv_evec);
  return (evec*eval*evec.H());
}

CMATRIX cos(CMATRIX& m1,complex<double> scl, double eps){
/****************************************************************************
  Computes cos(m1*scl)
  Works only for Hermitian m1: m1.H() = m1
*****************************************************************************/
  if(m1.n_rows!=m1.n_cols){ cout<<"Error in exp: Can not exponentiate non-square CMATRIX\n"; exit(0); }
  int n = m1.n_rows;
  CMATRIX evec(n,n),eval(n,n);//,inv_evec(n,n);

  m1.eigen(eps,eval,evec,2);
  for(int i=0;i<n;i++){  eval.M[i*n+i] = cos(eval.M[i*n+i].real()*scl); }
//  evec.direct_inverse(eps,inv_evec);
  return (evec*eval*evec.H());
}

CMATRIX pow(CMATRIX& m1,double nn, double eps){
/****************************************************************************
  Computes pow(m1,nn). In particular if nn = 1/2 result is sqrt(m1)
  Works only for Hermitian m1: m1.H() = m1
*****************************************************************************/
  if(m1.n_rows!=m1.n_cols){ cout<<"Error in exp: Can not exponentiate non-square CMATRIX\n"; exit(0); }
  int n = m1.n_rows;
  CMATRIX evec(n,n),eval(n,n);//,inv_evec(n,n);

  m1.eigen(eps,eval,evec,2);
  for(int i=0;i<n;i++){  eval.M[i*n+i] = std::pow(eval.M[i*n+i].real(),nn); }
//  evec.direct_inverse(eps,inv_evec);
  return (evec*eval*evec.H());
}




void CMATRIX::inverse(double EPS,CMATRIX& INV,int opt){

  if(opt==1){ direct_inverse(EPS,INV); }  // this is much faster way - works fine for ~n = 350 and more
  else if(opt==2){                        // actually this is slower version - works only up ~n = 100
    CMATRIX I(n_rows,n_cols); I.load_identity();
    solve_linsys(*this, I, INV, EPS, 10000, 1.4); // CX = D, so if D = I => X = C^-1  
  }

}



void CMATRIX::direct_inverse(double EPS,CMATRIX& INV){
  int num_of_rows = n_rows;
  int num_of_cols = n_cols;

  double zero = EPS;
  complex<double> *R_time;   R_time=new complex<double>[num_of_rows*num_of_cols];
  complex<double> *L_time;   L_time=new complex<double>[num_of_rows*num_of_cols];

  for(int k=0;k<num_of_rows*num_of_cols;k++){ R_time[k]=M[k]; }

  k=0;
  for(int i=0;i<num_of_cols;i++){
    for(int j=0;j<num_of_cols;j++){
      if(i==j) {L_time[k]=complex<double>(1.0,0.0);}
      else     {L_time[k]=complex<double>(0.0,0.0);}
      k++;
    }// for j
  }// for i


  complex<double> alpha;
  for(int row1=0;row1<num_of_rows-1;row1++){

    // Diagonal element is zero
    if(abs(R_time[row1*num_of_cols+row1])<=zero){

      // Find the row with the element in row1 column being not zero
      int row=row1+1;
      while(abs(R_time[row*num_of_cols+row1])<=zero){ row++;}
      complex<double> temp1,temp2;
      // Swap rows  <row> and <row1> in R_temp and L_temp matrices
      for(int col=0;col<num_of_cols;col++){
        temp1=R_time[row1*num_of_cols+col];
        R_time[row1*num_of_cols+col]=R_time[row*num_of_cols+col];
        R_time[row*num_of_cols+col]=temp1;

        temp2=L_time[row1*num_of_cols+col];
        L_time[row1*num_of_cols+col]=L_time[row*num_of_cols+col];
        L_time[row*num_of_cols+col]=temp2;
      }// for col
    }// if

    // The diagonal element is non-zero (originally or after swapping just above)
    if(abs(R_time[row1*num_of_cols+row1])>zero){

      // Eliminate all the elements unde this diagonal element
      for(int row2=row1+1;row2<num_of_rows;row2++){
        //if(abs(R_time[row2*num_of_cols+row1])>zero){

          alpha=-R_time[row2*num_of_cols+row1]/R_time[row1*num_of_cols+row1];

          for(int col=0;col<num_of_cols;col++){
            R_time[row2*num_of_cols+col]=R_time[row2*num_of_cols+col]+alpha*R_time[row1*num_of_cols+col];
            L_time[row2*num_of_cols+col]=L_time[row2*num_of_cols+col]+alpha*L_time[row1*num_of_cols+col];
          }// for col

        //}// if !=0
        //else continue;
      }// for row2
    }// if !=0
    else{ cout<<"Error in direct_inverse: The leading element is smaller than "<<zero<<endl; exit(0); }
  }// for row1

  for(row1=num_of_rows-1;row1>0;row1--){
    alpha=R_time[row1*num_of_cols+row1];
    R_time[row1*num_of_cols+row1]=complex<double>(1.0,0.0);

    for(int col=(num_of_cols-1);col>=0;col--){ 
      L_time[row1*num_of_cols+col]=L_time[row1*num_of_cols+col]/alpha;
    }// for col
    for(int row2=row1-1;row2>=0;row2--){
      alpha=-R_time[row2*num_of_cols+row1];
      for(int col=(num_of_cols-1);col>=0;col--){
        R_time[row2*num_of_cols+col]=R_time[row2*num_of_cols+col]+alpha*R_time[row1*num_of_cols+col];
        L_time[row2*num_of_cols+col]=L_time[row2*num_of_cols+col]+alpha*L_time[row1*num_of_cols+col];
      }// for col
    }// for row2
  }// for row1

  alpha=R_time[0];
  R_time[0]=complex<double>(1.0,0.0);
  for(int col=(num_of_cols-1);col>=0;col--){   L_time[col]=L_time[col]/alpha;  }

  k=0;
  for(int row=0;row<num_of_rows;row++){
    for(int col=0;col<num_of_cols;col++){
      INV.M[row*num_of_cols+col] = L_time[k];
      k++;
    }
  }
  delete [] R_time;
  delete [] L_time;
}


void solve_linsys(CMATRIX& C,CMATRIX& D, CMATRIX& X,double eps,int maxiter,double omega){
/*********************************************
 Here we solve the system of linear equations
      CX = D  --->     AX = D', where  A = C.T()*C   D' = C.T()*D
 using Gauss-Seidel iterative procedure

 Inputs: A, D - matrices
         eps  - precision criterion
         omega - parameter of overrelaxation - must be in range (1,2)
                 to accelerate convergence
 Output: X

 Some preliminary transformations are made in order
 to be able to use Gauss-Seidel method for any CMATRIX A

 More details:
 80.47 An iterative Algorithm for Matrix Inversion
 Which is Always Convergent
 Authors: S. Simons
 Source: The Mathematical Gazette, Vol. 80, No. 489
 (Nov., 1996), pp. 567-569

 url: http://www.jstor.org/stable/pdfplus/3618529.pdf

**********************************************/

// Do the transformations A = C^H * C and b = C^H * d
// If matrices d and c have more then 1 columns we do the
// procedure for each column


    int i,j,k;   // counters
    int n,m,p;   // dimetions
    complex<double> s;    // sums
    double error;// error
    int iter;    // number of iterations

    if(C.n_rows!=D.n_rows)
        {std::cout<<"Error: The number of rows of matrices C and D in equation CX = D must be equal\n"; exit(35); } // n
    if(C.n_cols!=X.n_rows)
        {std::cout<<"Error: The number of cols of CMATRIX C and num of rows in CMATRIX D in equation CX = D must be equal\n"; exit(35); } // m
    if(X.n_cols!=D.n_cols)
        {std::cout<<"Error: The number of cols of matrices X and D in equation CX = D must be equal\n"; exit(35); } // p

    // Set dimentions
    n = C.n_rows;
    m = C.n_cols;
    p = D.n_cols;  // this is just 1 in most of the cases

    cout<<"n = "<<n<<endl;
    cout<<"m = "<<m<<endl;
    cout<<"p = "<<p<<endl;

    CMATRIX A(m,m); A = C.H() * C;       
    eps = eps*eps;
    error = 2.0*eps;
    iter = 0;

    CMATRIX d(n,1);
    CMATRIX x(m,1);
    CMATRIX xprev(m,1);
    CMATRIX b(m,1);


    while((error>eps)&&(iter<maxiter)){

    error = 0.0;

    for( k = 0; k < p; k++ ){

        //------- Matrix preparation step -----------

        for(i = 0;i<n;i++){ d.M[i] = D.M[i*p+k]; }
        for(i = 0;i<m;i++){ x.M[i] = X.M[i*p+k]; }
        xprev = x;
        b = C.H() * d;

        //------- Gauss-Seidel step -----------------

        for( i = 0; i < m; i++ ){

            s = 0.0;

            for( j = 0; j < i; j++ ){

                s += A.M[i*m + j]*x.M[j];

            }// for j

            for( j = i+1; j < m; j++ ){

                s += A.M[i*m + j]*xprev.M[j];

            }// for j

            x.M[i] = (b.M[i] - s)/A.M[i*m + i];

        }// for i - all elements of vector x


        //-------- Now calculate the error and update X ---------

        for( i = 0; i < m; i++ ){

            int indx = i*p+k;

            error += ((std::conj(x.M[i] - xprev.M[i]))*(x.M[i] - xprev.M[i])).real();

            X.M[indx] = omega*x.M[i] + (1.0-omega)*X.M[indx]; 

        }// for i



    }// for k - all columns of D

    error = (error/double(m));

    iter++;

    }// loop over convergence


}

void solve_linsys1(CMATRIX& C, CMATRIX& X,double eps,int maxiter,double omega){
/*********************************************
 This is gonna be optimized version of the solve_linsys function
 for the case when D = 0

 omega - is a convergence parameter:
 x^(n+1) = omega * z^(n+1) + (1-omega)*x^n, where:
 
 z^(n+1) - is a normal Gauss-Seidel iterate (that is omega = 1)

 Here we solve the system of linear equations
      CX = D  --->     AX = D', where  A = C.T()*C   D' = C.T()*D
 using Gauss-Seidel iterative procedure

 Inputs: A, D - matrices
         eps  - precision criterion
 Output: X

 Some preliminary transformations are made in order
 to be able to use Gauss-Seidel method for any CMATRIX A

 More details:
 80.47 An iterative Algorithm for Matrix Inversion
 Which is Always Convergent
 Authors: S. Simons
 Source: The Mathematical Gazette, Vol. 80, No. 489
 (Nov., 1996), pp. 567-569

 url: http://www.jstor.org/stable/pdfplus/3618529.pdf

**********************************************/

// Do the transformations A = C^H * C and b = C^H * d
// If matrices d and c have more then 1 columns we do the
// procedure for each column


    int i,j,k;   // counters
    int n,m,p;   // dimetions
    complex<double> s;    // sums
    complex<double> diff;
    double error;// error
    int iter;    // number of iterations

    if(C.n_cols!=X.n_rows)
        {std::cout<<"Error: The number of cols of CMATRIX C and num of rows in CMATRIX D in equation CX = D must be equal\n"; exit(35); } // m
    if(X.n_cols!=1)
        {std::cout<<"Error: The number of cols of matrices X and D in equation CX = D must be equal\n"; exit(35); } // p

    // Set dimentions
    n = C.n_rows;
    m = C.n_cols;
    p = 1;  // this is just 1 in most of the cases

    CMATRIX d(n,1);
    CMATRIX x(m,1);
    CMATRIX xprev(m,1);
   
   
    CMATRIX A(m,m); A = C.H() * C;       
    eps = eps*eps;
    error = 2.0*eps;
    iter = 0;

    int im; // i*m
  

    while((error>eps)&&(iter<maxiter)){

    error = 0.0;

    for( k = 0; k < p; k++ ){

        //------- Matrix preparation step -----------

        for(i = 0;i<m;i++){ xprev.M[i] = x.M[i] = X.M[i*p+k]; }

        //------- Gauss-Seidel step -----------------
        im = 0;
        for( i = 0; i < m; i++ ){

            s = 0.0;

            for( j = 0; j < i; j++ ){

                s += A.M[im + j]*x.M[j];

            }// for j

            for( j = i+1; j < m; j++ ){

                s += A.M[im + j]*xprev.M[j];

            }// for j


            x.M[i] = (-s)/A.M[im + i];


            im += m;

        }// for i - all elements of vector x

        //-------- Now calculate the error and update X ---------

        for( i = 0; i < m; i++ ){

            error += ((std::conj(x.M[i] - xprev.M[i]))*(x.M[i] - xprev.M[i])).real();

            X.M[i] = omega*x.M[i] + (1.0-omega)*X.M[i];  // k takes value of only 0, p = 1, so i*p+k = i

        }// for i



    }// for k - all columns of D

    error = error/double(m);

    iter++;

    }// loop over convergence


}


void dft(CMATRIX& in,CMATRIX& out){
/***************************************
  Discrete Fourier Transform
  e.g. http://en.wikipedia.org/wiki/Fast_Fourier_transform
****************************************/

  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,mul;
  double argg;

  for(int k=0;k<N;k++){

    out.M[k] = 0.0;
    argg = 2.0*M_PI*k/((double)N);

    f = complex<double>(std::cos(argg),-std::sin(argg));
    mul = 1.0;

    for(int n=0;n<N;n++){
      out.M[k] += in.M[n]*mul;
      mul *= f;
    }// for j
  }// for i

}

void inv_dft(CMATRIX& in,CMATRIX& out){
/***************************************
  Inverse Discrete Fourier Transform
  e.g. http://en.wikipedia.org/wiki/Discrete_Fourier_transform
****************************************/

  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,mul;
  double argg;

  for(int k=0;k<N;k++){

    out.M[k] = 0.0;
    argg = 2.0*M_PI*k/((double)N);

    f = complex<double>(std::cos(argg),std::sin(argg));
    mul = 1.0;

    for(int n=0;n<N;n++){
      out.M[k] += in.M[n]*mul;
      mul *= f;
    }// for j
  }// for i

  argg = 1.0/((double)N);
  
  for(k=0;k<N;k++){ out.M[k] *= argg; }

}


void cft(CMATRIX& in,CMATRIX& out,double xmin,double dx){
/***************************************
  Continuous Fourier Transform
  f(k) = Integral ( f(r) * exp(-2*pi*i*k*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*k*r_n)) * dx =
     n
  = dx *exp(-2*pi*i*k*xmin) * sum ( in[n] * exp(-2*pi*k*dx*n)  )
                               n
  r_n = xmin + dx * n
****************************************/

  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,mul;
  double argg;
  double L = dx*N;

  for(int k=0;k<N;k++){

    double K = (k/L);
    complex<double> pref(0.0,-2.0*M_PI*K*xmin);
    mul = dx*std::exp(pref);
    out.M[k] = 0.0;
    argg = 2.0*M_PI*K*dx;
    f = complex<double>(std::cos(argg),-std::sin(argg));

    for(int n=0;n<N;n++){
      out.M[k] += in.M[n]*mul;
      mul *= f;
    }// for j
  }// for i

}

void cft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx){
/***************************************
  Continuous Fourier Transform
  f(kmin + k) = Integral ( f(r) * exp(-2*pi*i*(kmin+k)*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*(kmin+k)*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*(kmin+k)*xmin) * sum ( in[n] * exp(-2*pi*(kmin+k)*dx*n)  )
                                       n
  r_n = xmin + dx * n

****************************************/

  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,mul;
  double argg;
  double L = dx*N;

  for(int k=0;k<N;k++){

    double K = kmin + (k/L);
    complex<double> pref(0.0,-2.0*M_PI*K*xmin);
    mul = dx*std::exp(pref);
    out.M[k] = 0.0;
    argg = 2.0*M_PI*K*dx;
    f = complex<double>(std::cos(argg),-std::sin(argg));

    for(int n=0;n<N;n++){
      out.M[k] += in.M[n]*mul;
      mul *= f;
    }// for j

  }// for i

}

void cfft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx){
/***************************************
  Continuous Fast Fourier Transform

  f(kmin + k) = Integral ( f(r) * exp(-2*pi*i*(kmin+k)*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*(kmin+k)*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*(kmin+k)*xmin) * sum ( in[n] * exp(-2*pi*(kmin+k)*dx*n)  )
                                       n
  r_n = xmin + dx * n

  The size of the grid should be the power of 2
****************************************/

  // Initial check 
  int N = in.n_elts; // N must be 2^n, n - some integer
  if( (N>1) && (N%2!=0) ){  cout<<"N must be power of 2\n"; exit(0);}

  if(N>2){

    // Some constants and variables
    double dk = 1.0/(dx*N);
    int Nhalf = N/2;
    int kp;
    complex<double> one(0.0,1.0);
    complex<double> p1, p2,ea;

    ea = std::exp(-M_PI*one*xmin/dx); // alp(2*dx)
    p1 = std::exp(-2.0*M_PI*one*kmin*dx);
    p2 = std::exp(-2.0*M_PI*one*dk*dx);

    CMATRIX in_even(1,Nhalf);
    CMATRIX in_odd(1,Nhalf);
    CMATRIX out_even(1,Nhalf);
    CMATRIX out_odd(1,Nhalf);

    // Divide set on the even and odd part
    for(int m=0;m<Nhalf;m++){
      in_even.M[m] = in.M[2*m];
      in_odd.M[m] = in.M[2*m+1];
    }

    // Perform fft on smaller grids
    cfft1(in_even,out_even,xmin,kmin,2.0*dx);
    cfft1(in_odd,out_odd,xmin,kmin,2.0*dx);


    // Compute ft of the original grid using the fts for smaller grids
    for(int k=0;k<Nhalf;k++){    
      out.M[k] = 0.5*(out_even.M[k] + p1 * out_odd.M[k]);
      out.M[k+Nhalf] = 0.5*ea*(out_even.M[k] - p1 * out_odd.M[k]);
      p1 *= p2;
    }  


  }// if N>=2 

  else if(N==2){ // No more recursive levels below this one - do all explicitly 
    cft1(in,out,xmin,kmin,dx); 
  }

}


void cft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy){
/***************************************
  Continuous 2-D Fourier Transform
  f(kmin + k) = Integral ( f(r) * exp(-2*pi*i*(kmin+k)*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*(kmin+k)*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*(kmin+k)*xmin) * sum ( in[n] * exp(-2*pi*(kmin+k)*dx*n)  )
                                       n
  r_n = xmin + dx * n

  in and out - are now 2D matrices (columns and rows) - this is for convenience, although they both may be 
               reformulated as 1D matrices (colums or rows)

****************************************/

  // in and out are Nx x Ny matrices
  int Nx = in.n_rows;
  int Ny = in.n_cols; 

  complex<double> sumx,sumy,fx,fy,Px,Py,fnx,fny;
  complex<double> one(0.0,1.0);

  // Real space box [xmin, xmin+Lx] x [ymin, ymin+Ly]
  double argg;
  double Lx = dx*Nx;
  double Ly = dy*Ny;
  double dkx = 1.0/Lx;
  double dky = 1.0/Ly;
  double Kx,Ky;

  // Conceptually:  out.M[kx*Ny+ky] +=  in.M[nx*Ny+ny] * exp(-2.0*M_PI*one*(Kx*Rx + Ky*Ry));

  for(int kx=0;kx<Nx;kx++){
    Kx = (kxmin + kx * dkx);
    argg = -2.0*M_PI*Kx*xmin;
    Px = dx*complex<double>(std::cos(argg),std::sin(argg));

    argg = -2.0*M_PI*Kx*dx;
    fx = complex<double>(std::cos(argg),std::sin(argg));

    
    for(int ky=0;ky<Ny;ky++){
      Ky = (kymin + ky * dky);
      argg = -2.0*M_PI*Ky*ymin;
      Py = dy*complex<double>(std::cos(argg),std::sin(argg));

      argg = -2.0*M_PI*Ky*dy;
      fy = complex<double>(std::cos(argg),std::sin(argg));


      sumx = 0.0; fnx = 1.0;
      for(int nx=0;nx<Nx;nx++){ 

        sumy = 0.0; fny = 1.0;
        for(int ny=0;ny<Ny;ny++){ sumy += in.M[nx*Ny+ny]*fny; fny *= fy; }// for ny

        sumx += sumy*fnx; fnx *= fx;

      }// for nx

                                        
      out.M[kx*Ny+ky] = Px*Py*sumx;

    }// for ky
  }// kx


}

void cfft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy){
/***************************************
  Continuous Fast Fourier Transform for 2D
  
  see derivation in .doc file
 
  The size of the grid should be the power of 2
****************************************/

  // Initial check 
  // in and out are Nx x Ny matrices
  int Nx = in.n_rows;
  int Ny = in.n_cols; 

  
  complex<double> one(0.0,1.0);
  int Nxhalf,Nyhalf;
  int indx,INDX;

  //cout<<"cfft1_2D: Nx = "<<Nx<<endl;
  //cout<<"cfft1_2D: Ny = "<<Ny<<endl;

  if( (Nx>1) && (Nx%2!=0) ){  cout<<"Nx must be power of 2\n"; exit(0);}
  if( (Ny>1) && (Ny%2!=0) ){  cout<<"Ny must be power of 2\n"; exit(0);}

  // Case 1
  if(Nx>2 && Ny>2){

    // Some constants and variables
    Nxhalf = Nx/2;
    Nyhalf = Ny/2;

    complex<double> p1x,p2x,p1y,p2y,eax,eay;

    eax = std::exp(-M_PI*one*xmin/dx); // alpx(2*dx,2*dy)
    p1x = std::exp(-2.0*M_PI*one*kxmin*dx);
    p2x = std::exp(-2.0*M_PI*one/((double)Nx));

    eay = std::exp(-M_PI*one*ymin/dy); // alpy(2*dx,2*dy)
    p2y = std::exp(-2.0*M_PI*one/((double)Ny));


    CMATRIX in_ee(Nxhalf,Nyhalf);
    CMATRIX in_oe(Nxhalf,Nyhalf);
    CMATRIX in_eo(Nxhalf,Nyhalf);
    CMATRIX in_oo(Nxhalf,Nyhalf);

    CMATRIX out_ee(Nxhalf,Nyhalf);
    CMATRIX out_oe(Nxhalf,Nyhalf);
    CMATRIX out_eo(Nxhalf,Nyhalf);
    CMATRIX out_oo(Nxhalf,Nyhalf);

    // Divide set on the even and odd parts
    for(int m1=0;m1<Nxhalf;m1++){
      for(int m2=0;m2<Nyhalf;m2++){

        in_ee.M[m1*Nyhalf+m2] = in.M[(2*m1  )*Ny+(2*m2  )];
        in_oe.M[m1*Nyhalf+m2] = in.M[(2*m1+1)*Ny+(2*m2  )];
        in_eo.M[m1*Nyhalf+m2] = in.M[(2*m1  )*Ny+(2*m2+1)];
        in_oo.M[m1*Nyhalf+m2] = in.M[(2*m1+1)*Ny+(2*m2+1)];

      }// m2
    }// m1


    // Perform fft on smaller grids
    cfft1_2D(in_ee,out_ee,xmin,ymin,kxmin,kymin,2.0*dx,2.0*dy);
    cfft1_2D(in_oe,out_oe,xmin,ymin,kxmin,kymin,2.0*dx,2.0*dy);
    cfft1_2D(in_eo,out_eo,xmin,ymin,kxmin,kymin,2.0*dx,2.0*dy);
    cfft1_2D(in_oo,out_oo,xmin,ymin,kxmin,kymin,2.0*dx,2.0*dy);


    // Compute ft of the original grid using the fts for smaller grids
    for(int kx=0;kx<Nxhalf;kx++){  

      p1y = exp(-2.0*M_PI*one*kymin*dy);        
      for(int ky=0;ky<Nyhalf;ky++){   

        indx = kx*Nyhalf+ky;

        out.M[         kx*Ny + ky]          = 0.25*        (out_ee.M[indx] + p1x*out_oe.M[indx] + p1y*out_eo.M[indx] + p1x*p1y*out_oo.M[indx] );
        out.M[(kx+Nxhalf)*Ny + ky]          = 0.25*eax*    (out_ee.M[indx] - p1x*out_oe.M[indx] + p1y*out_eo.M[indx] - p1x*p1y*out_oo.M[indx] );
        out.M[         kx*Ny + (ky+Nyhalf)] = 0.25*eay*    (out_ee.M[indx] + p1x*out_oe.M[indx] - p1y*out_eo.M[indx] - p1x*p1y*out_oo.M[indx] );
        out.M[(kx+Nxhalf)*Ny + (ky+Nyhalf)] = 0.25*eax*eay*(out_ee.M[indx] - p1x*out_oe.M[indx] - p1y*out_eo.M[indx] + p1x*p1y*out_oo.M[indx] );

        p1y *= p2y;
      }// ky
      p1x *= p2x;
    }// kx

  }// if Nx>2  && Ny>2  : case 1


  // Case 2
  else if(Nx>2 && Ny==2){

    // Some constants and variables
    Nxhalf = Nx/2;

    complex<double> p1x,p2x,eax;

    eax = std::exp(-M_PI*one*xmin/dx); // alpx(2*dx,dy)
    p1x = std::exp(-2.0*M_PI*one*kxmin*dx);
    p2x = std::exp(-2.0*M_PI*one/((double)Nx));


    CMATRIX in_e(Nxhalf,Ny);
    CMATRIX in_o(Nxhalf,Ny);

    CMATRIX out_e(Nxhalf,Ny);
    CMATRIX out_o(Nxhalf,Ny);

    // Divide set on the even and odd parts
    for(int m1=0;m1<Nxhalf;m1++){
      for(int m2=0;m2<Ny;m2++){

        in_e.M[m1*Ny+m2] = in.M[(2*m1  )*Ny+m2];
        in_o.M[m1*Ny+m2] = in.M[(2*m1+1)*Ny+m2];

      }// m2
    }// m1


    // Perform fft on smaller grids
    cfft1_2D(in_e,out_e,xmin,ymin,kxmin,kymin,2.0*dx,dy);
    cfft1_2D(in_o,out_o,xmin,ymin,kxmin,kymin,2.0*dx,dy);


    // Compute ft of the original grid using the fts for smaller grids
    for(int kx=0;kx<Nxhalf;kx++){  
      for(int ky=0;ky<Ny;ky++){   

        indx = kx*Ny+ky;

        out.M[         kx*Ny + ky] = 0.5*    (out_e.M[indx] + p1x*out_o.M[indx] );
        out.M[(kx+Nxhalf)*Ny + ky] = 0.5*eax*(out_e.M[indx] - p1x*out_o.M[indx] );

      }// ky

      p1x *= p2x;
    }// kx

  }// if Nx>2  && Ny==2  : case 2

  // Case 3
  else if(Nx==2 && Ny>2){

    // Some constants and variables
    Nyhalf = Ny/2;

    complex<double> p1y,p2y,eay;

    eay = std::exp(-M_PI*one*ymin/dy); // alpy(dx,2*dy)
    p1y = std::exp(-2.0*M_PI*one*kymin*dy);
    p2y = std::exp(-2.0*M_PI*one/((double)Ny));


    CMATRIX in_e(Nx,Nyhalf);
    CMATRIX in_o(Nx,Nyhalf);

    CMATRIX out_e(Nx,Nyhalf);
    CMATRIX out_o(Nx,Nyhalf);

    // Divide set on the even and odd parts
    for(int m1=0;m1<Nx;m1++){
      for(int m2=0;m2<Nyhalf;m2++){

        in_e.M[m1*Nyhalf+m2] = in.M[m1*Ny+(2*m2)];
        in_o.M[m1*Nyhalf+m2] = in.M[m1*Ny+(2*m2+1)];

      }// m2
    }// m1


    // Perform fft on smaller grids
    cfft1_2D(in_e,out_e,xmin,ymin,kxmin,kymin,dx,2.0*dy);
    cfft1_2D(in_o,out_o,xmin,ymin,kxmin,kymin,dx,2.0*dy);


    // Compute ft of the original grid using the fts for smaller grids
    for(int ky=0;ky<Nyhalf;ky++){   
      for(int kx=0;kx<Nx;kx++){  

        indx = kx*Nyhalf+ky;

        out.M[kx*Ny + ky]          = 0.5*    (out_e.M[indx] + p1y*out_o.M[indx] );
        out.M[kx*Ny + (ky+Nyhalf)] = 0.5*eay*(out_e.M[indx] - p1y*out_o.M[indx] );

      }// kx

      p1y *= p2y;
    }// ky

  }// if Nx==2  && Ny>2  : case 3


  else if(Nx==2 && Ny==2){ // No more recursive levels below this one - do all explicitly 
    cft1_2D(in,out,xmin,ymin,kxmin,kymin,dx,dy);
  }// if Nx==2 && Ny==2 : case 4

  else{ cout<<"Not implemented\n";  exit(0); }


}




void cft2(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk){
/***************************************
  Continuous Fourier Transform
  f(k) = Integral ( f(r) * exp(-2*pi*i*k*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*k*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*k*xmin) * sum ( in[n] * exp(-2*pi*k*dx*n)  )
                                       n
  r_n = xmin + dx * n

  k =  kmin + dk * [k]    [k] - is integer equvalent of k

  Note: Number of real and reciprocal space grids is different
****************************************/

  int N_r =  in.n_elts; // <in> - is 1 x N_r or N_r x 1 CMATRIX
  int N_k = out.n_elts; // <out> - is 1 x N_k or N_k x 1 CMATRIX

  complex<double> f,mul;
  double K,argg;

  for(int k=0;k<N_k;k++){

    K = kmin + k*dk;
    argg = -2.0*M_PI*K*xmin;
    mul = complex<double>(dx*std::cos(argg),dx*std::sin(argg));


    argg = -2.0*M_PI*K*dx;
    f = complex<double>(std::cos(argg),std::sin(argg));

    out.M[k] = 0.0;
    for(int n=0;n<N_r;n++){
      out.M[k] += in.M[n]*mul;
      mul *= f;
    }// for j

  }// for i

}


void inv_cft(CMATRIX& in,CMATRIX& out,double xmin,double dx){
/***************************************
  Inverse Continuous Fourier Transform
  f(n) = Integral ( f(k) * exp(2*pi*i*k*r) * dk ) =

  = sum ( f(k') * exp(2*pi*i*(k'/L)*r_n)) * 1/L =
     k'

  = 1/(N*dx) * sum ( in[k'] * exp(2*pi*i*(k'/L)*(dx*n+ xmin))  ) =
                k'
  = 1/(N*dx) * sum ( in[k'] * exp(2*pi*i*(k'/N)*(n+ xmin/dx))  ) =
                k'

  k = k'/L, L = N*dx - reciprocal length
  r_n = xmin + dx*n
****************************************/

  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,mul;
  double arg;
  double L = N*dx;

  for(int n=0;n<N;n++){

    double r_n = xmin + dx*n;

    out.M[n] = 0.0;
    arg = 2.0*M_PI*r_n/L;

    f = complex<double>(std::cos(arg),std::sin(arg));
    mul = 1.0;

    for(int k=0;k<N;k++){
      out.M[n] += in.M[k]*mul;
      mul *= f;
    }// for j
  }// for i

  arg = (1.0/L);

  for(n=0;n<N;n++){ out.M[n] *= arg; }

}

void inv_cft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx){
/***************************************
  Inverse Continuous Fourier Transform
  f(n) = Integral ( f(k) * exp(2*pi*i*k*r) * dk ) =

  = sum ( f(k') * exp(2*pi*i*(kmin + (k'/L) )*r_n)) * 1/L =
     k'

  = 1/(N*dx) * sum ( in[k'] * exp(2*pi*i*( kmin + (k'/L) )*(dx*n+ xmin))  ) =
                k'

  k = k'/L, L = N*dx - reciprocal length
  r_n = xmin + dx*n

  k' - integer, discrete representation
  k  - double, continuous representation

****************************************/

  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,f1,mul;
  double argg;
  double L = N*dx;   // dk = 1/L,   kmin<= k < kmin + 1/dx

  for(int n=0;n<N;n++){

    double r_n = xmin + dx*n;
    out.M[n] = 0.0;
    
    argg = 2.0*M_PI*r_n*kmin;
    complex<double> f1 = complex<double>(std::cos(argg),std::sin(argg));

    argg = 2.0*M_PI*r_n/L;
    f = complex<double>(std::cos(argg),std::sin(argg));
    mul = f1;

    for(int k=0;k<N;k++){
      out.M[n] += in.M[k]*mul;
      mul *= f;
    }// for j
  }// for i

  argg = (1.0/L);

  for(n=0;n<N;n++){ out.M[n] *= argg; }

}

void inv_cfft1(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx){
/***************************************
  Inverse Continuous Fast Fourier Transform

  f(r) = Integral ( f(k) * exp(2*pi*i*(xmin+n*dr)*(kmin+k)) * dk )
                            

  The size of the grid should be the power of 2
****************************************/

  // Initial check 
  int N = in.n_elts; // N must be 2^n, n - some integer
  if( (N>1) && (N%2!=0) ){  cout<<"N must be power of 2\n"; exit(0);}

  if(N>2){

    // Some constants and variables
    double dk = 1.0/(dx*N);
    int Nhalf = N/2;
    int kp;
    complex<double> one(0.0,1.0);
    complex<double> p1, p2,ea;

    ea = std::exp(M_PI*one*kmin*dx*((double)N)); // alp(dx)
    p1 = std::exp(2.0*M_PI*one*xmin*dk);
    p2 = std::exp(2.0*M_PI*one*dk*dx);

    CMATRIX in_even(1,Nhalf);
    CMATRIX in_odd(1,Nhalf);
    CMATRIX out_even(1,Nhalf);
    CMATRIX out_odd(1,Nhalf);

    // Divide set on the even and odd part
    for(int m=0;m<Nhalf;m++){
      in_even.M[m] = in.M[2*m];
      in_odd.M[m] = in.M[2*m+1];
    }

    // Perform fft on smaller grids
    inv_cfft1(in_even,out_even,xmin,kmin,dx);
    inv_cfft1(in_odd,out_odd,xmin,kmin,dx);


    // Compute ft of the original grid using the fts for smaller grids
    for(int k=0;k<Nhalf;k++){    
      out.M[k] = 0.5*(out_even.M[k] + p1 * out_odd.M[k]);
      out.M[k+Nhalf] = 0.5*ea*(out_even.M[k] - p1 * out_odd.M[k]);
      p1 *= p2;
    }  


  }// if N>=2 

  else if(N==2){ // No more recursive levels below this one - do all explicitly 
    inv_cft1(in,out,xmin,kmin,dx); 
  }


}



void inv_cft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy){
/***************************************
  Inverse Continuous 2-D Fourier Transform
  f(r) = Integral ( f(k) * exp(2*pi*i*k*r) * dr ) =

  r_n = xmin + dx * n

  k = kmin + dk * k

  in and out - are now 2D matrices (columns and rows) - this is for convenience, although they both may be 
               reformulated as 1D matrices (colums or rows)
 
  in - k-space
  out - r-space

****************************************/

  // in and out are Nx x Ny matrices
  int Nx = in.n_rows;
  int Ny = in.n_cols; 

  complex<double> sumx,sumy,fx,fy,Px,Py,fnx,fny;
  complex<double> one(0.0,1.0);

  // Real space box [xmin, xmin+Lx] x [ymin, ymin+Ly]
  double argg;
  double Lx = dx*Nx;
  double Ly = dy*Ny;
  double dkx = 1.0/Lx;
  double dky = 1.0/Ly;
  double Kx,Ky,Rx,Ry;

  // Conceptually:  out.M[nx*Ny+ny] +=  in.M[kx*Ny+ky] * exp(2.0*M_PI*one*(Kx*Rx + Ky*Ry));

  for(int nx=0;nx<Nx;nx++){
    Rx = (xmin + nx * dx);
    argg = 2.0*M_PI*Rx*kxmin;
    Px = dkx*complex<double>(std::cos(argg),std::sin(argg));

    argg = 2.0*M_PI*Rx*dkx;
    fx = complex<double>(std::cos(argg),std::sin(argg));

    
    for(int ny=0;ny<Ny;ny++){
      Ry = (ymin + ny * dy);
      argg = 2.0*M_PI*Ry*kymin;
      Py = dky*complex<double>(std::cos(argg),std::sin(argg));

      argg = 2.0*M_PI*Ry*dky;
      fy = complex<double>(std::cos(argg),std::sin(argg));


      sumx = 0.0; fnx = 1.0;
      for(int kx=0;kx<Nx;kx++){ 

        sumy = 0.0; fny = 1.0;
        for(int ky=0;ky<Ny;ky++){ sumy += in.M[kx*Ny+ky]*fny; fny *= fy; }// for ky

        sumx += sumy*fnx; fnx *= fx;

      }// for kx

                                        
      out.M[nx*Ny+ny] = Px*Py*sumx;

    }// for ny
  }// nx


}


void inv_cfft1_2D(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy){
/***************************************
  Inverse Continuous Fast Fourier Transform for 2D
  
  see derivation in .doc file
 
  The size of the grid should be the power of 2
****************************************/

  // Initial check 
  // in and out are Nx x Ny matrices
  int Nx = in.n_rows;
  int Ny = in.n_cols; 

  
  complex<double> one(0.0,1.0);
  int Nxhalf,Nyhalf;
  int indx;


  if( (Nx>1) && (Nx%2!=0) ){  cout<<"Nx must be power of 2\n"; exit(0);}
  if( (Ny>1) && (Ny%2!=0) ){  cout<<"Ny must be power of 2\n"; exit(0);}

  // Case 1
  if(Nx>2 && Ny>2){

    // Some constants and variables
    Nxhalf = Nx/2;
    Nyhalf = Ny/2;

    complex<double> p1x,p2x,p1y,p2y,eax,eay;

    eax = std::exp(M_PI*one*kxmin*dx*((double)Nx)); // alpx(0.5*dx)
    p1x = std::exp(2.0*M_PI*one*xmin/(dx*(double)Nx));
    p2x = std::exp(2.0*M_PI*one/((double)Nx));

    eay = std::exp(M_PI*one*kymin*dy*((double)Ny)); // alpy(0.5*dy)
    p2y = std::exp(2.0*M_PI*one/((double)Ny));


    CMATRIX in_ee(Nxhalf,Nyhalf);
    CMATRIX in_oe(Nxhalf,Nyhalf);
    CMATRIX in_eo(Nxhalf,Nyhalf);
    CMATRIX in_oo(Nxhalf,Nyhalf);

    CMATRIX out_ee(Nxhalf,Nyhalf);
    CMATRIX out_oe(Nxhalf,Nyhalf);
    CMATRIX out_eo(Nxhalf,Nyhalf);
    CMATRIX out_oo(Nxhalf,Nyhalf);

    // Divide set on the even and odd parts
    for(int m1=0;m1<Nxhalf;m1++){
      for(int m2=0;m2<Nyhalf;m2++){

        in_ee.M[m1*Nyhalf+m2] = in.M[(2*m1  )*Ny+(2*m2  )];
        in_oe.M[m1*Nyhalf+m2] = in.M[(2*m1+1)*Ny+(2*m2  )];
        in_eo.M[m1*Nyhalf+m2] = in.M[(2*m1  )*Ny+(2*m2+1)];
        in_oo.M[m1*Nyhalf+m2] = in.M[(2*m1+1)*Ny+(2*m2+1)];

      }// m2
    }// m1


    // Perform fft on smaller grids
    inv_cfft1_2D(in_ee,out_ee,xmin,ymin,kxmin,kymin,dx,dy);
    inv_cfft1_2D(in_oe,out_oe,xmin,ymin,kxmin,kymin,dx,dy);
    inv_cfft1_2D(in_eo,out_eo,xmin,ymin,kxmin,kymin,dx,dy);
    inv_cfft1_2D(in_oo,out_oo,xmin,ymin,kxmin,kymin,dx,dy);


    // Compute ft of the original grid using the fts for smaller grids
    for(int nx=0;nx<Nxhalf;nx++){  

      p1y = exp(2.0*M_PI*one*ymin/(dy*(double)Ny));     
      for(int ny=0;ny<Nyhalf;ny++){   

        indx = nx*Nyhalf+ny;

        out.M[         nx*Ny + ny]          = 0.25*        (out_ee.M[indx] + p1x*out_oe.M[indx] + p1y*out_eo.M[indx] + p1x*p1y*out_oo.M[indx] );
        out.M[(nx+Nxhalf)*Ny + ny]          = 0.25*eax*    (out_ee.M[indx] - p1x*out_oe.M[indx] + p1y*out_eo.M[indx] - p1x*p1y*out_oo.M[indx] );
        out.M[         nx*Ny + (ny+Nyhalf)] = 0.25*eay*    (out_ee.M[indx] + p1x*out_oe.M[indx] - p1y*out_eo.M[indx] - p1x*p1y*out_oo.M[indx] );
        out.M[(nx+Nxhalf)*Ny + (ny+Nyhalf)] = 0.25*eax*eay*(out_ee.M[indx] - p1x*out_oe.M[indx] - p1y*out_eo.M[indx] + p1x*p1y*out_oo.M[indx] );

        p1y *= p2y;
      }// ky
      p1x *= p2x;
    }// kx

  }// if Nx>2  && Ny>2  : case 1


  // Case 2
  else if(Nx>2 && Ny==2){

    // Some constants and variables
    Nxhalf = Nx/2;

    complex<double> p1x,p2x,eax;
    eax = std::exp(M_PI*one*kxmin*dx*((double)Nx)); // alpx(0.5*dx)
    p1x = std::exp(2.0*M_PI*one*xmin/(dx*(double)Nx));
    p2x = std::exp(2.0*M_PI*one/((double)Nx));


    CMATRIX in_e(Nxhalf,Ny);
    CMATRIX in_o(Nxhalf,Ny);

    CMATRIX out_e(Nxhalf,Ny);
    CMATRIX out_o(Nxhalf,Ny);

    // Divide set on the even and odd parts
    for(int m1=0;m1<Nxhalf;m1++){
      for(int m2=0;m2<Ny;m2++){

        in_e.M[m1*Ny+m2] = in.M[(2*m1  )*Ny+m2];
        in_o.M[m1*Ny+m2] = in.M[(2*m1+1)*Ny+m2];

      }// m2
    }// m1


    // Perform fft on smaller grids
    inv_cfft1_2D(in_e,out_e,xmin,ymin,kxmin,kymin,dx,dy);
    inv_cfft1_2D(in_o,out_o,xmin,ymin,kxmin,kymin,dx,dy);


    // Compute ft of the original grid using the fts for smaller grids
    for(int kx=0;kx<Nxhalf;kx++){  
      for(int ky=0;ky<Ny;ky++){   

        indx = kx*Ny+ky;

        out.M[         kx*Ny + ky] = 0.5*    (out_e.M[indx] + p1x*out_o.M[indx] );
        out.M[(kx+Nxhalf)*Ny + ky] = 0.5*eax*(out_e.M[indx] - p1x*out_o.M[indx] );

      }// ky

      p1x *= p2x;
    }// kx

  }// if Nx>2  && Ny==2  : case 2



  // Case 3
  else if(Nx==2 && Ny>2){

    // Some constants and variables
    Nyhalf = Ny/2;

    complex<double> p1y,p2y,eay;

    eay = std::exp(M_PI*one*kymin*dy*((double)Ny)); // alpy(0.5*dy)
    p1y = std::exp(2.0*M_PI*one*ymin/(dy*(double)Ny));     
    p2y = std::exp(2.0*M_PI*one/((double)Ny));

    CMATRIX in_e(Nx,Nyhalf);
    CMATRIX in_o(Nx,Nyhalf);

    CMATRIX out_e(Nx,Nyhalf);
    CMATRIX out_o(Nx,Nyhalf);

    // Divide set on the even and odd parts
    for(int m1=0;m1<Nx;m1++){
      for(int m2=0;m2<Nyhalf;m2++){

        in_e.M[m1*Nyhalf+m2] = in.M[m1*Ny+(2*m2)];
        in_o.M[m1*Nyhalf+m2] = in.M[m1*Ny+(2*m2+1)];

      }// m2
    }// m1


    // Perform fft on smaller grids
    inv_cfft1_2D(in_e,out_e,xmin,ymin,kxmin,kymin,dx,dy);
    inv_cfft1_2D(in_o,out_o,xmin,ymin,kxmin,kymin,dx,dy);


    // Compute ft of the original grid using the fts for smaller grids
    for(int ky=0;ky<Nyhalf;ky++){   
      for(int kx=0;kx<Nx;kx++){  

        indx = kx*Nyhalf+ky;

        out.M[kx*Ny + ky]          = 0.5*    (out_e.M[indx] + p1y*out_o.M[indx] );
        out.M[kx*Ny + (ky+Nyhalf)] = 0.5*eay*(out_e.M[indx] - p1y*out_o.M[indx] );

      }// kx

      p1y *= p2y;
    }// ky

  }// if Nx==2  && Ny>2  : case 3


  else if(Nx==2 && Ny==2){ // No more recursive levels below this one - do all explicitly 
    inv_cft1_2D(in,out,xmin,ymin,kxmin,kymin,dx,dy);
  }// if Nx==2 && Ny==2 : case 4

  else{ cout<<"Not implemented\n";  exit(0); }


}






void inv_cft2(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk){
/***************************************
  Inverse Continuous Fourier Transform
  f(r) = Integral ( f(k) * exp(2*pi*i*k*r) * dk ) =

  = sum ( f(k') * exp(2*pi*i*k*r)) * dk =
     k

  = dk * exp(2*pi*i*( kmin * r_n)) * sum ( in[k] * exp(2*pi*i*dk*[k]*r_n)  ) = f[n]
                                     [k]

  k = kmin + dk * [k]  [k] - is integer equivalent of k
  r_n = xmin + dx*n

****************************************/

  int N_k =  in.n_elts; // <in> - is 1 x N_k or N_k x 1 CMATRIX
  int N_r = out.n_elts; // <out> - is 1 x N_r or N_r x 1 CMATRIX

  complex<double> f,f1,mul;
  double r_n,argg;

  for(int n=0;n<N_r;n++){

    r_n = xmin + dx*n;
    argg = 2.0*M_PI*r_n*kmin;
    mul = complex<double>(dk*std::cos(argg),dk*std::sin(argg));

    argg = 2.0*M_PI*r_n*dk;
    f = complex<double>(std::cos(argg),std::sin(argg));

    out.M[n] = 0.0;
    for(int k=0;k<N_k;k++){
      out.M[n] += in.M[k]*mul;
      mul *= f;

    }// for k
  }// for n


}



void convolve(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx){
// Convolve two Fourier transforms
//  conv(k) = integral(  f(k') * g(k-k') ) dk' = sum (n'/L) * f[n'] * g[n-n'] = conv[n]
//                                               n'
// n, n' - integers,  k = n/L, k' = n'/L
// L = dx * N, N - size of the set

  int N = f.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> G;
  double L = N*dx;  L = (1.0/L);  // dk = 1/(dx*N),  kmin <= k < kmin + 1/dx

  for(int n=0;n<N;n++){
    conv.M[n] = 0.0;

    for(int np=0;np<N;np++){

      if((n-np)>=0 ){ G = g.M[n-np]; }
      else{ G = g.M[n-np+N]; }

      conv.M[n] += f.M[np]*G*L;

    }// for np
  }// for n

}  


void convolve_2D(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx,double dy){
// Convolve two Fourier transforms. Each in 2D
//  conv(k) = integral(  f(k') * g(k-k') ) dk' = sum (n'/L) * f[n'] * g[n-n'] = conv[n]
//                                               n'
// n, n' - integers,  k = n/L, k' = n'/L
// L = dx * N, N - size of the set


  // in and out are Nx x Ny matrices
  int Nx = f.n_rows;
  int Ny = f.n_cols; 

  complex<double> sum;

  // Real space box [xmin, xmin+Lx] x [ymin, ymin+Ly]
  double Lx = dx*Nx;
  double Ly = dy*Ny;
  double dkx = 1.0/Lx; // kxmin <= kx < kxmin + 1/dx
  double dky = 1.0/Ly; // kymin <= ky < kymin + 1/dy
  double dV = dkx*dky;

  int nx,ny;



  for(int kx=0;kx<Nx;kx++){
    for(int ky=0;ky<Ny;ky++){

      sum = 0.0;      
      for(int kxp=0;kxp<Nx;kxp++){
        nx = (kx - kxp); if(nx<0){ nx += Nx; }

        for(int kyp=0;kyp<Ny;kyp++){
          ny = (ky - kyp); if(ny<0){ ny += Ny; }

          sum += f.M[kxp*Ny+kyp]*g.M[nx*Ny+ny]; // F*G(kx,ky)  =  F[kx',ky'] * G[kx-kx',ky-ky']

        } // ky'
      }// kx'

          conv.M[kx*Ny+ky] = sum * dV;

    }// for ky
  }// kx

}  

}// namespace liblinalg
}// namespace libmmath

