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

#include "FT.h"


/// liblibra 
namespace liblibra{

/// liblinalg namespace
namespace liblinalg{



void solve_linsys(CMATRIX& C,CMATRIX& D, CMATRIX& X,double eps,int maxiter,double omega){
/**
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
*/

/// Do the transformations A = C^H * C and b = C^H * d
/// If matrices d and c have more then 1 columns we do the
/// procedure for each column


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
      for( j = 0; j < i; j++ ){ s += A.M[i*m + j]*x.M[j]; }
      for( j = i+1; j < m; j++ ){ s += A.M[i*m + j]*xprev.M[j]; }
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
/**
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

*/

/// Do the transformations A = C^H * C and b = C^H * d
/// If matrices d and c have more then 1 columns we do the
/// procedure for each column


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
        for( j = 0; j < i; j++ ){ s += A.M[im + j]*x.M[j];  }
        for( j = i+1; j < m; j++ ){ s += A.M[im + j]*xprev.M[j]; }
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
/**
  Discrete Fourier Transform
  e.g. http://en.wikipedia.org/wiki/Fast_Fourier_transform

  \param[in] in The input matrix
  \param[out] out The output matrix

  "in" and "out: are the vectors with n elements: n x 1
*/

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
/**
  Inverse Discrete Fourier Transform
  e.g. http://en.wikipedia.org/wiki/Discrete_Fourier_transform

  \param[in] in The input matrix
  \param[out] out The output matrix

  "in" and "out: are the vectors with n elements: n x 1

*/

  int k;
  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,mul;
  double argg;

  for(k=0;k<N;k++){

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
/**
  Continuous Fourier Transform
  f(k) = Integral ( f(r) * exp(-2*pi*i*k*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*k*r_n)) * dx =
     n
  = dx *exp(-2*pi*i*k*xmin) * sum ( in[n] * exp(-2*pi*k*dx*n)  )
                               n
  r_n = xmin + dx * n

  \param[in] in The input matrix
  \param[out] out The output matrix
  \param[in] xmin The minimal boundary of the real-space grid (in)
  \param[in] dx Spacing between points in the real space 

  "in" and "out: are the vectors with n elements: n x 1

*/

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
/**
  Continuous Fourier Transform
  f(kmin + k) = Integral ( f(r) * exp(-2*pi*i*(kmin+k)*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*(kmin+k)*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*(kmin+k)*xmin) * sum ( in[n] * exp(-2*pi*(kmin+k)*dx*n)  )
                                       n
  r_n = xmin + dx * n

  \param[in] in The input matrix
  \param[out] out The output matrix
  \param[in] xmin The minimal boundary of the real-space grid (in)
  \param[in] kmin The minimal boundary of the reciprocal-space grid (out)
  \param[in] dx Spacing between points in the real space 

  "in" and "out: are the vectors with n elements: n x 1



*/

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
/**
  Continuous Fast Fourier Transform

  f(kmin + k) = Integral ( f(r) * exp(-2*pi*i*(kmin+k)*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*(kmin+k)*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*(kmin+k)*xmin) * sum ( in[n] * exp(-2*pi*(kmin+k)*dx*n)  )
                                       n
  r_n = xmin + dx * n

  The size of the grid should be the power of 2

  \param[in] in The input matrix
  \param[out] out The output matrix
  \param[in] xmin The minimal boundary of the real-space grid (in)
  \param[in] kmin The minimal boundary of the reciprocal-space grid (out)
  \param[in] dx Spacing between points in the real space 

  "in" and "out: are the vectors with n elements: n x 1

*/

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
/**
  Continuous 2-D Fourier Transform
  f(kmin + k) = Integral ( f(r) * exp(-2*pi*i*(kmin+k)*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*(kmin+k)*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*(kmin+k)*xmin) * sum ( in[n] * exp(-2*pi*(kmin+k)*dx*n)  )
                                       n
  r_n = xmin + dx * n

  in and out - are now 2D matrices (columns and rows) - this is for convenience, although they both may be 
               reformulated as 1D matrices (colums or rows)

  \param[in] in The input matrix
  \param[out] out The output matrix
  \param[in] xmin The minimal boundary of the real-space grid along X axis
  \param[in] ymin The minimal boundary of the real-space grid along Y axis
  \param[in] kxmin The minimal boundary of the reciprocal-space grid along X direction
  \param[in] kymin The minimal boundary of the reciprocal-space grid along Y direction
  \param[in] dx Spacing between points in the real space, along the X direction
  \param[in] dy Spacing between points in the real space, along the Y direction


  "in" and "out: are Nx x Ny matrices

*/

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
/**
  Continuous Fast Fourier Transform for 2D
  
  see derivation in .doc file
 
  The size of the grid should be the power of 2

  \param[in] in The input matrix
  \param[out] out The output matrix
  \param[in] xmin The minimal boundary of the real-space grid along X axis
  \param[in] ymin The minimal boundary of the real-space grid along Y axis
  \param[in] kxmin The minimal boundary of the reciprocal-space grid along X direction
  \param[in] kymin The minimal boundary of the reciprocal-space grid along Y direction
  \param[in] dx Spacing between points in the real space, along the X direction
  \param[in] dy Spacing between points in the real space, along the Y direction


  "in" and "out: are Nx x Ny matrices

*/

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
/**
  Continuous Fourier Transform
  f(k) = Integral ( f(r) * exp(-2*pi*i*k*r) * dr ) =

  = sum ( f(r_n) * exp(-2*pi*i*k*r_n)) * dx =
     n

  = dx * exp(-2*pi*i*k*xmin) * sum ( in[n] * exp(-2*pi*k*dx*n)  )
                                       n
  r_n = xmin + dx * n

  k =  kmin + dk * [k]    [k] - is integer equvalent of k

  Note: Number of real and reciprocal space grids is different

  \param[in] in The input matrix
  \param[out] out The output matrix
  \param[in] xmin The minimal boundary of the real-space grid along X axis
  \param[in] kxmin The minimal boundary of the reciprocal-space grid along X direction
  \param[in] dx Spacing between points in the real space, along the X direction
  \param[in] dk Spacing between points in the reciprocal space, along the X direction


  "in" is 1 x Nr or Nr x 1 matrix , Nr - the number of grid points in the real space
  "out" is 1 x Nk or Nk x 1 matrix, Nk - the number of grid points in the reciprocal space

*/

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
/**
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

  \param[in] in The input matrix
  \param[out] out The output matrix
  \param[in] xmin The minimal boundary of the real-space grid along X axis
  \param[in] dx Spacing between points in the real space, along the X direction

*/
  int n;
  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,mul;
  double arg;
  double L = N*dx;

  for(n=0;n<N;n++){

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
/**
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

*/
  int n;
  int N = in.n_elts; // <in> and <out> are the vectors with n elements: n x 1
  complex<double> f,f1,mul;
  double argg;
  double L = N*dx;   // dk = 1/L,   kmin<= k < kmin + 1/dx

  for(n=0;n<N;n++){

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
/**
  Inverse Continuous Fast Fourier Transform

  f(r) = Integral ( f(k) * exp(2*pi*i*(xmin+n*dr)*(kmin+k)) * dk )
                            

  The size of the grid should be the power of 2
*/

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
/**
  Inverse Continuous 2-D Fourier Transform
  f(r) = Integral ( f(k) * exp(2*pi*i*k*r) * dr ) =

  r_n = xmin + dx * n

  k = kmin + dk * k

  in and out - are now 2D matrices (columns and rows) - this is for convenience, although they both may be 
               reformulated as 1D matrices (colums or rows)
 
  in - k-space
  out - r-space

*/

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
/**
  Inverse Continuous Fast Fourier Transform for 2D
  
  see derivation in .doc file
 
  The size of the grid should be the power of 2
*/

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
/**
  Inverse Continuous Fourier Transform
  f(r) = Integral ( f(k) * exp(2*pi*i*k*r) * dk ) =

  = sum ( f(k') * exp(2*pi*i*k*r)) * dk =
     k

  = dk * exp(2*pi*i*( kmin * r_n)) * sum ( in[k] * exp(2*pi*i*dk*[k]*r_n)  ) = f[n]
                                     [k]

  k = kmin + dk * [k]  [k] - is integer equivalent of k
  r_n = xmin + dx*n

*/

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
/**
  Convolve two Fourier transforms
   conv(k) = integral(  f(k') * g(k-k') ) dk' = sum (n'/L) * f[n'] * g[n-n'] = conv[n]
                                                n'
  n, n' - integers,  k = n/L, k' = n'/L
  L = dx * N, N - size of the set
  <in> and <out> are the vectors with n elements: n x 1
*/

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
/**
  Convolve two Fourier transforms. Each in 2D
   conv(k) = integral(  f(k') * g(k-k') ) dk' = sum (n'/L) * f[n'] * g[n-n'] = conv[n]
                                                n'
  n, n' - integers,  k = n/L, k' = n'/L
  L = dx * N, N - size of the set

  in and out are Nx x Ny matrices
*/

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



} /// liblinalg namespace
} /// liblibra

