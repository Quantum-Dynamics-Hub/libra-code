#include "SpecialFunctions.h"
//using libmmath::liblinalg;

//================== Functions ==========================
// Teylor expansions for some formula are from:
// http://web.mit.edu/kenta/www/three/taylor.html

namespace libmmath{
namespace libspecialfunctions{


double FAST_POW(double x,int n){

    double res,res1;    
    if(n==0)      { res = 1.0; }
    else if(n==1) { res = x; }                                    // 0 mult
    else if(n==2) { res = x*x;}                                   // 1 mult
    else if(n==3) { res = x*x*x; }                                // 2 mult
    else if(n==4) { res = x*x; res = res*res; }                   // 2 mult
    else if(n==5) { res = x*x; res = res*res*x; }                 // 3 mult
    else if(n==6) { res = x*x*x; res = res*res; }                 // 3 mult
    else if(n==7) { res = x*x*x; res = res*res*x; }               // 4 mult
    else if(n==8) { res = x*x; res = res*res; res = res*res;}     // 3 mult
    else if(n==9) { res = x*x*x; res = res*res*res; }             // 4 mult
    else if(n==10){ res1= res = x*x; 
                    res = res*res;
                    res = res*res;
                    res = res*res1;
                  }                                               // 4 mult
    else if(n==11){ res1= res = x*x;
                    res = res*res;
                    res = res*res;
                    res = res*res1*x;
                  }                                               // 5 mult
    else if(n==12){ res = x*x*x;
                    res = res*res;
                    res = res*res;
                  }                                               // 4 mult
    else{   res = pow(x,n);   }

    return res;

}


double sinh_(double x){
/*  This is sinh(x)/x function  */

    double tol = 1e-12;
    double res,x2,x4,x6,x8;
    double c1 = (1.0/6.0);
    double c2 = (1.0/120.0);
    double c3 = (1.0/5040.0);
    double c4 = (1.0/362880);

    if(fabs(x)>tol){  res = (sinh(x)/x); }
    else{ // Use Maclaurin series
          x2 = x*x;
          x4 = x2*x2;
          x6 = x2*x4;
          x8 = x2*x6;

          res = 1.0 + c1*x2 + c2*x4 + c3*x6 + c4*x8;
    }

    return res;

}

double sin_(double x){
/*  This is sin(x)/x function */

    double tol = 1e-12;
    double res,x2,x4,x6,x8,x10;
    double c1 = (1.0/6.0);
    double c2 = (1.0/120.0);
    double c3 = (1.0/5040.0);
    double c4 = (1.0/362880.0);
    double c5 = (1.0/39916800.0);
    double c6 = (1.0/6227020800.0);

    if(fabs(x)>tol){  res = (sin(x)/x); }
    else{ // Use Maclaurin series
          x2 = x*x;
          x4 = x2*x2;
          x6 = x2*x4;
          x8 = x2*x6;
          x10 = x2*x6;

          res = 1.0 - c1*x2 + c2*x4 - c3*x6 + c4*x8 - c5*x10;
    }

    return res;

}

double ERF(double x){

/* Error function approximation:

Source:
http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf

Article:
S. Winitzki, Uniform approximations for transcendental functions, in Proc. ICCSA-
2003, LNCS, 2667/2003, p. 962

/Theory/erf-approx.pdf

*/

   int signx = 1;
   double x2,ax2;

   if(x<=0) { x = -x; signx = -1; }

   x2 = x * x;
   ax2 = 0.147*x2;

   double res = sqrt(1.0 - exp(-x2*(1.273239545 + ax2)/(1.0 + ax2)) );

   res = res*signx;

   return res;

}

double ERFC(double x){

/* Error function approximation:

Source:
http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf

Article:
S. Winitzki, Uniform approximations for transcendental functions, in Proc. ICCSA-
2003, LNCS, 2667/2003, p. 962

/Theory/erf-approx.pdf

*/

   int signx = 1;
   double x2,ax2;

   if(x<=0) { x = -x; signx = -1; }

   x2 = x * x;
   ax2 = 0.147*x2;

   double res = sqrt(1.0 - exp(-x2*(1.273239545 + ax2)/(1.0 + ax2)) );

   res = 1.0 - res*signx;

   return res;

}

double FACTORIAL(int n) {

        if(n<0){
                cout<<"Factorial of negative number is not defined!"<<endl;
        }
        else{

                if(n==1||n==0){
                        return 1.0;
                }
                else if(n>1){
                        return double(n*FACTORIAL(n-1));
                }
        }
}

int DFACTORIAL(int n) {

// This is n!! = n * (n-2)!!
  int res; 

  if(n<=1){ res = 1;}
  else if(n==3){ res = 3; }
  else if(n==5){ res = 15; }
  else if(n==7){ res = 105; }
  else if(n==9){ res = 945; }
  else if(n==11){ res = 10395; }
  else if(n==13){ res = 135135; }

  else{  res = n*DFACTORIAL(n-2); }

  return res;
}


double BINOM(int i,int n) {

  double a,b,c,res;

  if(n==0){ 
    if(i==0){ res = 1.0; }
  }
  else if(n==1){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 1.0; }
  }
  else if(n==2){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 2.0; }
    else if(i==2){ res = 1.0; }
  }
  else if(n==3){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 3.0; }
    else if(i==2){ res = 3.0; }
    else if(i==3){ res = 1.0; }
  }
  else if(n==4){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 4.0; }
    else if(i==2){ res = 6.0; }
    else if(i==3){ res = 4.0; }
    else if(i==4){ res = 1.0; }
  }
  else if(n==5){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 5.0; }
    else if(i==2){ res = 10.0; }
    else if(i==3){ res = 10.0; }
    else if(i==4){ res = 5.0; }
    else if(i==5){ res = 1.0; }
  }
  else if(n==6){
    if(i==0){ res = 1.0; }
    else if(i==1){ res = 6.0; }
    else if(i==2){ res = 15.0; }
    else if(i==3){ res = 20.0; }
    else if(i==4){ res = 15.0; }
    else if(i==5){ res = 6.0; }
    else if(i==6){ res = 1.0; }
  }

  else{

    a = FACTORIAL(n);
    b = FACTORIAL(i);
    c = FACTORIAL(n-i);

    res = (a/(b*c));
  }

  return res;
}

void LEGENDRE(int n,double x,double a,double b,double& p,double& q){

        // Legendre polynomials and its derivatives
        // Algorithm is taken from: Flowers, B. H. "An Introduction to Numerical Methods in C++", Clarendon Press, Oxford, 1995
        /*

        <> - subscript

                recursive relation: (n+1)P<n+1>(x) = (2n+1)xP<n>(x) - nP<n-1>(x)
                normaliztion:        = 2/(2n+1)
        recursive derivatives: P'<n+1>(x) = xP'<n>(x) + (n+1)P<n>(x)    
                weighting function: w(x) = 1
        
            defined for an interval [a,b]

        */

        x = 0.5*((b+a)+(b-a)*x);

        double r = 0, s = 0, t = 0;

        // r - value of polynomial one degree  less
        // t - value of polynomial two degrees less
        // s - differential coefficient corresponding to r

        p = 1; q = 0;

        for(float m=0;m<n;m++){

                r = p; s = q;
                p = (2.0*m+1)*x*r - m*t;
                p /= m+1;
                q = x*s + (m+1)*r;
                t = r;
        }

}

void CHEBYSHEV(int n,double x,double a,double b,double& p,double& q){

        // Chebyshev polynomials and its derivatives    
        /*

        <> - subscript

                recursive relation: T<n+1>(x) = 2xT<n>(x) - T<n-1>(x)
                normaliztion:        = PI (n==0);  PI/2 (n>0)
        recursive derivatives: T'<n+1>(x) = [(n+1)x/n]T'<n>(x) + (n+1)T<n>(x)   
                weighting function: w(x) = 1/sqrt(1-x^2)
        
            defined for an interval [a,b]

        */

        x = 0.5*((b+a)+(b-a)*x);

        double r = 0, s = 0, t = x;

        // r - value of polynomial one degree  less
        // t - value of polynomial two degrees less
        // s - differential coefficient corresponding to r


        p = 1; q = 0;

        for(float m=0;m<n;m++){

                r = p; s = q;
                p = 2.0*x*r - t;
                if(m==0){
                        q = (m+1)*r;
                }
                else{
                q = ((m+1)/m)*x*s + (m+1)*r;
                }
                t = r;
        }

}

void LAGUERRE(int n,double x,double& p,double& q){

        // Laguerre polynomials and its derivatives     
        /*

        <> - subscript

                recursive relation: (n+1)L<n+1>(x) = (2n+1-x)L<n>(x) - nL<n-1>(x)
                normaliztion:        = 1
        recursive derivatives: L'<n+1>(x) = L'<n>(x) - L<n>(x)  
                weighting function: w(x) = exp(-x)
        
            defined for an interval [0,+inf)

        */

        double r = 0, s = 0, t = 0;

        // r - value of polynomial one degree  less
        // t - value of polynomial two degrees less
        // s - differential coefficient corresponding to r

        p = 1; q = 0;

        for(float m=0;m<n;m++){

                r = p; s = q;
                p = (2.0*m+1-x)*r - m*t;
                p /= m+1;
                q = s - r;
                t = r;
        }

}

void HERMITE(int n,double x,double& p,double& q){

        // Hermite polynomials and its derivatives      
        /*

        <> - subscript

                recursive relation: H<n+1>(x) = 2xH<n>(x) - 2nH<n-1>(x)
                normaliztion:        = 2^n*n!*sqrt(PI)
        recursive derivatives: H'<n+1>(x) = 2(n+1)H<n>(x)       
                weighting function: w(x) = exp(-x^2)
        
            defined for an interval (-inf,+inf)

        */

        double r = 0, s = 0, t = 0;

        // r - value of polynomial one degree  less
        // t - value of polynomial two degrees less
        // s - differential coefficient corresponding to r

        p = 1; q = 0;

        for(float m=0;m<n;m++){

                r = p; s = q;
                p = 2.0*x*r - 2.0*m*t;
                q = 2.0*(m+1)*r;
                t = r;
        }

}

double Ellipe(double phi0,double k,int N)  // Elliptic integral of first kind (Jacobi elliptic function)
{

        double cp,sp,sp2;
        double k2 = k*k;
        double phi=phi0;
        double S,M,M1,ADD;

        cp = cos(phi);
        sp = sin(phi);
        sp2= sp*sp;

        S  = 0.5*(phi-cp*sp);
        M  = 0.5*k2;
        M1 = cp*sp*sp2;

        double sum=phi;

        ADD = M*S;

        for(int n=1;n<N;n++){


                sum+=ADD;

                S = (1.0/(2.0*n+2.0))*(-M1+(2.0*n+1.0)*S);

                M1 = M1*sp2;

                M *= ((2.0*(n+1.0)-1.0)*k2/(2.0*(n+1.0)));

                ADD = M*S;


        }

        //sum = phi/(div*a1);

        return sum;
}
void Ellipe2(double u,double k,double& sn,double& cn,double& dn){

        double k2,k4,k6;
        double u2,u4,u6,u8;
        double f2,f3,f4,f5,f6,f7,f8;

        f2 = 0.5;
        f3 = f2/3.0;
        f4 = f3/4.0;
        f5 = f4/5.0;
        f6 = f5/6.0;
        f7 = f6/7.0;
        f8 = f7/8.0;

        u2 = u*u;
        u4 = u2*u2;
        u6 = u2*u4;
        u8 = u2*u6;

        k2 = k*k;
        k4 = k2*k2;
        k6 = k2*k4;


        sn = u * (1.0 - (1.0+k2)*u2*f3 + (1.0+14.0*k2+k4)*u4*f5 - (1.0+135.0*k2+135.0*k4+k6)*u6*f7);

    cn = sqrt(1.0-sn*sn);

        dn = sqrt(1.0-k*k*sn*sn);

}
void Jacobi_Elliptic(double u,double m,double tol,double& am, double& sn, double& cn, double& dn){

        // [A. S.] = M. Abramowitz, I.A. Stegun, Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables, Dover, New
    //York, 1965

        double K,Kprime;

        double k = sqrt(m);

        K      = Km(m,tol);
        Kprime = Km(1-m,tol);

        double q  = exp(-M_PI*(Kprime/K)); // nome
        double q2 = q*q;

        double v = (M_PI*u/(2.0*K));

        double q_sqrt = sqrt(q);
        double q_n, q_2n;
        double term,term1,term2,term3;
        double mt1, mt2, mt3;

        am = v;

        q_n  = q;
        q_2n = q2;

        double n = 1;
        while(true){

                // [A. S. 17.3.24]
                term = 2.0*q_n*sin(2.0*n*v)/(n*(1+q_2n));

                am += term;

                n++;
                q_n *= q;
                q_2n *= q2;

                if(fabs(term)<tol){ break; }

        }

        sn = sin(am);          // [A.S. 17.2.2]
        cn = cos(am);          // [A.S. 17.2.3]
        dn = sqrt(1-m*sn*sn);  // [A.S. 17.2.4]

}
double Km(double m,double tol){

        // Return function:
        // K(m) = F(PI/2|m) = complete elliptic integral of first kind;

        double a0,b0,c0;
        double a1,b1,c1;
        double k   = sqrt(m); // m = k*k;


        a1 = a0 = 1;
        b1 = b0 = sqrt(1-m);
        c1 = c0 = sqrt(m);

        double n = 1;

        while(true){

                a1 = 0.5*(a0+b0);
                b1 = sqrt(a0*b0);
                c1 = 0.5*(a0-b0);

                if(fabs(c1*n)<tol) break;

                a0 = a1;
                b0 = b1;
                c0 = c1;

                n++;

        }

//        cout<<"m = "<<m<<endl;
//        cout<<"# of iterations = "<<n<<endl;

        return (0.5*M_PI/a1);

}

void Ellint(double m,double sinphi,double tol,double& Km_,double& value){
// This function calculates incomplete elliptic integral of first kind (returned in value)
// It also returns the complete elliptic integral of first kind (in Km_)
// Arguments: sinphi - sine of angle phi
// m - parameter, tol - accuracy 

//------ Precompute some quantities and the complete integral -------------
  if (m < 0.0) {
    fprintf(stderr, "ERROR: m < 0.0 in Ellint (m = %f): FIX IT!\n", m);
    exit(1);
  }

  int MAXN = 25;
  double* a; a = new double[MAXN];
  double* b; b = new double[MAXN];
  double* d; d = new double[MAXN];

  int N = 0;
  a[N] = 1.0;
  b[N] = sqrt(1.0-m);
  d[N] = b[N]/a[N];
  double err;

  while(N<MAXN){
    N++;
    a[N] = 0.5*(a[N-1] + b[N-1]);
    b[N] = sqrt(a[N-1]*b[N-1]);
    d[N] = b[N]/a[N];
    err = (a[N-1]-b[N-1]);
    
    if(fabs(err*(double)N)<tol) break;
  }

  Km_      = (0.5*M_PI/a[N]);
//  cout<<"m = "<<m<<endl;
//  cout<<"# of iterations = "<<N<<endl;

  double Kscale = (1.0/((1<<N)*a[N]));

//---------- Now calculate incomplete integral -----------
  if (sinphi <= -1.0){   value =  -Km_; }
  else if (sinphi >= 1.0){   value = Km_; }
  else{
    double t = (sinphi/sqrt(1.0-sinphi*sinphi));
    int s = 0;
    for (int i = 0; i < N; i++) {
      double c = 1.0 - d[i]*t*t;
      s *= 2;
      if (c < 0.0){ s += (t>0.0?1:-1);}
      t = (1.0 + d[i])*t/c;
    }
    value = Kscale*(atan(t)+s*M_PI);
  }// else

  delete [] a;
  delete [] b;
  delete [] d;

}


int randperm(int size,int of_size,vector<int>& result){
// Makes a random permutation of 'of_size' numbers and place 'size'
// of them into 'result' vector

  int* all_range;

  all_range = new int[of_size];

  // Initialize
  for(int i=0;i<of_size;i++){
     all_range[i] = i;
  }
  // Make permutations
  int indx;
  int lowest=0, highest=of_size-1;
  int range=(highest-lowest)+1;
  int tmp;

  for(int i=0;i<size;i++){

    indx = rand();
    indx = i+lowest+int((range-i)*(indx/(RAND_MAX + 1.0)));

    // Swap all_range[i] and all_range[indx]
    tmp = all_range[i];
    all_range[i] = all_range[indx];
    all_range[indx] = tmp;

  }
  // Now get 'size' first elements
  if(result.size()>0) {result.clear();}
  for(int i=0;i<size;i++){
     result.push_back(all_range[i]);
  }

  delete [] all_range;

  return 0;
}

MATRIX exp_(MATRIX& x,double dt){
  // This function calculates exp(x*dt)
  // C.T() * x * C = D
  MATRIX res,C,D; C = 0.0; D = 0.0;
  x.JACOBY_EIGEN(D, C, 1e-20);

  D.M[0] = exp(dt*D.M[0]);
  D.M[4] = exp(dt*D.M[4]);
  D.M[8] = exp(dt*D.M[8]);

  res = C * D * C.T();
  return res;
}

MATRIX exp1_(MATRIX& x,double dt){
  // This function calculates exp(x*dt)*sinh(x*t)/(x*t)
  // where 1/x is x^-1
  // C.T() * x * C = D
  MATRIX res,C,D,inv_x; C = 0.0; D = 0.0;
  x.JACOBY_EIGEN(D, C, 1e-20);

  D.M[0] = exp(dt*D.M[0])*sinh_(dt*D.M[0]);
  D.M[4] = exp(dt*D.M[4])*sinh_(dt*D.M[4]);
  D.M[8] = exp(dt*D.M[8])*sinh_(dt*D.M[8]);

  res = C * D * C.T();

  return res;
}

MATRIX3x3 exp_(MATRIX3x3& x,double dt){
  // This function calculates exp(x*dt)
  // C.T() * x * C = D
  MATRIX3x3 res,C,D; C = 0.0; D = 0.0;
  x.eigen(D, C, 1e-20);

  D.xx = exp(dt*D.xx);
  D.yy = exp(dt*D.yy);
  D.zz = exp(dt*D.zz);

  res = C * D * C.T();
  return res;
}

MATRIX3x3 exp1_(MATRIX3x3& x,double dt){
  // This function calculates exp(x*dt)*sinh(x*t)/(x*t)
  // where 1/x is x^-1
  // C.T() * x * C = D
  MATRIX3x3 res,C,D,inv_x; C = 0.0; D = 0.0;
  x.eigen(D, C, 1e-20);

  D.xx = exp(dt*D.xx)*sinh_(dt*D.xx);
  D.yy = exp(dt*D.yy)*sinh_(dt*D.yy);
  D.zz = exp(dt*D.zz)*sinh_(dt*D.zz);

  res = C * D * C.T();

  return res;
}



//============================ Random number generators =================
double RANDOM(double a,double b){

        double u;
               u = rand();
                   u =u/RAND_MAX;

        return (a+(b-a)*u);
}


void MATRIX_TO_QUATERNION(MATRIX& M,QUATERNION& L){
/***********************************************************
    Convert orientation matrix to corresponding quaternion 
                     representation
***********************************************************/

//----------- Quaternion initialization ---------
  double a,b,c,d,e,f;
  //****************************************************
  // from www.gamedev.ru/users/wat/articles/quaternions
  // *****************************************************/
                      
  double a00,a01,a02,a10,a11,a12,a20,a21,a22;
  double m[3][3];
  double ind_q[4];

  m[0][0] = a00 = M.M[0];
  m[0][1] = a01 = M.M[1];
  m[0][2] = a02 = M.M[2];
  m[1][0] = a10 = M.M[3];
  m[1][1] = a11 = M.M[4];
  m[1][2] = a12 = M.M[5];
  m[2][0] = a20 = M.M[6];
  m[2][1] = a21 = M.M[7];
  m[2][2] = a22 = M.M[8];

  int ind_i,ind_j,ind_k;
  int nxt[3] = {1,2,0};
  double  tr,s;

  tr = a00 + a11 + a22;
  
  if(tr>0.0){
     s = sqrt(tr + 1.0);
     L.Lt = s/2.0; 
  
     s = 0.5/s;

     L.Lx = (a12 - a21)*s;
     L.Ly = (a20 - a02)*s;
     L.Lz = (a01 - a10)*s;
  }else{
     ind_i = 0;
     if(a11>a00)              ind_i = 1;
     if(a22>m[ind_i][ind_i])  ind_i = 2;

     ind_j  = nxt[ind_i];
     ind_k  = nxt[ind_j];

     s = sqrt(m[ind_i][ind_i] - (m[ind_j][ind_j]+m[ind_k][ind_k])+1.0 );

     ind_q[ind_i] = s * 0.5;

     if(s != 0.0)  s = 0.5/s;

     ind_q[3]      = (m[ind_j][ind_k] - m[ind_k][ind_j]) * s;
     ind_q[ind_j]  = (m[ind_i][ind_j] + m[ind_j][ind_i]) * s;
     ind_q[ind_k]  = (m[ind_i][ind_k] + m[ind_k][ind_i]) * s;

     L.Lx = ind_q[0];
     L.Ly = ind_q[1];
     L.Lz = ind_q[2];
     L.Lt = ind_q[3];

  }


}


void MATRIX_TO_QUATERNION(MATRIX3x3& M,QUATERNION& L){
/***********************************************************
    Convert orientation matrix to corresponding quaternion
                     representation
***********************************************************/

//----------- Quaternion initialization ---------
  double a,b,c,d,e,f;
  //****************************************************
  // from www.gamedev.ru/users/wat/articles/quaternions
  // *****************************************************/

  double a00,a01,a02,a10,a11,a12,a20,a21,a22;
  double m[3][3];
  double ind_q[4];

  m[0][0] = a00 = M.xx;
  m[0][1] = a01 = M.xy;
  m[0][2] = a02 = M.xz;
  m[1][0] = a10 = M.yx;
  m[1][1] = a11 = M.yy;
  m[1][2] = a12 = M.yz;
  m[2][0] = a20 = M.zx;
  m[2][1] = a21 = M.zy;
  m[2][2] = a22 = M.zz;

  int ind_i,ind_j,ind_k;
  int nxt[3] = {1,2,0};
  double  tr,s;

  tr = a00 + a11 + a22;

  if(tr>0.0){
     s = sqrt(tr + 1.0);
     L.Lt = s/2.0;

     s = 0.5/s;

     L.Lx = (a12 - a21)*s;
     L.Ly = (a20 - a02)*s;
     L.Lz = (a01 - a10)*s;
  }else{
     ind_i = 0;
     if(a11>a00)              ind_i = 1;
     if(a22>m[ind_i][ind_i])  ind_i = 2;

     ind_j  = nxt[ind_i];
     ind_k  = nxt[ind_j];

     s = sqrt(m[ind_i][ind_i] - (m[ind_j][ind_j]+m[ind_k][ind_k])+1.0 );

     ind_q[ind_i] = s * 0.5;

     if(s != 0.0)  s = 0.5/s;

     ind_q[3]      = (m[ind_j][ind_k] - m[ind_k][ind_j]) * s;
     ind_q[ind_j]  = (m[ind_i][ind_j] + m[ind_j][ind_i]) * s;
     ind_q[ind_k]  = (m[ind_i][ind_k] + m[ind_k][ind_i]) * s;

     L.Lx = ind_q[0];
     L.Ly = ind_q[1];
     L.Lz = ind_q[2];
     L.Lt = ind_q[3];
  }
}


void QUATERNION_TO_MATRIX(QUATERNION& L,MATRIX& M){
/**************************************************************
   Convert quaternion to corresponding orientation matrix 
                        representation
***************************************************************/
   double q0,q1,q2,q3;

   q0 = L.Lt;
   q1 = L.Lx;
   q2 = L.Ly;
   q3 = L.Lz;

   M.M[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
   M.M[1] = 2.0*(q1*q2 + q0*q3);
   M.M[2] = 2.0*(q1*q3 - q0*q2);
   M.M[3] = 2.0*(q1*q2 - q0*q3);
   M.M[4] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
   M.M[5] = 2.0*(q2*q3 + q0*q1);
   M.M[6] = 2.0*(q1*q3 + q0*q2);
   M.M[7] = 2.0*(q2*q3 - q0*q1);
   M.M[8] = q0*q0 - q1*q1 - q2*q2 + q3*q3;


}

void QUATERNION_TO_MATRIX(QUATERNION& L,MATRIX3x3& M){
/**************************************************************
   Convert quaternion to corresponding orientation matrix
                        representation
***************************************************************/
   double q0,q1,q2,q3;

   q0 = L.Lt;
   q1 = L.Lx;
   q2 = L.Ly;
   q3 = L.Lz;

   M.xx = q0*q0 + q1*q1 - q2*q2 - q3*q3;
   M.xy = 2.0*(q1*q2 + q0*q3);
   M.xz = 2.0*(q1*q3 - q0*q2);
   M.yx = 2.0*(q1*q2 - q0*q3);
   M.yy = q0*q0 - q1*q1 + q2*q2 - q3*q3;
   M.yz = 2.0*(q2*q3 + q0*q1);
   M.zx = 2.0*(q1*q3 + q0*q2);
   M.zy = 2.0*(q2*q3 - q0*q1);
   M.zz = q0*q0 - q1*q1 - q2*q2 + q3*q3;


}


//=========================== Just some function ================================
    
    
void solve_linsys(MATRIX& C,MATRIX& D, MATRIX& X,double eps,int maxiter){
/*********************************************
 Here we solve the system of linear equations
              AX = D
 using Gauss-Seidel iterative procedure

 Inputs: A, D - matrices
         eps  - precision criterion
 Output: X

 Some preliminary transformations are made in order
 to be able to use Gauss-Seidel method for any matrix A

 More details:
 80.47 An iterative Algorithm for Matrix Inversion
 Which is Always Convergent
 Authors: S. Simons
 Source: The Mathematical Gazette, Vol. 80, No. 489
 (Nov., 1996), pp. 567-569

 url: http://www.jstor.org/stable/pdfplus/3618529.pdf

**********************************************/

// Do the transformations A = C^T * C and b = C^T * d
// If matrices d and c have more then 1 columns we do the
// procedure for each column


    int i,j,k;   // counters
    int n,m,p;   // dimetions
    double s;    // sums
    double error;// error
    int iter;    // number of iterations

    if(C.num_of_rows!=D.num_of_rows)
        {std::cout<<"Error: The number of rows of matrices C and D in equation CX = D must be equal\n"; exit(35); } // n
    if(C.num_of_cols!=X.num_of_rows)
        {std::cout<<"Error: The number of cols of matrix C and num of rows in matrix D in equation CX = D must be equal\n"; exit(35); } // m
    if(X.num_of_cols!=D.num_of_cols)
        {std::cout<<"Error: The number of cols of matrices X and D in equation CX = D must be equal\n"; exit(35); } // p

    // Set dimentions
    n = C.num_of_rows;
    m = C.num_of_cols;
    p = D.num_of_cols;

    MATRIX A(m,m); A = C.T() * C;
    X = 0.0;
    error = 2.0*eps;
    iter = 0;

    while((error>eps)&&(iter<maxiter)){

    error = 0.0;

    for( k = 0; k < p; k++ ){

        //------- Matrix preparation step -----------

        MATRIX d(n,1);  for(i = 0;i<n;i++){ d.M[i] = D.M[i*p+k]; }
        MATRIX x(m,1);  for(i = 0;i<m;i++){ x.M[i] = X.M[i*p+k]; }
        MATRIX xprev(m,1); xprev = x;
        MATRIX b(m,1);  b = C.T() * d;

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

            error += (x.M[i] - xprev.M[i])*(x.M[i] - xprev.M[i]);

            X.M[i*p + k] = x.M[i];

        }// for i


    }// for k - all columns of D

    error = sqrt(error/double(m));
//    cout<<"Iteration "<<iter<<" error = "<<error<<endl;

    iter++;

    }// loop over convergence

    if(error>eps){
      cout<<"Error in solve_linsys: convergence to eps= "<<eps<<" is not achieved for "<<iter<<" iterations\n";
      exit(0);
    }

//    cout<<"Convergence achieved after "<<iter<<" iterations\n";


}



int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out){

  if(out.size()>0){ out.clear(); }
  int sz = in.size();
  if(sz==0){ }
  else if(sz==1){ out = in; }
  else{
    // Divide in into 2 blocks of approximately same size each
    int half = sz/2;
    vector< pair<int,double> > in1,in2,out1,out2;
    for(int i=0;i<half;i++){ in1.push_back(in[i]); }
    merge_sort(in1,out1);
       
    for(int i=half;i<sz;i++){ in2.push_back(in[i]); }
    merge_sort(in2,out2); 
      
    // Now merge two parts
    int cl,cr; cl = 0; cr = half;
    while((cl<half) && (cr<sz)){
 
    /// EXTREMELY IMPORTANT !!!  This simple, slight difference - the use of < or <= makes HUGE differnece
    /// The "good" version maximally preserves the ordering of orbitals, so one does not run into trouble of alternating
    /// charges - this also leads to symmetric charge distribution in unrestricted formulations even without population smearing
    /// I think it is even more than that - this can lead to convergence (or faster convergence), while the wrong
    /// method may lead to either non-convergent scheme or to sifnificantly slower convergenc.

    if(out1[cl].second<=out2[cr-half].second){ out.push_back(out1[cl]); cl++; }  ///< <-- This is good
//    if(out1[cl].second<out2[cr-half].second){ out.push_back(out1[cl]); cl++; }  ///< <-- Try old

    /// The "bad" version will alternate order of nearby orbitals
    /// It is here only for the purpose of "demonstration of pathological implementation"
//    if(out1[cl].second<out2[cr-half].second){ out.push_back(out1[cl]); cl++; } <-- This is BAD
      else{ out.push_back(out2[cr-half]); cr++; }
    } 
    while(cl<half){ out.push_back(out1[cl]); cl++;}
    while(cr<sz)  { out.push_back(out2[cr-half]); cr++;}
  }
  return 0;

}// int merge_sort(vector< pair<int,double> >& in, vector< pair<int,double> >& out)



}// namespace libspecialfunctions
}// namespace libmmath



