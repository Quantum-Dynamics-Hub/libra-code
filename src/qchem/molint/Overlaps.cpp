#include "Overlaps.h"
#include "A_coefficients.h"

namespace libqchem{
namespace libmolint{




double gaussian_overlap_ref(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb ){
/****************************************************************************
 This function computes the 1D overlap: this is the reference method, compact, but not optimized
 this version also does not implement derivatives and normalizaition
 Gaussians here are not normalized

 <g_a(x)|g_b(x)>  - just the case of 1D Gaussians

 Remember also:

 exp(-alp_a*(x-Xa)) * exp(-alp_b*(x-Xb)) = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * exp(-gamma*(x-Xp)^2)

 where Xp = (alp_a*Xa + alp_b * Xb)/gamma

 gamma = alp_a + alp_b

 So:
                                         nx
 (x-Xa)^nxa = (x - Xp + Xp - Xa)^nxa =  sum  [  C_nxa^i * (x-Xp)^i * (Xp-Xa)^(nxa-i) ]
                                       i = 0

                                        nxb
 (x-Xb)^nxb = (x - Xp + Xp - Xb)^nxa =  sum  [  C_nxb^j * (x-Xp)^j * (Xp-Xb)^(nxb-j) ]
                                       j = 0

                                                 nxa     nxb
 <g_a|g_b> = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * sum    sum [ C_nxa^i * C_nxb^j * (Xp-Xa)^i * (Xp-Xb)^j * { (x-Xp)^(nxa-i+nxb-j) * exp(-gamma*(x-Xp)^2) } ]
                                                 i = 0   j=0

****************************************************************************/

  double gamma = alp_a + alp_b;
  double Xp = (alp_a*Xa + alp_b*Xb)/gamma;
  double Xpa = Xp - Xa;
  double Xpb = Xp - Xb;


  double res = 0.0;

  double pwi = 1.0;
  for(int i=0;i<=nxa;i++){
    double bini = BINOM(i,nxa);

    double pwj = 1.0;
    for(int j=0;j<=nxb;j++){
      double binj = BINOM(j,nxb);


      res += bini * binj * pwi * pwj * gaussian_int((nxa-i+nxb-j), gamma);


      pwj *= Xpb;
    }// for j

    pwi *= Xpa;
  }// for i

  res *= exp(-alp_a*alp_b*(Xa-Xb)*(Xa-Xb)/gamma);


  return res;

}// gaussian_overlap_ref





double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize,
                        int is_derivs, double& dI_dXa,double& dI_dXb,
                        vector<double*>& aux,int n_aux ){
/****************************************************************************
 This function computes the 1D overlap: and its derivatives w.r.t. input parameters (positions) - optional
 Computation of derivatives is controlled by is_derivs flag
 This version is optimized - by passing vector to memory pointer (so no need to do allocation every time the function is called)
 The function can also optionally assume the gaussians are normalized (is_normalize==1) or assume gaussians as they are

 <g_a(x)|g_b(x)>  - just the case of 1D Gaussians

 Remember also:

 exp(-alp_a*(x-Xa)) * exp(-alp_b*(x-Xb)) = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * exp(-gamma*(x-Xp)^2)

 where Xp = (alp_a*Xa + alp_b * Xb)/gamma

 gamma = alp_a + alp_b

 So:
                                         nx
 (x-Xa)^nxa = (x - Xp + Xp - Xa)^nxa =  sum  [  C_nxa^i * (x-Xp)^i * (Xp-Xa)^(nxa-i) ]
                                       i = 0

                                        nxb
 (x-Xb)^nxb = (x - Xp + Xp - Xb)^nxa =  sum  [  C_nxb^j * (x-Xp)^j * (Xp-Xb)^(nxb-j) ]
                                       j = 0

                                                 nxa     nxb
 <g_a|g_b> = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * sum    sum [ C_nxa^i * C_nxb^j * (Xp-Xa)^i * (Xp-Xb)^j * { (x-Xp)^(nxa-i+nxb-j) * exp(-gamma*(x-Xp)^2) } ]
                                                 i = 0   j=0

                                   +inf
  For the derivative:  if   G(n) = int  {  (x-Xp)^n * exp(-gamma*(x-Xp)^2) * dx
                                   -inf

  then:  dG(n)/dXp = -n*G(n-1) - 2.0*gamma*G(n+1)

****************************************************************************/


  int i;
  double gamma = alp_a + alp_b;
  double ag = alp_a/gamma;
  double bg = alp_b/gamma;

  double Xp = ag*Xa + bg*Xb;
  double Xpa = Xp - Xa;
  double Xpb = Xp - Xb;

  // Jacobian elements:
  double dXpa_dXa = ag - 1.0;
  double dXpb_dXa = ag;
  double dXpa_dXb = bg;
  double dXpb_dXb = bg - 1.0;
  double dXp_dXa = ag;
  double dXp_dXb = bg;


  // Aliaces
  double* f;  f = aux[0];
  double* dfdXpa; dfdXpa = aux[1];
  double* dfdXpb; dfdXpb = aux[2];
  double* g;  g = aux[3];       // will containg gaussian integral

  // Compute binomial expansion and derivatives, if necessary
  binomial_expansion(nxa, nxb, Xpa, Xpb, f, dfdXpa, dfdXpb, is_derivs); // (x+x1)^n1 * (x+x2)^n2 = summ_i { x^i * f_i (x1,x2,n1,n2) }

  // Compute Gaussian integrals
  for(i=0;i<=(nxa+nxb+1);i++){ 
    g[i] = gaussian_int(i, gamma);
  }


  // Now we are ready to put everything together and compute the overall integral, and its derivatives w.r.t. original variables
  double I = 0.0;
  dI_dXa = 0.0;
  dI_dXb = 0.0;


  for(int i=0;i<=(nxa+nxb);i++){
    I += f[i] * g[i];

    if(is_derivs){
      double dg_dXp = -2.0*gamma*g[i+1];  // dg[i]/dXp = -i*g[i-1] - 2.0*gamma*G[i+1]
      if(i>=1){ dg_dXp -= i*g[i-1]; }

      dI_dXa += ( dfdXpa[i] * dXpa_dXa + dfdXpb[i] * dXpb_dXa ) * g[i]; // + f[i] * dg_dXp * dXp_dXa;
      dI_dXb += ( dfdXpa[i] * dXpa_dXb + dfdXpb[i] * dXpb_dXb ) * g[i]; // + f[i] * dg_dXp * dXp_dXb;

    }// is_derivs
  }

  double pref = exp(-alp_a*alp_b*(Xa-Xb)*(Xa-Xb)/gamma);
  double res = I * pref;

  if(is_derivs){  
    dI_dXa *= pref;
    dI_dXb *= pref;

    dI_dXa += (-2.0*alp_a*alp_b*(Xa-Xb)/gamma)*res;
    dI_dXb -= (-2.0*alp_a*alp_b*(Xa-Xb)/gamma)*res;
  }

  // In case we need to normalize initial Gaussians
  if(is_normalize){

    double nrm = gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b);
    res *= nrm;
    dI_dXa *= nrm; 
    dI_dXb *= nrm; 

  }

  return res;

}// gaussian_overlap

double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize,
                        int is_derivs, double& dI_dXa,double& dI_dXb
                       ){

  // Allocate working memory
  int i;
  int n_aux = 20; //nxa+nxb+1;
  vector<double*> auxd(5);
  for(i=0;i<5;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_overlap(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb, auxd, n_aux);

  // Clean working memory
  for(i=0;i<5;i++){ delete [] auxd[i]; }  
  auxd.clear();
 
  return res;

}// gaussian_overlap


boost::python::list gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                                     int is_normalize, int is_derivs ){
  double dI_dXa, dI_dXb;
  double I = gaussian_overlap(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, is_derivs, dI_dXa, dI_dXb);

  boost::python::list res;

  res.append(I);
 
  if(is_derivs){
    res.append(dI_dXa);
    res.append(dI_dXb);
  }

  return res;
 
}



double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                        int is_normalize
                       ){
  double dI_dxa, dI_dxb;
  double res = gaussian_overlap(nxa,alp_a,Xa, nxb,alp_b,Xb, is_normalize, 0, dI_dxa, dI_dxb);
  return res;

}// gaussian_overlap


double gaussian_overlap(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb){

  double res = gaussian_overlap(nxa,alp_a,Xa, nxb,alp_b,Xb, 1);
  return res;

}// gaussian_overlap







double gaussian_overlap_fast1(int nxa,double alp_a, double Xa, int nxb,double alp_b, double Xb,
                             int is_normalize,
                             int is_derivs, double& dI_dxa,double& dI_dxb,
                             vector<double*>& aux,int n_aux ){

/****************************************************************************
 is_normalize = 1 - multiply results by the normalization constants of two Gaussians
              = 0 - just use the result as is

 This function computes the 1D overlap: - this is practical version

 <g_a(x)|g_b(x)>  - just the case of 1D Gaussians

 Remember also:

 exp(-alp_a*(x-Xa)) * exp(-alp_b*(x-Xb)) = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * exp(-gamma*(x-Xp)^2)

 where Xp = (alp_a*Xa + alp_b * Xb)/gamma

 gamma = alp_a + alp_b

 So:
                                         nx
 (x-Xa)^nxa = (x - Xp + Xp - Xa)^nxa =  sum  [  C_nxa^i * (x-Xp)^i * (Xp-Xa)^(nxa-i) ]
                                       i = 0

                                        nxb
 (x-Xb)^nxb = (x - Xp + Xp - Xb)^nxa =  sum  [  C_nxb^j * (x-Xp)^j * (Xp-Xb)^(nxb-j) ]
                                       j = 0

                                                 nxa     nxb
 <g_a|g_b> = exp(-alp_a*alp_b*(Xa-Xb)^2/gamma) * sum    sum [ C_nxa^i * C_nxb^j * (Xp-Xa)^i * (Xp-Xb)^j * { (x-Xp)^(nxa-i+nxb-j) * exp(-gamma*(x-Xp)^2) } ]
                                                 i = 0   j=0


  The Gaussian integral is directly expressed:

   +inf
  int   {   x^2m * exp(-alp*x^2) dx }   = ( (2m-1)! /( a^(m/2) * 2^(m/2) ) ) * sqrt(pi/a) (zero for 2m+1 power)
   -inf

****************************************************************************/

  double* f;  f = aux[0];
  double* dfdxa; dfdxa = aux[1];
  double* dfdxb; dfdxb = aux[2];
  double* PWa; PWa = aux[3];
  double* dPWa; dPWa = aux[4];
  double* PWb; PWb = aux[5];
  double* dPWb; dPWb = aux[6];

  double gamma = alp_a + alp_b;
  double Xp = (alp_a*Xa + alp_b*Xb)/gamma;
  double Xpa = Xp - Xa;
  double Xpb = Xp - Xb;


  for(int k=0;k<=(nxa+nxb);k++){ f[k] = 0.0; dfdxa[k] = 0.0; dfdxb[k] = 0.0; }

  // Pre-compute PW1 and PW2 arrays independently:
  double pw,dpw,bin;
  
  for(int i=nxa;i>=0;i--){

    int nia = nxa - i;

    if(nia==0){ pw = 1.0; dpw = 0.0;    }
    else{  
      if(is_derivs){ dpw = pw + Xpa*dpw; }
      pw *= Xpa;
    } // order matters here!!!
   
    bin = BINOM(i,nxa);

    PWa[i] = bin * pw;  
    if(is_derivs){ dPWa[i] = bin * dpw;}

  }

  for(int i=nxb;i>=0;i--){

    int nib = nxb - i;

    if(nib==0){ pw = 1.0; dpw = 0.0;    }
    else{  
      if(is_derivs){dpw = pw + Xpb*dpw;}
      pw *= Xpb;
    } // order matters here!!!
   
    bin = BINOM(i,nxb);

    PWb[i] = bin * pw;  
    if(is_derivs){ dPWb[i] = bin * dpw;}

  }

  // Now do the double loop
  for(int i=0;i<=nxa;i++){
    for(int j=0;j<=nxb;j++){

      f[i+j] += PWa[i]*PWb[j];
      if(is_derivs){
        dfdxa[i+j] += dPWa[i]*PWb[j];
        dfdxb[i+j] += PWa[i]*dPWb[j];
      }

    }// for j
  }// for i

  

  double I = 0.0;
  dI_dxa = 0.0;
  dI_dxb = 0.0;


  double rat = (0.5/gamma);
  double w = sqrt(M_PI/gamma);
  for(int k=0; k<=(nxa+nxb); k+=2){   // this implies that k%2 = 0 - always

    if(k>0){   w *= (k-1) * rat;  }   // w = gaussian_int((nxa-i+nxb-j), gamma);

    I += f[k] * w;
    if(is_derivs){
      dI_dxa += dfdxa[k] * w;
      dI_dxb += dfdxb[k] * w;
    }
  }

  I *= exp(-alp_a*alp_b*(Xa-Xb)*(Xa-Xb)/gamma);



  if(is_normalize){
    double nrm = gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b);
    dI_dxa *= nrm;
    dI_dxb *= nrm;
    I *= nrm;
  }

  return I;

}// gaussian_overlap_fast



double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb
                       ){

  // This works even if Gaussians are not normalized
  VECTOR dIdA, dIdB;
  double res = gaussian_overlap(nxa,nya,nza,alp_a,Ra,nxb,nyb,nzb,alp_b,Rb,1,dIdA,dIdB);

  return res;
}

double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        VECTOR& dIdA, VECTOR& dIdB
                       ){
  // Allocate working memory
  int i;
  int n_aux = 40;
  vector<double*> auxd(20);
  for(i=0;i<20;i++){ auxd[i] = new double[n_aux]; }

  // Do computations
  double res = gaussian_overlap(nxa,nya,nza,alp_a,Ra, nxb,nyb,nzb,alp_b,Rb, is_normalize, dIdA, dIdB, auxd, n_aux);

  // Clean working memory
  for(i=0;i<20;i++){ delete [] auxd[i]; }  
  auxd.clear();
 

  return res;

}


double gaussian_overlap(int nxa,int nya, int nza, double alp_a, VECTOR& Ra,
                        int nxb,int nyb, int nzb, double alp_b, VECTOR& Rb,
                        int is_normalize, 
                        VECTOR& dIdA, VECTOR& dIdB,
                        vector<double*>& auxd,int n_aux
                       ){
// Compute the overlap integral of two Gaussians GA and GB (each with arbitrary angular momentum)
// is_normalize - controls if the Gaussians should be normalized:
//   is_normalize = 0 - use Gaussians as they are
//   is_normalize = 1 - use normalize the overlap as if the Gaussians are normalized


  double Ix = gaussian_overlap(nxa, alp_a, Ra.x, nxb, alp_b, Rb.x);
  double Iy = gaussian_overlap(nya, alp_a, Ra.y, nyb, alp_b, Rb.y);
  double Iz = gaussian_overlap(nza, alp_a, Ra.z, nzb, alp_b, Rb.z);


  if(is_normalize){
    Ix *= (gaussian_norm(nxa,alp_a) * gaussian_norm(nxb,alp_b));
    Iy *= (gaussian_norm(nya,alp_a) * gaussian_norm(nyb,alp_b));
    Iz *= (gaussian_norm(nza,alp_a) * gaussian_norm(nzb,alp_b));
  }

  double I = Ix * Iy * Iz;
  if(fabs(I)<1e-15){ I = 0.0; }

  dIdA = 0.0;
  dIdB = 0.0;

  return I;

}// gaussian_overla


double sto_norm(int n, double alp){
// In fact, this is not the norm, but the normalization coefficient
/*
  Meaning: if N - is a result of this function and
  S(n,l,m,alpha) = (r-center)^(n-1) * exp(-alp * (r - center))  is a STO
  then
  s = N * S(l,m,n,alp) - is normaized: integral(s,s) = 1.0
 
  e.g. see: http://en.wikipedia.org/wiki/Slater-type_orbital
*/


  double res = sqrt(2.0*alp/FACTORIAL(2*n))*FAST_POW(2.0*alp,n);  // this is normalization factor

  return res;
}


double sto_overlap(int na, int la, int ma, double alp_a, int nb, int lb, int mb, double alp_b, 
                   double R, double R_cutoff){
// Compute the overlap integral of two STO: (na, la, ma) and (nb, lb, mb) separated by distance R

  double I;

  if(ma == mb){
 
    int m = abs(ma);

    if(R<R_cutoff){  

      // Exponents
      double rhoA = R * alp_a;
      double rhoB = R * alp_b;
      double sum = rhoA + rhoB;
      double dif = rhoA - rhoB;
    
      double ra = 2.0*alp_a/(alp_a + alp_b); //(rhoA/sum);
      double rb = 2.0*alp_b/(alp_a + alp_b); //(rhoB/sum);
    
      double pref = sqrt(ra*rb) * FAST_POW(ra,na) * FAST_POW(rb,nb);
    
      // Allocate memory for coefficients
      double* f;  f = new double[na+nb+1];
      double* g;  g = new double[na+nb+1];
    
    
      Aux_F1(rhoA, rhoB, f, na+nb);       
        
      I = 0.0;      
      for(int mu=0;mu<=(na+nb);mu++){
    
        g[mu] = 0.0;
        int v1 = max(max(0,abs(la-lb)-mu),mu-(la+lb));
        int v2 = na+nb - m;
    
        for(int v=v1;v<=v2;v++){
          g[mu] += FAST_POW(sum,v)*A_coefficient_general(mu,v, na, nb, la, lb, m);
        }
    
        I += FAST_POW(dif,mu) * g[mu] * f[mu];
    
      }// for mu
    
      // NOT ENTIRELY SURE HOW TO DEAL WITH THE SIGN
      if(R==0.0){
        I *= pow(-1.0,la+ma); // According to the beginning of paper
      }

      // Deallocate memory
      delete [] f;
      delete [] g;

    } // R<R_cutoff


  }// if ma==mb
  else{ I = 0.0; }

  return I;

}// sto_overlap



double sto_overlap_fast(int na, int la, int ma, double alp_a, int nb, int lb, int mb, double alp_b, 
                   double R, double R_cutoff){
// Compute the overlap integral of two STO: (na, la, ma) and (nb, lb, mb) separated by distance R

  double I;

  if(ma == mb){
 
    int m = abs(ma);

    if(R<R_cutoff){  

      // Exponents
      double rhoA = R * alp_a;
      double rhoB = R * alp_b;
      double sum = rhoA + rhoB;
      double dif = rhoA - rhoB;
    
      double ra = 2.0*alp_a/(alp_a + alp_b); //(rhoA/sum);
      double rb = 2.0*alp_b/(alp_a + alp_b); //(rhoB/sum);
    
      double pref = sqrt(ra*rb) * FAST_POW(ra,na) * FAST_POW(rb,nb);
    
      // Allocate memory for coefficients
      double** A_coeff; 
      A_coeff = new double*[10];
      for(int i=0;i<10;i++){  A_coeff[i] = new double[10]; }
      double* f;  f = new double[na+nb+1];
      double* g;  g = new double[na+nb+1];
    
    
      Aux_F1(rhoA, rhoB, f, na+nb);       
      A_coefficients_fast(0,0, na, nb, la, lb, m, A_coeff);
    
    
      I = 0.0;
    
    
      for(int mu=0;mu<=(na+nb);mu++){
    
        g[mu] = 0.0;
        int v1 = max(max(0,abs(la-lb)-mu),mu-(la+lb));
        int v2 = na+nb - m;
    
        for(int v=v1;v<=v2;v++){
          g[mu] += FAST_POW(sum,v)*A_coeff[mu][v];
        }
    
        I += FAST_POW(dif,mu) * g[mu] * f[mu];
    
      }// for mu
    
      // NOT ENTIRELY SURE HOW TO DEAL WITH THE SIGN
      if(R==0.0){
        I *= pow(-1.0,la+ma); // According to the beginning of paper
      }

      // Deallocate memory
      delete [] f;
      delete [] g;
      for(int i=0;i<10;i++){ delete [] A_coeff[i]; }
      delete [] A_coeff;

    } // R<R_cutoff


  }// if ma==mb
  else{ I = 0.0; }

  return I;

}// sto_overlap_fast





}// namespace libmolint
}// namespace libqchem


