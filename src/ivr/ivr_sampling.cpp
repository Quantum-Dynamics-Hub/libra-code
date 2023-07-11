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
  \file ivr_sampling.cpp
  \brief These are the C++ re-implementation of the Fortran codes from the Ananth group.
         The original codes can be found here:
         https://github.com/AnanthGroup/SC-IVR-Code-Package    

  According to original documentation:

! This file contains MonteCarlo subroutines that
! generate an array of initial phase space points
! for the multidimensional SC-IVRs.

! The form of the provided sampling distribution is 
! commented in each subroutine, but we define operator A
! as the projection of an initial coherent state,

!       A = |pIn,qIn><pIn,qIn|

! We use a Gaussian RNG



*/

#include "ivr.h"


/// liblibra namespace
namespace liblibra{


namespace libivr{

MATRIX ivr_Husimi(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd){
/**
  \brief Multi-dimensional  Husimi-IVR MonteCarlo
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] Random - random number generator object
      
  Samples one point from the distribution:     w = |<p0,q0|pIn,qIn>|^2

  Return value: a Ndof x 2 matrix with q (column 0) and p (column 1)

*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_Husimi: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  MATRIX res(Ndof,2);  ///< first column is position, second is momentum

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));
    res.set(i, 0, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q
    res.set(i, 1,    s    * rnd.normal() + pIn.get(i,0) );  // p

  }// for i

  return res;
}


vector<MATRIX> ivr_Husimi(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size){
/**
  \brief Multi-dimensional  Husimi-IVR MonteCarlo
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] Random - random number generator object
  \param[in] sample_size - how many points we want to sample
      
  Samples many points from the distribution:     w = |<p0,q0|pIn,qIn>|^2

  Return value: a vector of Ndof x 2 matrices with q (column 0) and p (column 1)
*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_Husimi: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  vector<MATRIX> res(sample_size, MATRIX(Ndof,2));  ///< first column is position, second is momentum

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));

    for(int j=0; j<sample_size; j++){
        res[j].set(i, 0, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q
        res[j].set(i, 1,    s    * rnd.normal() + pIn.get(i,0) );  // p

    }// for j - all sampling points
  }// for i - all dofs

  return res;

}



MATRIX ivr_LSC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd){
/**
  \brief Multi-dimensional  LSC-IVR MonteCarlo 
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] Random - random number generator object
      
  Samples one point from the distribution:     w = |<p0,q0|pIn,qIn>|^2

  Return value: a Ndof x 2 matrix with q (column 0) and p (column 1)

  This implementation is exactly the same as Husimi-IVR

*/

  return ivr_Husimi(qIn, pIn, Width0, rnd);

}

vector<MATRIX> ivr_LSC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size){
/**
  \brief Multi-dimensional  LSC-IVR MonteCarlo
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] Random - random number generator object
  \param[in] sample_size - how many points we want to sample
      
  Samples many points from the distribution:     w = |<p0,q0|pIn,qIn>|^2

  Return value: a vector of Ndof x 2 matrices with q (column 0) and p (column 1)

  This implementation is exactly the same as Husimi-IVR
*/

  return ivr_Husimi(qIn, pIn, Width0, rnd, sample_size);
}



MATRIX ivr_DHK(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd){
/**
  \brief Multi-dimensional  DHK-IVR MonteCarlo 
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] Random - random number generator object
      
  Samples one point from the distribution:  w = |<p0,q0|pIn,qIn><pIn,qIn|p0',q0'>|^2

  Return value: a Ndof x 4 matrix with q (column 0) and p (column 1), q' (column 2), p' (column 3)

*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_DHK: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  MATRIX res(Ndof,4);  

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));
    res.set(i, 0, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q0
    res.set(i, 1,    s    * rnd.normal() + pIn.get(i,0) );  // p0
    res.set(i, 2, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q0'
    res.set(i, 3,    s    * rnd.normal() + pIn.get(i,0) );  // p0'


  }// for i

  return res;

}


vector<MATRIX> ivr_DHK(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size){
/**
  \brief Multi-dimensional  DHK-IVR MonteCarlo 
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] Random - random number generator object
  \param[in] sample_size - how many points we want to sample
      
  Samples many points from the distribution:     w = |<p0,q0|pIn,qIn><pIn,qIn|p0',q0'>|^2

  Return value: a vector of Ndof x 4 matrices with q (column 0) and p (column 1), q' (column 2), p' (column 3)
*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_DHK: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  vector<MATRIX> res(sample_size, MATRIX(Ndof,4));  ///< (q0, p0, q0', p0')

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));

    for(int j=0; j<sample_size; j++){
        res[j].set(i, 0, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q
        res[j].set(i, 1,    s    * rnd.normal() + pIn.get(i,0) );  // p
        res[j].set(i, 2, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q
        res[j].set(i, 3,    s    * rnd.normal() + pIn.get(i,0) );  // p


    }// for j - all sampling points
  }// for i - all dofs

  return res;

}





MATRIX ivr_FB_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd){
/**
  \brief Multi-dimensional  FB-MQC-IVR MonteCarlo: for general B = B(q,p)
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] TuningQ - a Ndof x Ndof matrix , cq,  see theory
  \param[in] TuningP - a Ndof x Ndof matrix , cp, see theory
  \param[in] flag - selection of special cases:
             0 - general,
             1 - B = B(q), do jump only in momentum
             2 - B = B(p), do jump only in position
  \param[in] Random - random number generator object
      
  Sample forward trajetory as well as position/momentum jump at time t

       w = |<p0,q0|pIn,qIn>|^2 exp(-cq Dq**2 /2) exp(-cp Dp**2 /2)

  Return value: a Ndof x 4 matrix with q0, p0, and jumps in position and momentum t Dq = qt - q0, Dp = pt - p0

*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_FB_MQC: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  MATRIX res(Ndof,4);  

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));
    double sq = 0.0;
    double sp = 0.0;
    if(flag==0 || flag==1){   sp = 1.0/sqrt(TuningP.get(i,i));   }
    if(flag==0 || flag==2){   sq = 1.0/sqrt(TuningQ.get(i,i));   }

    res.set(i, 0, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q0
    res.set(i, 1,    s    * rnd.normal() + pIn.get(i,0) );  // p0
    res.set(i, 2,   sq    * rnd.normal() );                 // Dq = qt' - qt
    res.set(i, 3,   sp    * rnd.normal() );                 // Dp = pt' - pt

  }// for i

  return res;

}


vector<MATRIX> ivr_FB_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd, int sample_size){
/**
  \brief Multi-dimensional  FB-MQC-IVR MonteCarlo: for general B = B(q,p)
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] TuningQ - a Ndof x Ndof matrix , cq,  see theory
  \param[in] TuningP - a Ndof x Ndof matrix , cp, see theory
  \param[in] flag - selection of special cases:
             0 - general,
             1 - B = B(q), do jump only in momentum
             2 - B = B(p), do jump only in position
  \param[in] Random - random number generator object
  \param[in] sample_size - how many points we want to sample
      
  Sample many points for forward trajetory as well as position/momentum jump at time t

       w = |<p0,q0|pIn,qIn>|^2 exp(-cq Dq**2 /2) exp(-cp Dp**2 /2)

  Return value: a vector of Ndof x 4 matrices with q0, p0, and jumps in position and momentum t Dq = qt - q0, Dp = pt - p0

*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_FB_MQC: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  vector<MATRIX> res(sample_size, MATRIX(Ndof,4));  

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));
    double sq = 0.0;
    double sp = 0.0;
    if(flag==0 || flag==1){   sp = 1.0/sqrt(TuningP.get(i,i));   }
    if(flag==0 || flag==2){   sq = 1.0/sqrt(TuningQ.get(i,i));   }

    for(int j=0; j<sample_size; j++){

      res[j].set(i, 0, (1.0/s) * rnd.normal() + qIn.get(i,0) );  // q0
      res[j].set(i, 1,    s    * rnd.normal() + pIn.get(i,0) );  // p0
      res[j].set(i, 2,   sq    * rnd.normal() );                 // Dq = qt' - qt
      res[j].set(i, 3,   sp    * rnd.normal() );                 // Dp = pt' - pt

    }// for j - all sampling points
  }// for i - all DOFs

  return res;

}



MATRIX ivr_FF_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd){
/**
  \brief Multi-dimensional  FF-MQC-IVR MonteCarlo: for general B
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] TuningQ - a Ndof x Ndof matrix , cq,  see theory
  \param[in] TuningP - a Ndof x Ndof matrix , cp, see theory
  \param[in] Random - random number generator object


  Sample initial average and difference variables for pair of forward trajectories

       w = |<pb,qb|pIn,qIn>|^2 exp(-cq Dq0**2 /2) exp(-cp Dp0**2 /2)

  Use average and difference to generated correlated initial conditions
     
  Return value: a Ndof x 4 matrix with q0, p0, and jumps in position and momentum t Dq = qt - q0, Dp = pt - p0

*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_FF_MQC: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  MATRIX res(Ndof,4);  

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));
    double sq = 1.0/sqrt(TuningQ.get(i,i));
    double sp = 1.0/sqrt(TuningP.get(i,i));

    double qav = (1.0/s) * rnd.normal() + qIn.get(i,0);
    double pav = s * rnd.normal() + pIn.get(i,0);
    double Dq0 = sq    * rnd.normal();
    double Dp0 = sp    * rnd.normal();

    res.set(i, 0,  qav - 0.5 * Dq0);  // q0
    res.set(i, 1,  pav - 0.5 * Dp0);  // p0
    res.set(i, 2,  qav + 0.5 * Dq0);  // q0'
    res.set(i, 3,  pav + 0.5 * Dp0);  // p0'


  }// for i

  return res;

}


vector<MATRIX> ivr_FF_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd, int sample_size){
/**
  \brief Multi-dimensional  FF-MQC-IVR MonteCarlo: for general B
  \param[in] qIn - initial coordinates (Ndof x 1 matrix)
  \param[in] pIn - initial momenta (Ndof x 1 matrix)
  \param[in] Width0 - a Ndof x Ndof matrix with the width parameters of coherent states (CS) 
  \param[in] TuningQ - a Ndof x Ndof matrix , cq,  see theory
  \param[in] TuningP - a Ndof x Ndof matrix , cp, see theory
  \param[in] Random - random number generator object
  \param[in] sample_size - how many points we want to sample

  Sample many points for initial average and difference variables for pair of forward trajectories

       w = |<pb,qb|pIn,qIn>|^2 exp(-cq Dq0**2 /2) exp(-cp Dp0**2 /2)

  Use average and difference to generated correlated initial conditions
     
  Return value: a vector of Ndof x 4 matrices with q0, p0, and jumps in position and momentum t Dq = qt - q0, Dp = pt - p0

*/

  if(qIn.n_rows!=qIn.n_rows){
    cout<<"Error in ivr_FF_MQC: the input matrices qIn and pIn have different # of rows\n";
    exit(0);
  }
  
  int Ndof = qIn.n_rows;

  vector<MATRIX> res(sample_size, MATRIX(Ndof,4));  

  for(int i=0; i<Ndof; i++){

    double s = sqrt(Width0.get(i,i));
    double sq = 1.0/sqrt(TuningQ.get(i,i));
    double sp = 1.0/sqrt(TuningP.get(i,i));

    for(int j=0;j<sample_size;j++){

      double qav = (1.0/s) * rnd.normal() + qIn.get(i,0);
      double pav = s * rnd.normal() + pIn.get(i,0);
      double Dq0 = sq    * rnd.normal();
      double Dp0 = sp    * rnd.normal();

      res[j].set(i, 0,  qav - 0.5 * Dq0);  // q0
      res[j].set(i, 1,  pav - 0.5 * Dp0);  // p0
      res[j].set(i, 2,  qav + 0.5 * Dq0);  // q0'
      res[j].set(i, 3,  pav + 0.5 * Dp0);  // p0'

    }// for j - all samples
  }// for i - all DOFs

  return res;

}


}// namespace libivr
}// liblibra

