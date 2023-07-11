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

#include "Potentials_angles.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{


double Angle_Harmonic(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                      VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                      double k_theta,double theta_0     /* Parameters*/
                     ){
//******************** double HARMONIC_BENDING **************************
//*                                                                     *
//*         E =        k_theta*(theta_ijk-theta_0)^2;                   *
//*                                                                     *
//***********************************************************************
  VECTOR r12 = r1 - r2;
  VECTOR r32 = r3 - r2;
  VECTOR f12,f32;
  double d12,d32,cos_theta,sin_theta,theta,diff;
  double energy;
  d12 = r12.length();
  d32 = r32.length();

  // Calculate angle
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta >  1.0) { cos_theta = 1.0;}
  else if (cos_theta < -1.0) { cos_theta = -1.0;}
  sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  theta = acos(cos_theta);

  // Energy
  diff = theta-theta_0;
  energy = k_theta * diff * diff;

  // Forces
  // Attention!!!
  // sin_theta = 0, when we have a linear molecules!!!
  // sign included into next expressions!!!
  diff *= (2.0*k_theta) / sin_theta;
  f12 = -(diff/d12)*(r12.unit()*cos_theta - r32.unit());
  f1 = f12;
  f32 = -(diff/d32)*(r32.unit()*cos_theta - r12.unit());
  f3 = f32;
  f2 = -f12 - f32;

  return energy;
}

double Angle_Fourier(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                     VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                     double k_ijk,double C0,double C1,
                     double C2,int coordination        /* Parameters*/
                     ){
  double energy = 0.0;
  if((coordination==1)||(coordination==3)||(coordination==4)){
    energy = Angle_Fourier_Special(r1,r2,r3,f1,f2,f3,k_ijk,coordination);
  }
  else{
    energy = Angle_Fourier_General(r1,r2,r3,f1,f2,f3,k_ijk,C0,C1,C2);
  }
  return energy;
}

double Angle_Fourier_General(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                             double k_ijk,double C0,double C1,
                             double C2                        /* Parameters*/
                     ){
//******************** double FOURIER_SERIES ****************************
//*                                                                     *
//*        u_theta = K_ijk*[C_0 + C_1*cos(theta) + C_2*cos(2*theta)]    *
//*                                                                     *
//*  Force Fields: UFF                                                  *
//*                                                                     *
//***********************************************************************
  VECTOR r12 = r1 - r2;
  VECTOR r32 = r3 - r2;
  VECTOR f12,f32;
  double d12,d32,cos_theta,sin_theta,theta,diff;
  double energy,cos_3theta,cos_4theta;
  d12 = r12.length();
  d32 = r32.length();

  // Calculate angle
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta >  1.0) { cos_theta = 1.0;}
  else if (cos_theta < -1.0) { cos_theta = -1.0;}
  sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  theta = acos(cos_theta);

  // Energy and forces
  //diff = -dE/dcos(theta)
  energy = k_ijk*(C0 + C1*cos_theta + C2*(2.0*cos_theta*cos_theta - 1.0));
  diff = k_ijk*(-C1 - 4.0*C2*cos_theta);

  f12 = -(diff/d12)*(r12.unit()*cos_theta - r32.unit());
  f1 = f12;
  f32 = -(diff/d32)*(r32.unit()*cos_theta - r12.unit());
  f3 = f32;
  f2 = -f12 - f32;

  return energy;
}

double Angle_Fourier_Special(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                             double k_ijk,int coordination     /* Parameters*/
                            ){
//******************** double FOURIER_SERIES ****************************
//*                                                                     *
//*   coordination == n == 1, 3, 4                                      *
//*                                                                     *
//*                u_theta =  (K_ijk/n^2)*[1 - cos(n*theta)]            *
//*                                                                     *
//*  Force Fields: UFF                                                  *
//*                                                                     *
//***********************************************************************
  VECTOR r12 = r1 - r2;
  VECTOR r32 = r3 - r2;
  VECTOR f12,f32;
  double d12,d32,cos_theta,sin_theta,theta,diff;
  double energy,cos_3theta,cos_4theta;
  d12 = r12.length();
  d32 = r32.length();

  // Calculate angle
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta >  1.0) { cos_theta = 1.0;}
  else if (cos_theta < -1.0) { cos_theta = -1.0;}
  sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  theta = acos(cos_theta);

//  cout<<endl<<"=============== Special angle ===========================\n";
//  cout<<"r1 = "<<r1<<" r2 = "<<r2<<" r3 = "<<r3<<endl;
//  cout<<"theta = "<<theta<<" cos_theta = "<<cos_theta<<" coordination = "<<coordination;

  // Energy and forces
  //diff = -dE/dcos(theta)
  switch(coordination){
    case 1: {
      energy = k_ijk*(1.0+cos_theta);  // Note! in UFF article there is "-" in front of cos(theta)
                                       // this is inapproriate. Look for example DREIDING atricle
                                       // the same case uses "+" sign
      diff = -k_ijk;
    } break;
    case 3: {
      // Keep original "-" sign for these cases!!!
      // According to formula: cos(3x)=4(cos(x))^3-3cos(x)
      cos_3theta = ((4.0*cos_theta*cos_theta)-3.0)*cos_theta;
      energy = (1.0/9.0)*k_ijk*(1.0-cos_3theta);
      diff = (1.0/9.0)*k_ijk*(12.0*cos_theta*cos_theta-3.0);
    } break;
    case 4: {
      // According to formula: cos(4x)=8*(cos(x))^4-8*(cos(x))^2+1
      cos_4theta = 8.0*(cos_theta*cos_theta)*(cos_theta*cos_theta-1.0) + 1.0;
      energy = (1.0/16.0)*k_ijk*(1-cos_4theta);
      diff = (1.0/16.0)*k_ijk*16.0*cos_theta*(2.0*cos_theta-1.0);
    } break;
  }// switch coordination

//  cout<<" energy = "<<energy<<endl;

  f12 = -(diff/d12)*(r12.unit()*cos_theta - r32.unit());
  f1 = f12;
  f32 = -(diff/d32)*(r32.unit()*cos_theta - r12.unit());
  f3 = f32;
  f2 = -f12 - f32;

  return energy;
}

double Angle_Harmonic_Cos(VECTOR& r1,VECTOR& r2,VECTOR& r3,                  /* Inputs */
                          VECTOR& f1,VECTOR& f2,VECTOR& f3,                  /* Outputs*/
                          double k_theta,double cos_theta_0,int coordination /* Parameters*/
                          ){
  double energy = 0.0;
  if(coordination==1){
    energy = Angle_Fourier_Special(r1,r2,r3,f1,f2,f3,k_theta,coordination);
  }
  else{
    k_theta = 0.5*k_theta/(1.0-cos_theta_0*cos_theta_0);
    energy = Angle_Harmonic_Cos_General(r1,r2,r3,f1,f2,f3,k_theta,cos_theta_0);
  }
  return energy;
}


double Angle_Harmonic_Cos_General(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                                  VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                                  double k_theta,double cos_theta_0 /* Parameters*/
                                  ){
//******************** double HARMONIC_COS_BENDING **********************
//*                                                                     *
//*         E =        k_theta*(cos(theta_ijk)-cos(theta_0))^2;         *
//*                                                                     *
//***********************************************************************
  VECTOR r12 = r1 - r2;
  VECTOR r32 = r3 - r2;
  VECTOR f12,f32;
  double d12,d32,cos_theta,sin_theta,theta,diff;
  double energy;
  d12 = r12.length();
  d32 = r32.length();

  // Calculate angle
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta >  1.0) { cos_theta = 1.0;}
  else if (cos_theta < -1.0) { cos_theta = -1.0;}
  sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  theta = acos(cos_theta);

  // Energy
  diff = cos_theta-cos_theta_0;
  energy = k_theta * diff * diff;

  // Forces
  // Attention!!!
  // sin_theta = 0, when we have a linear molecules!!!
  // sign included into next expressions!!!
  diff *= (2.0*k_theta);
  f12 = -(diff/d12)*(r12.unit()*cos_theta - r32.unit());
  f1 = f12;
  f32 = -(diff/d32)*(r32.unit()*cos_theta - r12.unit());
  f3 = f32;
  f2 = -f12 - f32;

  return energy;
}

double Angle_Cubic(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                   VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                   double k_theta,double theta_0     /* Parameters*/
                   ){
//******************** double HARMONIC_BENDING **************************
//*                                                                     *
//*         E =        k_theta*(theta_ijk-theta_0)^2;                   *
//*                                                                     *
//***********************************************************************
  double cb = -0.4;
  VECTOR r12 = r1 - r2;
  VECTOR r32 = r3 - r2;
  VECTOR f12,f32;
  double d12,d32,cos_theta,sin_theta,theta,diff;
  double energy;
  d12 = r12.length();
  d32 = r32.length();

  // Calculate angle
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta >  1.0) { cos_theta = 1.0;}
  else if (cos_theta < -1.0) { cos_theta = -1.0;}
  sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  theta = acos(cos_theta);

  // Energy
  diff = theta-theta_0;
  energy = k_theta * diff * diff*(1.0 + cb * diff);

  // Forces
  // Attention!!!
  // sin_theta = 0, when we have a linear molecules!!!
  // this is -dE/dcos(theta) = -dE/dtheta *(-1/sin(theta))
  diff *= (k_theta*(1.0 + 1.5*cb*diff)) / sin_theta;
  // - dE/dcos(theta) *(dcos(theta)/drx) x = i,j,k
  f12 = -(diff/d12)*(r12.unit()*cos_theta - r32.unit());
  f1 = f12;
  f32 = -(diff/d32)*(r32.unit()*cos_theta - r12.unit());
  f3 = f32;
  f2 = -f12 - f32;

  return energy;
}


//================================= Exported functions ==================
boost::python::list Angle_Harmonic(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                   double k_theta,double theta_0  /* Parameters*/
                                  ){
  boost::python::list res;
  double en = 0.0;
  VECTOR f1, f2, f3;

  en = Angle_Harmonic(r1,r2,r3,f1,f2,f3,k_theta,theta_0);

  res.append(en);
  res.append(f1);
  res.append(f2);
  res.append(f3);
 
  return res;
}


boost::python::list Angle_Fourier(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                  double k_ijk,double C0,double C1,double C2, int coordination
                                 ){
  boost::python::list res;
  double en = 0.0;
  VECTOR f1, f2, f3;

  en = Angle_Fourier(r1,r2,r3,f1,f2,f3,k_ijk,C0,C1,C2,coordination);

  res.append(en);
  res.append(f1);
  res.append(f2);
  res.append(f3);
 
  return res;
}


boost::python::list Angle_Fourier_General(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                  double k_ijk,double C0,double C1,double C2
                                 ){
  boost::python::list res;
  double en = 0.0;
  VECTOR f1, f2, f3;

  en = Angle_Fourier_General(r1,r2,r3,f1,f2,f3,k_ijk,C0,C1,C2);

  res.append(en);
  res.append(f1);
  res.append(f2);
  res.append(f3);
 
  return res;
}

boost::python::list Angle_Fourier_Special(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                  double k_ijk,int coordination
                                 ){
  boost::python::list res;
  double en = 0.0;
  VECTOR f1, f2, f3;

  en = Angle_Fourier_Special(r1,r2,r3,f1,f2,f3,k_ijk,coordination);

  res.append(en);
  res.append(f1);
  res.append(f2);
  res.append(f3);
 
  return res;
}


boost::python::list Angle_Harmonic_Cos(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                   double k_theta,double cos_theta_0,int coordination  /* Parameters*/
                                  ){
  boost::python::list res;
  double en = 0.0;
  VECTOR f1, f2, f3;

  en = Angle_Harmonic_Cos(r1,r2,r3,f1,f2,f3,k_theta,cos_theta_0,coordination);

  res.append(en);
  res.append(f1);
  res.append(f2);
  res.append(f3);
 
  return res;
}

boost::python::list Angle_Harmonic_Cos_General(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                   double k_theta,double cos_theta_0  /* Parameters*/
                                  ){
  boost::python::list res;
  double en = 0.0;
  VECTOR f1, f2, f3;

  en = Angle_Harmonic_Cos_General(r1,r2,r3,f1,f2,f3,k_theta,cos_theta_0);

  res.append(en);
  res.append(f1);
  res.append(f2);
  res.append(f3);
 
  return res;
}

boost::python::list Angle_Cubic(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                double k_theta,double theta_0  /* Parameters*/
                               ){
  boost::python::list res;
  double en = 0.0;
  VECTOR f1, f2, f3;

  en = Angle_Cubic(r1,r2,r3,f1,f2,f3,k_theta,theta_0);

  res.append(en);
  res.append(f1);
  res.append(f2);
  res.append(f3);
 
  return res;
}




} // namespace libpot
}// liblibra

