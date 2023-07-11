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

#ifndef POTENTIALS_ANGLES_H
#define POTENTIALS_ANGLES_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{

//-------------------- Angle potentials ---------------------------

double Angle_Harmonic(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                      VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                      double k_theta,double theta_0     /* Parameters*/
                     );
boost::python::list Angle_Harmonic(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                   double k_theta,double theta_0  /* Parameters*/
                                  );




double Angle_Fourier(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                     VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                     double k_ijk,double C0,double C1,
                     double C2,int coordination);
boost::python::list Angle_Fourier(VECTOR r1,VECTOR r2,VECTOR r3, double k_ijk,double C0,double C1, double C2,int coordination); 

double Angle_Fourier_General(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                             double k_ijk,double C0,double C1,
                             double C2                        /* Parameters*/
                     );
boost::python::list Angle_Fourier_General(VECTOR r1,VECTOR r2,VECTOR r3, double k_ijk,double C0,double C1, double C2); 

double Angle_Fourier_Special(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                             double k_ijk,int coordination     /* Parameters*/
                            );
boost::python::list Angle_Fourier_Special(VECTOR r1,VECTOR r2,VECTOR r3, double k_ijk,int coordination); 





double Angle_Harmonic_Cos(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                          VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                          double k_theta,double cos_theta_0,int coordination /* Parameters*/
                         );
boost::python::list Angle_Harmonic_Cos(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                   double k_theta,double cos_theta_0,int coordination  /* Parameters*/
                                  );

double Angle_Harmonic_Cos_General(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                          VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                          double k_theta,double cos_theta_0 /* Parameters*/
                         );
boost::python::list Angle_Harmonic_Cos_General(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                   double k_theta,double cos_theta_0  /* Parameters*/
                                  );




double Angle_Cubic(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                   VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                   double k_theta,double theta_0     /* Parameters*/
                   );
boost::python::list Angle_Cubic(VECTOR r1,VECTOR r2,VECTOR r3, /* Inputs */
                                double k_theta,double theta_0  /* Parameters*/
                               );


}// namespace libpot
}// liblibra

#endif // POTENTIALS_ANGLES_H

