#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "../mmath/libmmath.h"
#include "../cell/libcell.h"

using namespace libmmath;
using namespace libcell;

namespace libpot{

//------------------ Switching functions --------------------------
void SWITCH(VECTOR& r1,VECTOR&r2,
            double R_on,double R_off,
            double& SW,VECTOR& dSW);

void DOUBLE_SWITCH(double x,double a,double eps,double& SW,double& dSW);

//------------------ Bond potentials ------------------------------

double Bond_Harmonic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                     VECTOR& fi,VECTOR& fj,  /*Outputs*/
                      double K, double r0);  /*Parameters*/
double Bond_Quartic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                    VECTOR& fi,VECTOR& fj,  /*Outputs*/
                    double K, double r0);   /*Parameters*/
double Bond_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                  VECTOR& fi,VECTOR& fj,            /*Outputs*/
                  double D, double r0,double alp);  /*Parameters*/


//-------------------- Angle potentials ---------------------------

double Angle_Harmonic(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                      VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                      double k_theta,double theta_0     /* Parameters*/
                     );
double Angle_Fourier(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                     VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                     double k_ijk,double C0,double C1,
                     double C2,int coordination);

double Angle_Fourier_General(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                             double k_ijk,double C0,double C1,
                             double C2                        /* Parameters*/
                     );
double Angle_Fourier_Special(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                             double k_ijk,int coordination     /* Parameters*/
                            );
double Angle_Harmonic_Cos(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                          VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                          double k_theta,double cos_theta_0,int coordination /* Parameters*/
                         );
double Angle_Harmonic_Cos_General(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                          VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                          double k_theta,double cos_theta_0 /* Parameters*/
                         );
double Angle_Cubic(VECTOR& r1,VECTOR& r2,VECTOR& r3, /* Inputs */
                   VECTOR& f1,VECTOR& f2,VECTOR& f3, /* Outputs*/
                   double k_theta,double theta_0     /* Parameters*/
                   );

//------------------- Stretch-bend potentials --------------------------------
double Stretch_Bend_Harmonic(VECTOR& r1,VECTOR& r2,VECTOR& r3,         /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3,         /* Outputs*/
                             double k_ijk,double k_kji, double theta_0,
                             double r_ij0,double r_kj0                 /* Parameters*/
                            );


//------------------ Dihedral/Torsion potentials ------------------------------

double Dihedral_General(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl, /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl, /*Outputs*/
                        double Vphi,double phi0,int n,int opt        /*Parameters*/
                        );
double Dihedral_Fourier(VECTOR& ri,VECTOR& rj,VECTOR& rk,VECTOR& rl,    /*Inputs*/
                        VECTOR& fi,VECTOR& fj,VECTOR& fk,VECTOR& fl,    /*Outputs*/
                        double Vphi1,double Vphi2,double Vphi3,int opt  /*Parameters*/
                        );

//------------------ Out-of-plane (OOP) potentials -------------------------------

double OOP_Fourier(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,        /*Inputs*/
                   VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,        /*Outputs*/
                   double Kijkl,double C0,double C1,double C2,int opt  /*Parameters*/
                   );

double OOP_Wilson(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,   /*Inputs*/
                  VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,   /*Outputs*/
                  double Kijkl,double xi_0                       /*Parameters*/
                   );

double OOP_Harmonic(VECTOR& r0,VECTOR& r1,VECTOR& r2,VECTOR& r3,  /*Inputs*/
                    VECTOR& f0,VECTOR& f1,VECTOR& f2,VECTOR& f3,  /*Outputs*/
                    double Kijkl                                  /*Parameters*/
                   );

//--------------------------- Vdw potentials -------------------------------------------

double Vdw_LJ(VECTOR& ri,VECTOR& rj,          /*Inputs*/
              VECTOR& fi,VECTOR& fj,          /*Outputs*/
              double sigma, double espilon);  /*Parameters*/

double Vdw_Buffered14_7(VECTOR& ri,VECTOR& rj,          /*Inputs*/
                        VECTOR& fi,VECTOR& fj,          /*Outputs*/
                        double sigma, double espilon);  /*Parameters*/

double Vdw_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                 VECTOR& fi,VECTOR& fj,            /*Outputs*/
                 double D, double r0,double alp);  /*Parameters*/


//------------------------- Electrostatic potentials ----------------------------------------

double Elec_Coulomb(VECTOR& ri,VECTOR& rj,     /*Inputs*/
                    VECTOR& fi,VECTOR& fj,     /*Outputs*/
                    double qi,double qj,
                    double eps,double delta);  /*Parameters*/


//------------------------- Fragment-Fragment potentials ------------------------------------
double Gay_Berne(VECTOR& ri,VECTOR& rj,VECTOR& ui,VECTOR& uj,          /*Inputs*/
                 VECTOR& fi,VECTOR& fj,VECTOR& ti,VECTOR& tj,          /*Outputs*/
                 double di, double dj,double li,double lj,
                 double e0,double rat,double dw,double mu,double nu);  /*Parameters*/

double Girifalco12_6(VECTOR& ri,VECTOR& rj,
                     VECTOR& fi,VECTOR& fj,
                     double a,double alp,double bet
                    );

//--------------------- Many-body potentials ----------------------------------
double Elec_Ewald3D(VECTOR* r,         /* Inputs */
                    VECTOR* g,
                    VECTOR* m,
                    VECTOR* f,
                    MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs */
                    int sz,double* q,
                    int nexcl, int* excl1, int* excl2, double* scale,
                    MATRIX3x3* box,int rec_deg,int pbc_deg,
                    double etha,int is_cutoff, double R_on, double R_off,
                    int& time,vector< vector<triple> >& images, vector<triple>& central_translation,
                    double* dr2,double dT, int& is_update);  /* Parameters */                  

double Vdw_LJ(VECTOR* r,                                               /* Inputs */
              VECTOR* g,
              VECTOR* m,
              VECTOR* f,
              MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
              int sz,double* epsilon, double* sigma,
              int nexcl, int* excl1, int* excl2, double* scale,
              MATRIX3x3* box,int rec_deg,int pbc_deg,
              double etha,int is_cutoff, double R_on, double R_off,
              int& time, vector< vector<triple> >& images, vector<triple>& central_translation,
              double* dr2, double dT, int& is_update     /* Parameters */
             );

double Vdw_LJ1(VECTOR* r,                                               /* Inputs */
               VECTOR* g,
               VECTOR* m,
               VECTOR* f,
               MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
               int sz,double* epsilon, double* sigma,
               int nexcl, int* excl1, int* excl2, double* scale,
               MATRIX3x3* box,int rec_deg,int pbc_deg,
               double etha,int is_cutoff, double R_on, double R_off,
               int& time, vector< vector<quartet> >& images, vector<triple>& central_translation,
               double* dr2, double dT, int& is_update     /* Parameters */
             );

double Vdw_LJ2_no_excl(VECTOR* r,                                               /* Inputs */
                       VECTOR* g,
                       VECTOR* m,
                       VECTOR* f,
                       MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                       int sz,double* epsilon, double* sigma,
                       int nexcl, int* excl1, int* excl2, double* scale,
                       MATRIX3x3* box,int rec_deg,int pbc_deg,
                       double etha,int is_cutoff, double R_on, double R_off,
                       int& time,vector< vector<excl_scale> >& excl_scales
                      );

double Vdw_LJ2_excl(VECTOR* r,                                               /* Inputs */
                    VECTOR* g,
                    VECTOR* m,
                    VECTOR* f,
                    MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                    int sz,double* epsilon, double* sigma,
                    int nexcl, int* excl1, int* excl2, double* scale,
                    MATRIX3x3* box,int rec_deg,int pbc_deg,
                    double etha,int is_cutoff, double R_on, double R_off,
                    int& time,vector< vector<excl_scale> >& excl_scales
                   );


double LJ_Coulomb(VECTOR* r,                                               /* Inputs */
                  VECTOR* g,
                  VECTOR* m,
                  VECTOR* f,
                  MATRIX3x3& at_stress, MATRIX3x3& fr_stress, MATRIX3x3& ml_stress, /* Outputs*/
                  int sz,double* epsilon, double* sigma,double* charge,int is_cutoff, double R_on, double R_off,
                  int nexcl, int* excl1, int* excl2, double* scale       /* Parameters */
                 );

} // namespace libpot



#endif // POTENTIALS_H
