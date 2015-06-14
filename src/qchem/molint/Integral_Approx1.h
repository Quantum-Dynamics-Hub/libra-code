#ifndef INTEGRAL_APPROX1_H
#define INTEGRAL_APPROX1_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{



//============================================================
//=== auxiliary functions according to Rosen =================
double A_plus_1(int n,double alpha);
double A_minus_1(int n,double alpha);
double B(int n,double alpha);
double D(int m,int n,int p);
double K2ab(int m,int n,double alpha,double betha,double R);
double K2aa(int m,int n,double alpha,double betha,double R);

//============================================================
//=== According to O. Kitao and T. Ogawa =====================
double Jab(int m,int n,double alpha,double betha,double R);
double Jab(int m,int n,double alpha,double betha);
double dJab_dqa(int m,int n,double alpha,double betha,double R);

//=============================================================
//==== Parameterized way to calculate Coulomb integral ========
double Coulomb_Integral(double R,int n_i, double Jii, double ksi_i, std::string type_i, double q_i,
                                 int n_j, double Jjj, double ksi_j, std::string type_j, double q_j,
                                 double epsilon, int mode);


}// namespace libmolint
}// namespace libqchem


#endif // INTEGRAL_APPROX1_H

