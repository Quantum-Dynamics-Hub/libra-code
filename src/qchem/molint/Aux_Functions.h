#ifndef AUXFUNCT_H
#define AUXFUNCT_H

#include "../../mmath/libmmath.h"
using namespace libmmath;


namespace libqchem{
namespace libmolint{

// Auxiliary functions
//----- For GTO quadratures ----------

/*double Aux_Integral2(int n1,int n2,double x1,double x2,double alp,double& dI_dx1,double& dI_dx2,vector<double*>& aux,int n_aux);

void Aux_Function4(int n1,int n2,double PA,double PB,double PC,double gamma,
                   double* G,double* dGdA, double* dGdB, double* dGdC,double* f, double* dfda, double* dfdb, int n_aux);

*/
void Aux_Function5(int ,int,double,double,double,    double*, double*, double*,   double*, double*, double*, int);

void Aux_Function6(int ,int,int,int,double,double,double,double,double, double,double,
 double* , double* , double* , double* , double* , double* ,  double* , double* , double* , 
 double* , double* , double* , double* , double* , double* ,
 int );





}// namespace libmolint
}// namespace libqchem


#endif // AUXFUNCT_H

