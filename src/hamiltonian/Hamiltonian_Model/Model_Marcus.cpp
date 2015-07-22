#include "Model_Marcus.h"

namespace libhamiltonian{
namespace libhamiltonian_model{


void Marcus_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params){ 
// Marcus spin-boson Hmiltonian in diabatic representation

  if(H->num_of_elems!=4){ std::cout<<"Error in Marcus_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH->num_of_elems!=4){ std::cout<<"Error in Marcus_Ham: dH matrix must be allocated\n"; exit(0);}
  if(d2H->num_of_elems!=4){ std::cout<<"Error in Marcus_Ham: d2H matrix must be allocated\n"; exit(0);}


  // Marcus spin-boson potetnial
  // Default parameters: from Landry-Subotnik paper (a.u.)
  double mass = 1.0;       double omega = 3.5e-4;
  double E_r = 2.39e-2;    double V = 5.0e-5; 
  double eps_0 = 2.0e-2;

  if(params.size()>=5){
    mass = params[0];    omega = params[1];
    E_r = params[2];     V = params[3];
    eps_0 = params[4];
  }


  // Derived (auxiliary) parameters
  double mo2 = 0.5*mass * omega * omega;
  double M = sqrt(mo2*E_r);

  // H00
  H->M[0] = mo2*x*x + M*x;
  dH->M[0] = 2.0*mo2*x + M;
  d2H->M[0] = 2.0*mo2;

  // H11
  H->M[3] = mo2*x*x - M*x - eps_0; 
  dH->M[3] = 2.0*mo2*x - M;
  d2H->M[3] = 2.0*mo2;

  // H01 and H10
  H->M[1] = H->M[2] = V;
  dH->M[1] = dH->M[2] = 0.0;
  d2H->M[1] = d2H->M[2] = 0.0;


}

boost::python::list Marcus_Ham(double x, boost::python::list params_){ 
// Marcus spin-boson hamiltonian in diabatic representation

  MATRIX H(2,2);
  MATRIX dH(2,2);
  MATRIX d2H(2,2);

  int sz = boost::python::len(params_);
  vector<double> params(sz,0.0);
  for(int i=0;i<sz;i++){ params[i] = boost::python::extract<double>(params_[i]);  }

  Marcus_Ham(x,&H,&dH,&d2H,params);

  boost::python::list res;
  res.append(x);
  res.append(H);
  res.append(dH);
  res.append(d2H);
 
  return res;

}

}// namespace libhamiltonian_model
}// namespace libhamiltonian

