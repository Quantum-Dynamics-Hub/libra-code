#ifndef ENSEMBLE_H
#define ENSEMBLE_H


#include "../../mmath/libmmath.h"
#include "../../hamiltonian/libhamiltonian.h"
#include "../nuclear/libnuclear.h"
#include "../electronic/libelectronic.h"

using namespace libmmath;
using namespace libmmath::librandom;
using namespace libhamiltonian;


namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;

namespace libensemble{


class Ensemble{
  // Object of Ensemble type essentially represents a time-dependent electron-nuclear wavefunction
  // el - represents electronic component
  // mol - represent nuclear component
  // This is a semiclassical representation of a wavefunction via a swarm (ensemble) of trajectories
  //
  // PSI(t) = sum_k {  delta((R,P)-mol[k]) * sum_i { el[k][i] * |i;k> } }
  //
  // mol[k] - k-th copy of the nuclear system (k-th trajectory)
  // el[k][i] - coefficients of i-th electronic basis state associated with k-th copy of the system
  // |i;k> - is i-th electronic basis state parameterized by nuclear degrees of freedom mol[k]
  //

  public:

  // Propagated variables
  int ntraj;                // number of trajectories in ensemble 
  int nnucl;                // number of nuclear degrees of freedom
  int nelec;                // number of electronic DOFs

  vector<int> is_active;    // flag stating if the i-th trajectory is active, if not - it is assumed to be adsorbed
  vector<Nuclear*>    mol;  // nuclear subsystems
  vector<Electronic*>  el;  // electronic subsystems
  vector<Hamiltonian*> ham; // Hamiltonian "handlers" - unique for each copy
  

  // Properties, countable, statistics
  // So far it is specific for 2D systems (and 1D too)

  vector<double> ave_q;    // averages of coordinates q
  vector<double> ave_p;    // averages of momenta p
  vector<double> sigma_q;  // spread in q
  vector<double> sigma_p;  // spread in p


  //----------- Methods -----------
  void init(int, int, int);

  Ensemble(){ ntraj = 0; }
  Ensemble(int,int,int);

  void se_pop(vector<double>&);
  void se_pop(vector<double>&,double,double);
  void sh_pop(vector<double>&);
  void sh_pop(vector<double>&,double,double);
  void sh_pop1(vector<double>&);
  void sh_pop1(vector<double>&,double,double);


//  void print_map(std::string prefix, double Xmin, double Xmax, double dx, double Ymin, double Ymax, double dy, int snap); // for 2D projections on XY plane
//  void integral_flux(vector< vector<double> >& Int_flx, double Xmin, double Xmax, double dx, double Ymin, double Ymax, double dy);

  void compute_averages();


};

}// namespace libensemble
}// namespace libdyn


#endif // ENSEMBLE
