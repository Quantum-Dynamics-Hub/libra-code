#include "Dynamics_Ensemble.h"

namespace libdyn{


void propagate_ensemble(double dt,Ensemble* ens,int opt){

  // Ensemble of independent trajectories
  for(int i=0;i<ens->ntraj;i++){

    if(ens->is_active[i]){

      propagate_electronic(0.5*dt,&ens->el[i], ens->ham[i]);
      propagate_nuclear(dt, &ens->mol[i], &ens->el[i], ens->ham[i],opt);
      propagate_electronic(0.5*dt,&ens->el[i], ens->ham[i]);

    }
       
  }// for i

}

void propagate_ensemble(double dt,Ensemble& ens,int opt){

  propagate_ensemble(dt, &ens, opt);
}






}// namespace libdyn


