#include "Ensemble.h"


namespace libdyn{
namespace libensemble{


void Ensemble::init(int _ntraj, int _nelec, int _nnucl){
// Allocate memory for an ensemble of ntraj trajectories, with all electronic components 
// represented in a basis of nstates electronic state and with nuclear component having the
// dimensionality of nnucl

  int i;
  ntraj = _ntraj;
  nelec = _nelec;
  nnucl = _nnucl;

  // Allocate electronic part
  el = vector<Electronic*>(ntraj);
  for(i=0;i<ntraj;i++){ el[i] = new Electronic(nelec,0); }  // all electronic states are the ground states

  // Allocate nuclear part
  mol = vector<Nuclear*>(ntraj);
  for(i=0;i<ntraj;i++){ mol[i] = new Nuclear(nnucl); }

  // Allocate Hamiltonian handlers
  ham = vector<Hamiltonian*>(ntraj);
  for(i=0;i<ntraj;i++){ ham[i] = new Hamiltonian(); } // just a generic one

  // Activate all trajectories
  is_active = vector<int>(ntraj,1);


  // Allocate memory for statistical data  
  ave_q = vector<double>(nnucl, 0.0);
  ave_p = vector<double>(nnucl, 0.0);
  sigma_q = vector<double>(nnucl, 0.0);
  sigma_p = vector<double>(nnucl, 0.0);


}

Ensemble::Ensemble(int _ntraj, int _nstates, int _nnucl){
// Constructor
   init(_ntraj,_nstates,_nnucl);
}


void Ensemble::se_pop(vector<double>& pops,double xmin, double xmax){

  int i,j;

  // In case array is of wrong size
  if(pops.size()!=nelec){
    pops.reserve(nelec);
    pops.resize(nelec,0.0);  
  }

  // pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  // of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  // and for all nuclear degrees of freedom
  
  for(i=0;i<nelec;i++){ 
    int ntr = 0;
    pops[i] = 0.0;

/*
    for(j=0;j<ntraj;j++){
      if(mol[j]->q[j]>xmin && mol[j]->q[j] < xmax){ // only those trajectories that are in given range

        double q = el[j]->q[i];
        double p = el[j]->p[i];
        pops[i] += (q*q + p*p);
        ntr += 1;

      }// if
    }// for j - all trajectories
*/
    if(ntr>0){    pops[i] /= (double)ntr; }  // normalize

  }// for i - all electronic states
  

}

void Ensemble::se_pop(vector<double>& pops){
  // Wavefunction population of all states, without regard to nuclear wavefunction localization
  // This is done by taking very large box
  // If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  se_pop(pops,-1000000.0,1000000.0);

}


void Ensemble::sh_pop(vector<double>& pops,double xmin, double xmax){

  int i,j;

  // In case array is of wrong size
  if(pops.size()!=nelec){
    pops.reserve(nelec);
    pops.resize(nelec,0.0);  
  }

  // pops[i] - will contain the fraction of the wavefunction (as given by the average weight)
  // of i-th electronic basis state that is enclosed by a box [xmin, xmax] in all dimensions (x,y,z)
  // and for all nuclear degrees of freedom
  
  for(i=0;i<el[0]->nstates;i++){  // all electronic states
    pops[i] = 0.0;

/*
    for(j=0;j<size;j++){
      if(mol[j]->R[0].x>xmin && mol[j]->R[0].x<xmax && el[j]->istate == i){ // only those trajectories that are in given range

        pops[i] += 1.0;

      }// if
    }// for j - all trajectories
*/
    pops[i] /= (double)ntraj;   // normalize

  }// for i - all electronic states
  

}

void Ensemble::sh_pop(vector<double>& pops){
  // Wavefunction population of all states, without regard to nuclear wavefunction localization
  // This is done by taking very large box
  // If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  sh_pop(pops,-1000000.0,1000000.0);

}


void Ensemble::sh_pop1(vector<double>& pops,double xmin, double xmax){
/*
// This is specially for Marcus spin-boson problem

  int i,j;

  // In case array is of wrong size
  if(pops.size()!=el[0]->nstates){

    pops.reserve(el[0]->nstates);
    pops.resize(el[0]->nstates,0.0);  

  }



  for(i=0;i<el[0]->nstates;i++){  // all electronic states
    pops[i] = 0.0;
  }

 
  for(j=0;j<size;j++){ // for all trajectories

    ham[j]->set_status(0);
    ham[j]->compute(mol[j]);

    double E0, E1, H0, H1, V;
    E0 = ham[j]->H(mol[j],0,0).real();
    E1 = ham[j]->H(mol[j],1,1).real();
    H0 = ham[j]->H(mol[j],0,0,0).real(); // diabatic
    H1 = ham[j]->H(mol[j],1,1,0).real();
    V  = ham[j]->H(mol[j],0,1,0).real();



    for(i=0;i<el[0]->nstates;i++){  // all electronic states

      if(mol[j]->R[0].x>xmin && mol[j]->R[0].x<xmax){ // only those trajectories that are in given range

        if(el[j]->istate == 0){ // we are in 0-th adiabatic state
        

          if(i==0){             // Probability to be on 0-th (left) diabat - f0
           
            pops[0] += V*V/((H0 - E0)*(H0 - E0) + V*V);  // f0^2

          }
          else if(i==1){             // Probability to be on 1-th (right) diabat - g0
           
            pops[1] += (H0 - E0)*(H0 - E0)/((H0 - E0)*(H0 - E0) + V*V);  // g0^2

          }


        }// 0-th adiabatic state
        else if(el[j]->istate == 1){

          if(i==0){             // Probability to be on 0-th (left) diabat - f1
         
            pops[0] += V*V/((H0 - E1)*(H0 - E1) + V*V);  // f1

          }
          else if(i==1){             // Probability to be on 1-th (right) diabat - g1

            pops[1] += (H0 - E1)*(H0 - E1)/((H0 - E1)*(H0 - E1) + V*V);  // g1

          }


        }// 1-th adiabatic state

      }// if min < x < max

    }// for i - all electronic states
  }// for j - all trajectories



  for(i=0;i<el[0]->nstates;i++){  // all electronic states

    pops[i] /= (double)size;   // normalize

  }
  
*/
}

void Ensemble::sh_pop1(vector<double>& pops){
  // Wavefunction population of all states, without regard to nuclear wavefunction localization
  // This is done by taking very large box
  // If you are out of this box - this is kinda strange (you need to reduce time of simulation, maybe)

  sh_pop1(pops,-100000000.0,100000000.0);

}


void Ensemble::print_map(std::string prefix, double Xmin, double Xmax, double dx, double Ymin, double Ymax, double dy, int snap){
/*
 // for 2D projections on XY plane

  std::string snaps,filename;
  stringstream ss(stringstream::in | stringstream::out);
  stringstream ss1(stringstream::in | stringstream::out);
  stringstream ss3(stringstream::in | stringstream::out);

  ss << snap;  ss >> snaps;
  
  int Nx = floor((Xmax-Xmin)/dx);
  int Ny = floor((Ymax-Ymin)/dy);

  filename = prefix+".frame"+snaps;

  ofstream out(filename.c_str(),ios::out);

  for(int nx=0;nx<Nx;nx++){
    double X = Xmin + nx*dx;

    for(int ny=0;ny<Ny;ny++){
      double Y = Ymin + ny*dy;

      // Compute density of trajectories: number of trajectories in the box [X,X+dx] x [Y,Y+dy]
      double dens = 0.0;
      double cnt = 0.0;
      for(int traj=0;traj<size;traj++){
        for(int n=0;n<mol[traj]->Nnucl;n++){

          if( (X<mol[traj]->R[n].x) && (mol[traj]->R[n].x<=(X+dx))  &&
              (Y<mol[traj]->R[n].y) && (mol[traj]->R[n].y<=(Y+dy))
            ){
             dens += 1.0;
           }
          cnt += 1.0;
        }// for n
      }// for traj
      dens /= cnt;

      out<<X<<"  "<<Y<<"  "<<dens<<endl;
    }// for ny
    out<<"\n";
  }// for nx

  out.close();

*/
}// print map


void Ensemble::integral_flux(vector< vector<double> >& Int_flx, double Xmin, double Xmax, double dx, 
                                                                double Ymin, double Ymax, double dy){
 // for 2D projections on XY plane
/*
  int Nx = Int_flx.size();
  int Ny = Int_flx[0].size();

  //
  double denom = 0.0;
  for(int traj=0;traj<size;traj++){
    for(int n=0;n<mol[traj]->Nnucl;n++){
      denom += 1.0;
    }
  }
  denom = 1.0/denom;



  for(int traj=0;traj<size;traj++){
    for(int n=0;n<mol[traj]->Nnucl;n++){

      int nx = floor((mol[traj]->R[n].x - Xmin)/dx);
      int ny = floor((mol[traj]->R[n].y - Ymin)/dy);

      if(nx>=0 && nx<Nx && ny>=0 && ny<Ny){
        Int_flx[nx][ny] += denom;
      }

    }
  }

*/

}// print map





void Ensemble::compute_averages(){

/*
  ave_x = 0.0; 
  ave_x2 = 0.0; 
  ave_y = 0.0; 
  ave_y2 = 0.0; 

  ave_px = 0.0; 
  ave_px2 = 0.0; 
  ave_py = 0.0; 
  ave_py2 = 0.0; 

  ave_xpx = 0.0; 
  ave_ypy = 0.0; 


  for(int traj=0;traj<size;traj++){

    ave_x  += mol[traj]->R[0].x;
    ave_x2 += (mol[traj]->R[0].x * mol[traj]->R[0].x);
    ave_y  += mol[traj]->R[0].y;
    ave_y2 += (mol[traj]->R[0].y * mol[traj]->R[0].y);

    ave_px  += mol[traj]->P[0].x;
    ave_px2 += (mol[traj]->P[0].x * mol[traj]->P[0].x);
    ave_py  += mol[traj]->P[0].y;
    ave_py2 += (mol[traj]->P[0].y * mol[traj]->P[0].y);

    ave_xpx = mol[traj]->R[0].x * mol[traj]->P[0].x; 
    ave_ypy = mol[traj]->R[0].y * mol[traj]->P[0].y; 

  }

  ave_x  /= (float(size));
  ave_x2 /= (float(size));
  ave_y  /= (float(size));
  ave_y2 /= (float(size));

  ave_px  /= (float(size));
  ave_px2 /= (float(size));
  ave_py  /= (float(size));
  ave_py2 /= (float(size));

  ave_xpx /= (float(size));
  ave_ypy /= (float(size));


  sx = sqrt(ave_x2 - ave_x*ave_x); 
  sy = sqrt(ave_y2 - ave_y*ave_y); 

  psx = (ave_xpx - ave_x * ave_px)/sx;
  psy = (ave_ypy - ave_y * ave_py)/sy;
*/

}


/*
double Ensemble::Epot(){
  
  double res = 0.0; 
  for(int traj=0;traj<size;traj++){
    double ep = ham[traj]->H(mol[traj],0,0).real();
    res += ep;
  } 
  res /= ((double)size);

  return res;

}

double Ensemble::Ekin(){

  double res = 0.0; 
  for(int traj=0;traj<size;traj++){
    double ek = compute_kinetic_energy(mol[traj]);
    res += ek;
  }
  res /= ((double)size);

  return res;
}

double Ensemble::Etot(){

  return (Ekin() + Epot());

}
*/


}// namespace libensemble
}// namespace libdyn

