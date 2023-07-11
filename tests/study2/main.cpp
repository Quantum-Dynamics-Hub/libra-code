
#include "libdyn.h"
#include "libmmath.h"

using namespace libdyn;

void test1(){

  // This test computes diabatic or adiabatic energies of NACs along the coordinate

  // This is a convenient option to set customary Hamiltonian
  // The base "Hamiltonian" class simply provides an interface to other functions (dynamics)
  Nuclear* mol; mol = new Nuclear(1);  
  Hamiltonian* ham;
  ham = new Hamiltonian_Tully(1, 0, mol);


  for(double x=-5;x<=5.0;x+=0.1){
    mol->R[0].x = x;
    mol->P[0].x = 1.0;
//    ham->update_nuclear(mol);
    ham->compute(mol);
    std::cout<<x<<"  "<<ham->H(mol,0,0).real()<<" "<<ham->H(mol,1,1).real()<<"  "<<ham->H(mol,0,1).imag()<<"\n";
//    std::cout<<x<<"  "<<ham->dHdRx(0,0,0).real()<<" "<<ham->dHdRx(1,1,0).real()<<"  "<<ham->dHdRx(0,1,0).imag()<<"\n";
  }

}

void test3(){

  // This test computes diabatic or adiabatic energies of NACs along the coordinate

  // This is a convenient option to set customary Hamiltonian
  // The base "Hamiltonian" class simply provides an interface to other functions (dynamics)
  Nuclear* mol; mol = new Nuclear(1);  
  mol->mass[0] = 1.0;
  Hamiltonian* ham;
  int rep = 1; // 0 - diabatic, 1 - adiabatic
  int Ham_indx = 3; // potential = Marcus

  ham = new Hamiltonian_Tully(rep, Ham_indx, mol);

  ham->set_param("C",5.0e-5);

  double dE = (3.0e-2 - 1.5e-2)/16.0;
  ham->set_param("D",1.5e-2+0*dE);



  ofstream out("marcus_profile.txt",ios::out);

  for(double x=-500;x<=250.0;x+=1.0){
    mol->R[0].x = x;
    mol->P[0].x = 1.0;
    ham->set_status(0);    // because the geometry has changed
    ham->compute(mol);

    if(rep==0){
      std::cout<<x<<"  "<<ham->H(mol,0,0).real()<<" "<<ham->H(mol,1,1).real()<<"  "<<ham->H(mol,0,1).real()<<"\n";
//    std::cout<<x<<"  "<<ham->dHdRx(0,0,0).real()<<" "<<ham->dHdRx(1,1,0).real()<<"  "<<ham->dHdRx(0,1,0).imag()<<"\n";
    }
    else if(rep==1){
      std::cout<<x<<"  "<<ham->H(mol,0,0).real()<<" "<<ham->H(mol,1,1).real()<<"  "<<ham->H(mol,0,1).imag()<<"\n";
//    std::cout<<x<<"  "<<ham->dHdRx(0,0,0).real()<<" "<<ham->dHdRx(1,1,0).real()<<"  "<<ham->dHdRx(0,1,0).imag()<<"\n";
    }



        double E0, E1, H0, H1, V;
        E0 = ham->H(mol,0,0).real();
        E1 = ham->H(mol,1,1).real();
        H0 = ham->H(mol,0,0,0).real(); // diabatic
        H1 = ham->H(mol,1,1,0).real();
        V  = ham->H(mol,0,1,0).real();

        double f0 = V*V/((H0 - E0)*(H0 - E0) + V*V);  // f0
        double g0 = (H0 - E0)*(H0 - E0)/((H0 - E0)*(H0 - E0) + V*V); // g0           
        double f1 = V*V/((H0 - E1)*(H0 - E1) + V*V);  // f1
        double g1 = (H0 - E1)*(H0 - E1)/((H0 - E1)*(H0 - E1) + V*V); // g1


        out<<x<<"  "<<E0<<"  "<<E1<<"  "<<H0<<"  "<<H1<<"  "<<V<<"  "<<f0<<"  "<<g0<<"  "<<f1<<"  "<<g1<<endl;


  }// for x

  out.close();

}



void run_semiclassical(){

// Tully models 

  // Init random seed
//  srand(time(0));

  // Simulation parameters
  int nstates = 2;       // number of electronic states
  int Nnucl = 1;         // number of nuclei
  int nsnap = 500;       // number of snapshots
  int Ham_indx = 2;      // SAC = 0, DAC = 1, ECWR = 2
  int nstep = 10;        // number of steps per 1 snap;
  int ntraj = 250; 
  int do_rescaling = 1;  // 0 - CPA, 1 - noraml electron back-reaction
  int method = 2;        // 0 - FSSH,  1 - MSSH, 2 - FSSH + FC

  double xstart = -15.0; // initial nuclear coordinate
  int isurface = 0;      // initial electronic state
  int rep = 1;           // 0 - diabatic, 1 = adiabatic

  int i,j,k;

  vector< vector<double> > Pij(nstates,vector<double>(nstates,0.0));


//  Hamiltonian* ham;
//  ham = new Hamiltonian_Tully(1, 0, -5.0, 1.0);

//  Nuclear* mol; mol = new Nuclear(1);  

  
  ofstream f("collect.txt",ios::out);
  f.close();



  //############# Loop over initial k #############
  for(double kstart= 1; kstart<=50; kstart+=0.25){


    double dt = 0.1*(2000.0/kstart);  //prms.dt; a.u.

    // Simulation variables  
    Ensemble* ens; ens = new Ensemble(ntraj, nstates, Nnucl);

    vector< vector< vector<double> > > FCij(ntraj,vector< vector<double> >(nstates, vector<double>(nstates, 1.0) ));

//

//    for(i=0;i<prms.num_sh_traj;i++){   
//      cout<<"i= "<<i<<" q[0]= "<<ens.el[0]->q[0]<<" p[0]= "<<ens.el[0]->p[0]<<" q[1]= "<<ens.el[0]->q[1]<<" p[1]= "<<ens.el[0]->p[1]<<endl;
//    }


    // Initialization
    for(i=0;i<ntraj;i++){

        // Init nuclear variables - beyond default constructor
      ens->mol[i]->mass[0] = 2000.0; //prms.m0;
      ens->mol[i]->R[0].x = xstart + 0.0*normal();  
      ens->mol[i]->R[0].y = 0.0;
      ens->mol[i]->R[0].z = 0.0;

      ens->mol[i]->P[0].x = kstart + 0.0*normal();
      ens->mol[i]->P[0].y = 0.0;
      ens->mol[i]->P[0].z = 0.0;

//      cout<<"kstart = "<<ens->mol[i]->P[0].x<<endl;
      ens->mol[i]->F[0] = 0.0;

      ens->el[i]->istate = isurface;

      ens->ham[i] = new Hamiltonian_Tully(rep, Ham_indx, ens->mol[i]);



    }// i

//    exit(0);


    // Observables - reflection and transition probabilities on all states
    vector<double> R(nstates,0.0);
    vector<double> T(nstates,0.0);
                     
    //################# Loop over snapshots ######################
    for(int isnap=0; isnap<nsnap; isnap++){
      for(int istep=0; istep<nstep; istep++){

        propagate_ensemble(dt, ens, "FSSH");

        //---------- Now do surface hopping stuff -------------
        // Manually
        for(int traj=0;traj<ens->size;traj++){

          if(ens->is_active[traj]){

            if(method==0){
              compute_hopping_probabilities(ens->mol[traj], ens->el[traj], ens->ham[traj], Pij, dt, 0, 300.0 );
            }
            else if(method==1){
              compute_hopping_probabilities_mssh(ens->mol[traj], ens->el[traj], ens->ham[traj], Pij, dt, 0, 300.0 );
            }
            if(method==2){
              compute_hopping_probabilities_decoh(ens->mol[traj], ens->el[traj], ens->ham[traj], Pij, FCij[traj], dt, 0, 300.0 );
            }


            double ksi = uniform(0.0,1.0);

/*
          if(ens->mol[traj]->R[0].x<0.1 && ens->mol[traj]->R[0].x>-0.1){

          cout<<"traj= "<<traj<<" ksi= "<<ksi<<" istate= "<<ens->el[traj]->istate<<endl;
          cout<<"x = "<<ens->mol[traj]->R[0].x<<endl;
          for(i=0;i<Pij.size();i++){
            for(j=0;j<Pij[i].size();j++){
              cout<<Pij[i][j]<<" ";
            }
            cout<<endl;
          }
          }
*/
 
          if(method==0||method==1){

            hop(ens->el[traj]->istate, ens->mol[traj], ens->ham[traj], ksi, Pij, do_rescaling, rep);

          }
          else if(method==2){

            hop_decoh(ens->el[traj]->istate, ens->mol[traj], ens->ham[traj], ksi, Pij, FCij[traj], do_rescaling, rep);

          }

//          cout<<" istate(after)= "<<ens->el[traj]->istate<<endl;

            // Deactivate trajectories:
            if(ens->mol[traj]->R[0].x <-20 || ens->mol[traj]->R[0].x > 10.0) { ens->is_active[traj] = 0; }

          }// if active
         
        }// for traj


      }// for istep
     
 
      double Epot = compute_potential_energy(ens->mol[0], ens->el[0], ens->ham[0], "FSSH");
      double Ekin = compute_kinetic_energy(ens->mol[0]);

//      cout<<"px = "<<kstart<<" isnap= "<<isnap<<" F= "<<ens->mol[0]->F[0]<<"  R= "<<ens->mol[0]->R[0]<<" P= "<<ens->mol[0]->P[0]<<endl;
/*
      cout<<"px = "<<kstart<<" isnap= "<<isnap<<" x= "<<ens->mol[0]->R[0].x
          <<" Ekin= "<<Ekin<<" Epot= "<<Epot<<" Etot= "<<Ekin+Epot
          <<" istate= "<<ens->el[0]->istate
          <<" Fx= "<<ens->mol[0]->F[0].x<<" P01= "<<Pij[0][1]
          <<" q[0]= "<<ens->el[0]->q[0]<<" p[0]= "<<ens->el[0]->p[0]
          <<" q[1]= "<<ens->el[0]->q[1]<<" p[1]= "<<ens->el[0]->p[1]<<endl;
*/


    }// for isnap

    ens->sh_pop(R,-1000000.0, 0.0);
    ens->sh_pop(T, 0.0, 1000000.0);



    // Print probabilities 
    ofstream f("collect.txt",ios::app);    
    f<<"kstart = "<<kstart<<"  ";


    for(i=0;i<nstates;i++){
      f<<" R["<<i<<"]= "<<R[i];
      f<<" T["<<i<<"]= "<<T[i];
      f<<" Err["<<i<<"]= "<<0.0;

    }
    f<<" R= "<<ens->mol[0]->R[0].x;
    f<<endl;
    f.close();



    delete ens;


  }// for kstart



}// run_semiclassical



void boltz(double& x_,double& p_){


  double Er = 2.39e-2;
  double omega = 3.5e-4;
  double kT = 9.5e-4;

 
  double mo2 = 0.5*omega*omega; // mass = 1
  double M = sqrt(mo2*Er);

  x_ = -M/mo2;           // minimum
  p_ = sqrt(kT);         // momentum that corresponds to temperatures

  double Eold = mo2*x_*x_ + M*x_ + 0.5*p_*p_;    // energy
  double Enew = 0.0;


  for(int i=0;i<1000;i++){
    double x = 3.0*(M/mo2)*normal();  // proposed state x
    double p = 3.0*sqrt(kT)*normal(); // proposed state p

    Enew = mo2*x*x + M*x + 0.5*p*p;

    double dE = Enew - Eold;

    double ksi = uniform(0.0,1.0);
    double prob = min(1.0, exp(-dE/kT)); 

    if(ksi<prob){  // accept new state with Metropolis scheme
      Eold = Enew; 
      x_ = x;
      p_ = p;
    }
   
  }// for i



}

void run_Marcus(){
// Marcus spin-boson model

  // Init random seed
  srand(time(0));


  // Simulation parameters
  int nstates = 2;       // number of electronic states
  int Nnucl = 1;         // number of nuclei
  int nsnap = 1000;      // number of snapshots
  int Ham_indx = 3;      // SAC = 0, DAC = 1, ECWR = 2, Marcus = 3
  int nstep = 200;       // number of steps per 1 snap;
  int ntraj = 1000; 
  int do_rescaling = 1;  // 0 - CPA, 1 - noraml electron back-reaction
  int method = 0;        // 0 - FSSH,  1 - MSSH, 2 - FSSH + FC


//  double xstart = -1.0; // initial nuclear coordinate
  int isurface = 0;      // initial electronic state
  int rep = 1;           // 0 - diabatic, 1 = adiabatic
  double dt = 2.5;       // a.u.

  int i,j,k;

  vector< vector<double> > Pij(nstates,vector<double>(nstates,0.0));



  ofstream f("collect.txt",ios::out);
  f.close();


  // Model parameters
  double Er = 2.39e-2;
  double omega = 3.5e-4;
  double kT = 9.5e-4;


  const double kb = 3.166811429e-6; // Hartree/K 

  double T = kT / kb;  // ~300 K
  


  //############# ??? #############
  for(int ie=0;ie<=16;ie++){  // fake parameter loop

    cout<<"ie= "<<ie<<endl;

    // Simulation variables  
    Ensemble* ens; ens = new Ensemble(ntraj, nstates, Nnucl);

    vector< vector< vector<double> > > FCij(ntraj,vector< vector<double> >(nstates, vector<double>(nstates, 1.0) ));

    Thermostat* therm;  therm = new Thermostat();

// NVT
    therm->set_gamma(3.0e-4);
    therm->set_sigma(sqrt(2.0*1.0*3.0e-4*kT/dt));

// NVE
//    therm->set_gamma(0.0);                    
//    therm->set_sigma(0.0);

    // Initialization
    for(i=0;i<ntraj;i++){

//      cout<<"trajectory# "<<i<<endl;

        // Init nuclear variables - beyond default constructor
      ens->mol[i]->mass[0] = 1.0;
      
      boltz(ens->mol[i]->R[0].x, ens->mol[i]->P[0].x);

      ens->mol[i]->R[0].y = 0.0;
      ens->mol[i]->R[0].z = 0.0;

      ens->mol[i]->P[0].y = 0.0;
      ens->mol[i]->P[0].z = 0.0;

      ens->mol[i]->F[0] = 0.0;


      ens->ham[i] = new Hamiltonian_Tully(rep, Ham_indx, ens->mol[i]);

      // A = omega
      // B = E_r
      // C = V
      // D = eps_0
      ens->ham[i]->set_param("C",5.0e-5);

      double dE = (3.0e-2 - 1.5e-2)/16.0;
      ens->ham[i]->set_param("D",1.5e-2+ie*dE);

      

      // Now we start in left diabatic well 
      if(rep==0){
        ens->el[i]->istate = isurface;
      }
      else if(rep==1){
        // We start in left diabatic state - this corresponds to a mixture of two adiabatic states
        // The TD-SE populations are set from the transformation coefficients:

        // Later - set Hamiltonian parameters, if vary them
        ens->ham[i]->set_status(0);
        ens->ham[i]->compute(ens->mol[i]);


        double E0, E1, H0, H1, V;
        E0 = ens->ham[i]->H(ens->mol[i],0,0).real();
        E1 = ens->ham[i]->H(ens->mol[i],1,1).real();
        H0 = ens->ham[i]->H(ens->mol[i],0,0,0).real(); // diabatic
        H1 = ens->ham[i]->H(ens->mol[i],1,1,0).real();
        V  = ens->ham[i]->H(ens->mol[i],0,1,0).real();

//        cout<<"x= "<<ens->mol[i]->R[0].x<<endl;
//        cout<<"E0= "<<E0<<endl;
//        cout<<"E1= "<<E1<<endl;
//        cout<<"H0= "<<H0<<endl;
//        cout<<"H1= "<<H1<<endl;
//        cout<<"V= "<<V<<endl;

        double f0 = V/sqrt((H0 - E0)*(H0 - E0) + V*V);         // f0 - left diabat on 0
        double g0 = (H0 - E0)/sqrt((H0 - E0)*(H0 - E0) + V*V); // g0 - right diabat on 0
        double f1 = V/sqrt((H0 - E1)*(H0 - E1) + V*V);         // f1 - left diabat on 1
        double g1 = (H0 - E1)/sqrt((H0 - E1)*(H0 - E1) + V*V); // g1 - right diabat on 1

//        cout<<"Probability matrix:\n";
//        cout<<" f0(left on 0)= "<<f0*f0<<"    f1(left on 1)= "<<f1*f1<<"\n";
//        cout<<" g0(right on 0)= "<<g0*g0<<"    g1(right on 1)= "<<g1*g1<<"\n";

        // Inverse C:
        // It is easy to prove that transpose is equal to inverse, so we don't need explicit inverse (below)
        // But still, lets use direct inverse:
        double den = g1*f0 - g0*f1;
        double a00 = g1/den;
        double a01 = -f1/den;
        double a10 = -g0/den;
        double a11 = f0/den;

        double ksi = uniform(-2.0*M_PI, 2.0*M_PI);

        // We start on left diabat, so use only f (a01 and a11) coefficients, phases may be arbitrary
//        ens->el[i]->q[0] = a01 * cos(ksi);    ens->el[i]->p[0] = a01 * sin(ksi); 
//        ens->el[i]->q[1] = a11 * cos(ksi);    ens->el[i]->p[1] = a11 * sin(ksi); 

        ens->el[i]->q[0] = f0 * cos(ksi);    ens->el[i]->p[0] = f0 * sin(ksi); 
        ens->el[i]->q[1] = f1 * cos(ksi);    ens->el[i]->p[1] = f1 * sin(ksi); 


//        cout<<"Inverse: Left diabatic state contains:\n";
//        cout<<"f0= "<<f0<<" 0-th adiabatic state\n";
//        cout<<"f1= "<<f1<<" 1-th adiabatic state\n";
//        cout<<"a01= "<<a01<<" 0-th adiabatic state\n";
//        cout<<"a11= "<<a11<<" 1-th adiabatic state\n";

        // Set initial discrete state randomly:
        ksi = uniform(0.0, 1.0);
       
//        ens->el[i]->istate = 0;
        if(ksi<f0*f0){  ens->el[i]->istate = 0; } // 0-th adiabatic state
        else{ ens->el[i]->istate = 1; }

      }


    }// i


//    exit(0);



    // Observables - reflection and transition probabilities on all states
    vector<double> pops(nstates,0.0);


    ens->sh_pop1(pops);
    cout<<"Initial populations: "<<pops[0]<<"  "<<pops[1]<<endl;

//    exit(0);

                     
    //################# Loop over snapshots ######################
    for(int isnap=0; isnap<nsnap; isnap++){

/*
      cout<<"snap = "<<isnap<<endl;

      cout<<" isnap= "<<isnap<<" x= "<<ens->mol[0]->R[0].x
          <<" istate= "<<ens->el[0]->istate
          <<" Fx= "<<ens->mol[0]->F[0].x<<" P01= "<<Pij[0][1]
          <<" q[0]= "<<ens->el[0]->q[0]<<" p[0]= "<<ens->el[0]->p[0]
          <<" q[1]= "<<ens->el[0]->q[1]<<" p[1]= "<<ens->el[0]->p[1]<<endl;
*/

//    exit(0);

      for(int istep=0; istep<nstep; istep++){

        propagate_ensemble_lang(dt, ens, therm, "FSSH"); // lang = Langevin

        //---------- Now do surface hopping stuff -------------
        // Manually
        for(int traj=0;traj<ens->size;traj++){


            if(method==0){
              compute_hopping_probabilities(ens->mol[traj], ens->el[traj], ens->ham[traj], Pij, dt, 0, 300.0 );
            }
            else if(method==1){
              compute_hopping_probabilities_mssh(ens->mol[traj], ens->el[traj], ens->ham[traj], Pij, dt, 0, 300.0 );
            }
            if(method==2){
              compute_hopping_probabilities_decoh(ens->mol[traj], ens->el[traj], ens->ham[traj], Pij, FCij[traj], dt, 0, 300.0 );
            }


            double ksi = uniform(0.0,1.0);

 
          if(method==0||method==1){
            hop(ens->el[traj]->istate, ens->mol[traj], ens->ham[traj], ksi, Pij, do_rescaling, rep);
          }
          else if(method==2){
            hop_decoh(ens->el[traj]->istate, ens->mol[traj], ens->ham[traj], ksi, Pij, FCij[traj], do_rescaling, rep);
          }


        }// for traj
      }// for istep


      
//      ens->sh_pop(pops);
      ens->compute_averages();
      ens->sh_pop1(pops); // specifically for Marcus spin-boson model in adiabatic representation with projection on diabatic states


      // Print probabilities
      stringstream ss(stringstream::in|stringstream::out);  ss << ie;
      std::string ss1;  ss >> ss1;
      ofstream f(("relax"+ss1+".txt").c_str(),ios::app);    
      f<<isnap*nstep*dt<<"  "<<pops[0]<<"  "<<pops[1]<<"  "<<ens->ave_x<<endl;  
      f.close();


    }// for isnap

//    exit(0);

    delete ens;


  }// for i- fake parameter loop


}// run_Marcus




void print_ens(std::string filename, Ensemble* ens){

  ofstream f(filename.c_str(), ios::out);

  for(int i=0;i<ens->size;i++){
    for(int n=0;n<ens->mol[i]->Nnucl;n++){

      f<<ens->mol[i]->R[n]<<endl;

    }
  }

  f.close();

}


void run_double_slit(){

// Brooksby & Prezhdo (from Miller)

  // Simulation parameters
  int nstates = 1;       // number of electronic states
  int Nnucl = 1;         // number of nuclei
  int nsnap = 100;      // number of snapshots
  int Ham_indx = 0;      // 0 - double slit
  int nstep = 5000;       // number of steps per 1 snap;
  int ntraj = 150; 
  int do_rescaling = 1;  // 0 - CPA, 1 - noraml electron back-reaction
  int method = 0;        // 0 - FSSH,  1 - MSSH, 2 - FSSH + FC


  double xstart = -350.0; // initial nuclear coordinate
  double ystart = 0.0;

  // 1 cm^-1 = 4.5564e-6
  double kstart = 0.1366 * 1.0; // sqrt(2048.0 * 4.5564e-6 * 2.0 * 1.0);   // initial nuclear momentum - corresponds to Ekin = 2048 cm^-1

  int isurface = 0;      // initial electronic state
  int rep = 1;           // 0 - diabatic, 1 = adiabatic
  double dt = 0.1;  //prms.dt; a.u.

  int i,j,k;

  double V0 = 0.0364512;     // 8000 cm^-1 
  double omega = 0.00273384; // 600 cm^-1  
  double alpha = 50.0;       // in a.u.





  // Simulation variables  
  Ensemble* ens; ens = new Ensemble(ntraj, nstates, Nnucl);


  // Initialization
  for(i=0;i<ntraj;i++){

  // Note: if  <q^2> - <q>^2 = s^2 ==>  <p^2> - <p>^2 = hbar^2 / (4*s^2)

    // Init nuclear variables - beyond default constructor
    ens->mol[i]->mass[0] = 1.0; // electron mass
    ens->mol[i]->R[0].x = xstart + 50.0*normal();  
    ens->mol[i]->R[0].y = ystart + 279.34*normal();  // sqrt(16.0*V0/(1.0*omega*omega));  
    ens->mol[i]->R[0].z = 0.0;

    ens->mol[i]->P[0].x = kstart + (0.5/50.0)*normal();
    ens->mol[i]->P[0].y = 0.0 + 0.0*(0.5/279.34)*normal();  // It is essentiall to make some distribution in this direction
    ens->mol[i]->P[0].z = 0.0;

    ens->mol[i]->F[0] = 0.0;
    ens->el[i]->istate = isurface;


    ens->ham[i] = new Hamiltonian_Double_Slit(ens->mol[i]);



  }// i

  ens->compute_averages();

  // "Screen"
  double Xscreen = 1250; 
  double Ymin = -1200.0; // 
  double phi_min = atan2(Ymin,Xscreen);
  double Ymax =  1200.0;
  double phi_max = -phi_min;// atan2(Ymax,Xscreen);
  double dy = 25.0; 
  double dphi = 0.05; // radian

//  int binN = floor((Ymax-Ymin)/dy);
  int binN = floor((phi_max-phi_min)/dphi);
  vector<double> binY(binN,0.0);
  vector<double> binTeta(binN,0.0);
  vector<double> binC(binN,0.0);  // count

  for(i=0;i<binN;i++){ 
//    binY[i] = Ymin + i*dy;
//    binTeta[i] = atan(binY[i]/Xscreen);
    binTeta[i] = phi_min + i*dphi;
    binY[i] = Ymin + Xscreen*tan(binTeta[i]);

  }

  vector<int> traj_state(ntraj,0); // 0 - before screen, 1 - after screen

  double fl_Xmax = 10000.0;
  double fl_Xmin = -500.0;
  double fl_dx = 20.0;
  double fl_Ymax = 500.0;
  double fl_Ymin = -500.0;
  double fl_dy = 10.0;

  int Nx = floor((fl_Xmax - fl_Xmin)/fl_dx);
  int Ny = floor((fl_Ymax - fl_Ymin)/fl_dy);
  vector< vector<double> > flux(Nx,vector<double>(Ny,0.0)); // Total (integrated over time flux density on the grid)

  

  std::string prefix = "res/wfc";

  ofstream en("energy.txt",ios::out);


  double Ecorr = 0.0;
  double counter = 0.0;
                     
  //################# Loop over snapshots ######################
  for(int isnap=0; isnap<nsnap; isnap++){

    cout<<"Printing snap "<<isnap<<endl;
//    ens->print_map(prefix, -500.0, 500.0, 10.0, -500.0, 500.0, 10.0, isnap);
    ens->print_map(prefix, -500.0, 2500.0, 20.0, -1000.0, 1000.0, 20.0, isnap);


    for(int istep=0; istep<nstep; istep++){

      counter += 1.0;

      Ecorr = 0.0;
      propagate_ensemble(dt, ens, "FSSH"); // classical description
//      propagate_ensemble4(dt, ens, "ENT", Ecorr); // classical description

      ens->integral_flux(flux, fl_Xmin, fl_Xmax, fl_dx, fl_Ymin, fl_Ymax, fl_dy);


/*
          if(ens->mol[traj]->R[0].x<0.1 && ens->mol[traj]->R[0].x>-0.1){

          cout<<"traj= "<<traj<<" ksi= "<<ksi<<" istate= "<<ens->el[traj]->istate<<endl;
          cout<<"x = "<<ens->mol[traj]->R[0].x<<endl;
          for(i=0;i<Pij.size();i++){
            for(j=0;j<Pij[i].size();j++){
              cout<<Pij[i][j]<<" ";
            }
            cout<<endl;
          }
          }
*/

      // Count the trajectories that end up crossing the screen
      for(int traj=0;traj<ens->size;traj++){

        if(ens->mol[traj]->R[0].x > Xscreen && traj_state[traj]==0){  
          if(ens->mol[traj]->R[0].y > Ymin && ens->mol[traj]->R[0].y < Ymax){ 

//            int indx = floor((ens->mol[traj]->R[0].y - Ymin)/dy);
            double teta = atan(ens->mol[traj]->R[0].y/Xscreen);
            int indx = floor((teta - phi_min)/dphi);

            binC[indx] += 1.0;

          }

          traj_state[traj] = 1;
        }

      }// for traj


      // Deactivate trajectories:
      for(int traj=0;traj<ens->size;traj++){
        if(ens->mol[traj]->R[0].x <-500 || ens->mol[traj]->R[0].x > 500.0) { ens->is_active[traj] = 0; }
      }


    }// for istep
     
//    Epot = compute_potential_energy(ens->mol[0], ens->el[0], ens->ham[0], "FSSH");
//    Ekin = compute_kinetic_energy(ens->mol[0]);

//      cout<<"px = "<<kstart<<" isnap= "<<isnap<<" F= "<<ens->mol[0]->F[0]<<"  R= "<<ens->mol[0]->R[0]<<" P= "<<ens->mol[0]->P[0]<<endl;
/*
      cout<<"px = "<<kstart<<" isnap= "<<isnap<<" x= "<<ens->mol[0]->R[0].x
          <<" Ekin= "<<Ekin<<" Epot= "<<Epot<<" Etot= "<<Ekin+Epot
          <<" istate= "<<ens->el[0]->istate
          <<" Fx= "<<ens->mol[0]->F[0].x<<" P01= "<<Pij[0][1]
          <<" q[0]= "<<ens->el[0]->q[0]<<" p[0]= "<<ens->el[0]->p[0]
          <<" q[1]= "<<ens->el[0]->q[1]<<" p[1]= "<<ens->el[0]->p[1]<<endl;
*/


    double Epot,Ekin,Etot;
    compute_energies(ens, Epot,Ekin, Etot, "FSSH");
   
    en<<isnap<<"  "<<Ekin<<"   "<<Epot<<"  "<<Etot<<"  "<<Ecorr<<"  "<<Etot+Ecorr<<endl;


    ofstream flux_map("flux_map.txt",ios::out);
    for(int nx=0;nx<Nx;nx++){
      for(int ny=0;ny<Ny;ny++){
          flux_map<<(fl_Xmin + nx*fl_dx)<<"  "<<(fl_Ymin + ny*fl_dy)<<"  "<<(flux[nx][ny]/counter)<<"\n";
      }
      flux_map<<"\n";
    }
    flux_map.close();
    


  }// for isnap

  en.close();


  // Now renormalize distributions:
  double n_hits = 0.0;
  for(i=0;i<binN;i++){
    n_hits += binC[i];
  }
  if(n_hits<1){ cout<<"Warning: No trajectory reached the screen\n"; }

  for(i=0;i<binN;i++){
    binC[i] /= n_hits;
  }


 
  // Compute and print probabilities 
  
  ofstream f("slit.txt",ios::out);      

  for(i=0;i<binN;i++){
    f<<i<<"  "<<binY[i]<<"  "<<binTeta[i]<<"  "<<binC[i]<<endl;
  }

  f.close();


  print_ens("final_snap.txt", ens);


  delete ens;



}// run_double_slit


void run_double_slit_exact(){

// Brooksby & Prezhdo (from Miller)

  // Simulation parameters
  int nstates = 1;       // number of electronic states
  int nsnap = 100;      // number of snapshots
  int nstep = 50;       // number of steps per 1 snap;


  double xstart = -350.0; // initial nuclear coordinate
  double ystart = 0.0;

  // 1 cm^-1 = 4.5564e-6 a.u.
  double kstart = 0.1366 * 1.0; // sqrt(2048.0 * 4.5564e-6 * 2.0 * 1.0);   // initial nuclear momentum - corresponds to Ekin = 2048 cm^-1

  int isurface = 0;      // initial electronic state
  double dt = 10.0;  //prms.dt; a.u.

  int i,j,k;

  double V0 = 0.0364512;     // 8000 cm^-1 
  double omega = 0.00273384; // 600 cm^-1  
  double alpha = 50.0;       // in a.u.




  // Simulation variables  
  Grid_wfc* wfc;
  wfc = new Grid_wfc(-2500.0, 2500.0, 50.0, -500.0, 500.0, 10.0,  1);

  cout<<"Properties of the grid:\n";
  cout<<"Nx= "<<wfc->Nx<<"  xmin= "<<wfc->xmin<<"  xmax= "<<wfc->xmax<<"  dx=  "<<wfc->dx<<endl;
  cout<<"Ny= "<<wfc->Ny<<"  ymin= "<<wfc->ymin<<"  ymax= "<<wfc->ymax<<"  dy=  "<<wfc->dy<<endl;
//  exit(0);

  // dx = 50.0
  // dy = 279.346
  wfc->init_wfc(xstart, ystart, 
                kstart, 0.0,  
                alpha,  1.0*sqrt(16.0*V0/(1.0*omega*omega)),
              0);

  // This does not make much sene, only allocation of memory and assignment of mass
  wfc->mol = new Nuclear(1);  
  wfc->mol->mass[0] = 1.0;

  // Allocation of generic Hamiltonian - for internal calculations
  wfc->ham = new Hamiltonian_Double_Slit(wfc->mol);
  wfc->ham->set_status(0);

  wfc->update_potential();
  wfc->update_propagator(0.5*dt,wfc->mol->mass[0]);
  wfc->update_propagator_K(0.5*dt,wfc->mol->mass[0]); // !!! NOTE: dt/2 not dt !!! - depends on splitting scheme


//  exit(0);

  // "Screen"
  double Xscreen = 250; 
  double Ymin = -800.0; // 
  double phi_min = atan2(Ymin,Xscreen);
  double Ymax =  800.0;
  double phi_max = -phi_min;// atan2(Ymax,Xscreen);
  double dy = 25.0; 
  double dphi = 0.05; // radian

  int binN = floor((phi_max-phi_min)/dphi);
  vector<double> binY(binN,0.0);
  vector<double> binTeta(binN,0.0);
  vector<double> binC(binN,0.0);  // count

  for(i=0;i<binN;i++){ 
    binTeta[i] = phi_min + i*dphi;
    binY[i] = Ymin + Xscreen*tan(binTeta[i]);

  }


  std::string prefix = "res/wfc";
  ofstream en("energy.txt",ios::out);
                     
  //################# Loop over snapshots ######################
  for(int isnap=0; isnap<nsnap; isnap++){

    cout<<"Printing snap "<<isnap<<endl;
    wfc->print_map(prefix, isnap, 0);

//    cout<<wfc->PSI[0]<<endl;

    for(int istep=0; istep<nstep; istep++){

    
      propagate_exact_2D(dt, wfc, 1);


/*
      // Count the trajectories that end up crossing the screen
      for(int traj=0;traj<ens->size;traj++){

        if(ens->mol[traj]->R[0].x > Xscreen && traj_state[traj]==0){  
          if(ens->mol[traj]->R[0].y > Ymin && ens->mol[traj]->R[0].y < Ymax){ 

            double teta = atan(ens->mol[traj]->R[0].y/Xscreen);
            int indx = floor((teta - phi_min)/dphi);

            binC[indx] += 1.0;

          }

          traj_state[traj] = 1;
        }

      }// for traj
*/


    }// for istep
     

    double Epot,Ekin,Etot,Ecorr;   
    en<<isnap<<"  "<<Ekin<<"   "<<Epot<<"  "<<Etot<<"  "<<Ecorr<<"  "<<Etot+Ecorr<<endl;


  }// for isnap

  en.close();


  // Now renormalize distributions:
  double n_hits = 0.0;
  for(i=0;i<binN;i++){
    n_hits += binC[i];
  }
  if(n_hits<1){ cout<<"Warning: No trajectory reached the screen\n"; }

  for(i=0;i<binN;i++){
    binC[i] /= n_hits;
  }


 
  // Compute and print probabilities 
  
  ofstream f("slit.txt",ios::out);      

  for(i=0;i<binN;i++){
    f<<i<<"  "<<binY[i]<<"  "<<binTeta[i]<<"  "<<binC[i]<<endl;
  }

  f.close();


}// run_double_slit_exact





int main(){

  srand(time(0));

//  test1();
//  test3();
//  run_semiclassical();
//  run_double_slit();
//  run_double_slit_exact();

  run_Marcus();


  return 0;
}
