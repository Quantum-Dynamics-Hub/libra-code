


double compute_kinetic_energy(Nuclear* mol){
/**
  \brief Compute kinetic energy of Nuclear object 
  \param[in] mol The pointer to the Nuclear DOF object

  This is the classical nuclear kinetic energy
*/

  double Ekin = 0.0;

  for(int i=0;i<mol->nnucl;i++){  Ekin += mol->p[i]*mol->p[i]/mol->mass[i];   }
  Ekin *= 0.5;

  return Ekin;

}

double compute_kinetic_energy(Nuclear& mol){  
/**
  \brief Compute kinetic energy of Nuclear object - Python-friendly
  \param[in] mol The reference to the Nuclear DOF object

  This is the classical nuclear kinetic energy
*/

  return compute_kinetic_energy(&mol); 
}

double compute_kinetic_energy(Ensemble& ens){
/**
  \brief Compute average kinetic energy of an ensemble of Nuclear objects - Python-friendly
  \param[in] ens The reference to the Ensemble object

  This is the classical nuclear kinetic energy
*/


  double res = 0.0;
  for(int traj=0;traj<ens.ntraj;traj++){
    res += compute_kinetic_energy(&ens.mol[traj]);
  }
  res /= (double)ens.ntraj;
  
  return res;
}



double compute_potential_energy(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt){
/**
  \brief Compute potential energy of a quantum-classical system
  \param[in,out] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
*/

  int i,j;
  double Heff = 0.0;
  double Epot = 0.0;


  // Calculate all surfaces,
  ham->set_q(mol->q);
  ham->compute();


  // Mean-field/Ehrenfest mixing
  if(opt==0){

    // Potential energy - electronic Hamiltonian
    // working in a.u., so hbar = 1
    for(i=0;i<el->nstates;i++){
      for(j=0;j<el->nstates;j++){
        Heff += 0.5*ham->Hvib(i,j).real() * (el->q[i] * el->q[j] + el->p[i] * el->p[j]);
        Heff += ham->Hvib(i,j).imag() * el->p[i] * el->q[j]; 
      }// for j
    }// for i

    Epot = 2.0*Heff;  // in a.u. - hbar is assumed = 1
  }// algo = "MF"

  // FSSH mixing 
  else if(opt==1){   
    i = el->istate; // current electornic state
    Epot = ham->H(i,i).real();  
  }

  return Epot; 

}

double compute_potential_energy(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt){
/**
  \brief Compute potential energy of a quantum-classical system - Python-friendly
  \param[in,out] mol Describes the nuclear DOF
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0

*/

  return compute_potential_energy(&mol, &el, &ham, opt);

}

double compute_potential_energy(Ensemble& ens, int opt){
/**
  \brief Compute potential energy of a quantum-classical system - Python-friendly
  \param[in,out] ens The reference to the ensemble of trajectories
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0

  The returned potential energy is also averaged over all trajectories included in the ensemble object ens.

*/

  double res = 0.0;
  for(int traj=0;traj<ens.ntraj;traj++){
    res += compute_potential_energy(&ens.mol[traj], &ens.el[traj], ens.ham[traj], opt);
  }
  res /= (double)ens.ntraj;
  
  return res;
}



double compute_forces(Nuclear* mol, Electronic* el, Hamiltonian* ham, int opt){
/**
  \brief Compute potential energy and forces (an all classical DOF) of a quantum-classical system
  \param[in,out] mol Describes the nuclear DOF, is modified during calculations to update forces
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
  The computed forces are stored in the mol object
*/


  int i,j,k;
  double Heff = 0.0;
  double Epot = 0.0;


  // Calculate all surfaces, if needed
  ham->set_q(mol->q);
  ham->compute();

  // Zero all forces in mol
  for(i=0;i<mol->nnucl;i++){ mol->f[k] = 0.0; }

  
  // Mean-field/Ehrenfest mixing
  if(opt==0){

    // Potential energy - electronic Hamiltonian
    // working in a.u., so hbar = 1
    for(i=0;i<el->nstates;i++){
      for(j=0;j<el->nstates;j++){

        double cij_re = (el->q[i] * el->q[j] + el->p[i] * el->p[j]);
        double cij_im = el->p[i] * el->q[j];

        Heff += 0.5*ham->Hvib(i,j).real() * cij_re;
        Heff += ham->Hvib(i,j).imag() * cij_im; 

        for(k=0;k<mol->nnucl;k++){
          
          mol->f[k] -= 2.0 * 0.5*ham->dHdq(i,j,k).real() * cij_re;
          mol->f[k] -= 2.0 * ham->dHdq(i,j,k).imag() * cij_im; 

        }// for k

      }// for j
    }// for i

    Epot = 2.0*Heff;  // in a.u. - hbar is assumed = 1

  }// algo = "MF"

  // FSSH mixing 
  else if(opt==1){

    i = el->istate; // current electornic state
    Epot = ham->H(i,i).real();  

    for(k=0;k<mol->nnucl;k++){  
      mol->f[k] = -ham->dHdq(i,i,k).real(); 
    }// for k

  }// opt == 1

  return Epot;

}

double compute_forces(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt){
/**
  \brief Compute potential energy and forces (an all classical DOF) of a quantum-classical system - Python-friendly
  \param[in,out] mol Describes the nuclear DOF, is modified during calculations to update forces
  \param[in] el Describes electronic DOF
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
  The computed forces are stored in the mol object
*/

  return compute_forces(&mol, &el, &ham, opt);

}

double compute_forces(Ensemble& ens, int opt){
/**
  \brief Compute potential energy and forces (an all classical DOF) of a quantum-classical system
  \param[in,out] ens The reference to the ensemble of trajectories
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  Returns the potential energy computed according to the selected option: either the poetial energy of an active electronic level
  (as given by el variable) - this is for opt = 1 - or the state averaged potential energy (the averaging also accounts for
  nonadiabatic couplings or off-diagonal Hamiltonian terms) - this is for opt = 0
  The computed forces are stored in the mol object

  The returned potential energy is also averaged over all trajectories included in the ensemble object ens.
  No averaging of forces is needed.

*/


  double epot = 0.0;

  for(int traj=0;traj<ens.ntraj;traj++){
    epot += compute_forces(&ens.mol[traj], &ens.el[traj], ens.ham[traj], opt);
  }
  epot = epot/float(ens.ntraj);

  return epot;
  
}



void compute_energies(Ensemble* ens, double& Epot, double& Ekin, double& Etot,int opt){
/**
  \brief Compute ensemble-averaged potential, kinetic and total energy
  \param[in,out] ens The pointer to Ensemble object for which we want to compute properties
  \param[out] Epot The computed averaged potential energy will be stored here
  \param[out] Ekin The computed averaged kinetic energy will be stored here
  \param[out] Etot The computed averaged total energy will be stored here
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

  
*/
  Epot = 0.0;
  Ekin = 0.0;
  Etot = 0.0;
  
  for(int traj=0;traj<ens->ntraj;traj++){

    double ek = compute_kinetic_energy(ens->mol[traj]);
    double ep = compute_potential_energy(&ens->mol[traj], &ens->el[traj], ens->ham[traj], opt);

    Epot += ep;
    Ekin += ek;
  } 

  Epot /= ((double)ens->ntraj);
  Ekin /= ((double)ens->ntraj);
  Etot = Ekin + Epot;

}



