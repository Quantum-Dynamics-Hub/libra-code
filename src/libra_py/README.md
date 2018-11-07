# Description of the modules

  * normal_modes.py  - Normal modes analysis 

   - covariance_matrix(X, M) 
   - visualize_modes(E, R, U, w, params)
   - compute_cov(R, V, A, M, E, params)
     [Ref: Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics 
      or Monte Carlo Simulations. J. Chem. Phys. 2003, 120, 1–4]
   - compute_cov2(R, A, E, params)
     [Ref: Pereverzev, A.; Sewell, T. D. Obtaining the Hessian from the Force Covariance Matrix:
      Application to Crystalline Explosives PETN and RDX. J. Chem. Phys. 2015, 142, 134110]
   - compute_dynmat(R, D, E, params)


  * LAMMSP_methods.py  - Utility functions to use with LAMMPS
  
   - compute_dynmat(lmp, filename, atoms, dr)