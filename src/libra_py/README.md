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

  * acf_matrix.py  - A file containing various schemes to compute the ACF of a given data set 
                     with data elements of the type MATRIX.
   - acf(data)
     Computes the autocovariance function using the method with the least bias
     [Ref: https://www.itl.nist.gov/div898/handbook/eda/section3/eda331.htm]
   - acf2(data)
     Computes the autocovariance function using the method with more bias, but is less
     prone to noise for large values of "k"
     [Ref: https://www.itl.nist.gov/div898/handbook/eda/section3/eda331.htm]

