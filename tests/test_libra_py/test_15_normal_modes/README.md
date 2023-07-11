# Tutorials on Normal Modes

## Example 1

### Instructions

  1.  Just run the script ```python test.py ```

  2.  Load any of the produced output .xyz files to the VMD to visualize
      You can use the pre-defined tcl scripts, e.g. ```source load_Si8.tcl``` in the 
      TkConsole of VMD

### Explanations 

  1. We use the methods ```R, V, A, M, E = QE_methods.read_md_data("x0.xml") ``` to read the 
  data of the MD simulations performed with QE and stored in the XML file. The matrices R, V, A, M
  with coordinates, velocities, accelerations, and masses are generated. In addition, we obtain a list
  E of the atom type labels (element names). This is all the info we need for further calculations

  2. We use two types of approaches to compute the normal modes and corresponding frequencies.
  One is according to [(1) Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics 
  or Monte Carlo Simulation. J. Chem. Phys. 2003, 120, 1-4.] which is invoked with
  ```normal_modes.compute_cov( R, V, A, M, E, params)``` and the other is according to 
  [(2) Pereverzev, A.; Sewell, T. D. Obtaining the Hessian from the Force Covariance Matrix:
  Application to Crystalline Explosives PETN and RDX. J. Chem. Phys. 2015, 142, 134110.] which is invoked
  with ```normal_modes.compute_cov2( R, A, M, E, T, params)```  
  Each function does the following:
  - computes the frequencies
  - prints out the selected normal modes as .xyz for further visualization

  3. Note that the function *compute_cov* tries to compute the normal modes/frequencies according to
  two approaches discussed in the corresponding paper - one based on the covariances of V and R, and another
  based on covariances of A and R. The results are printed in the files with "_velocities" and "_accelerations"
  correspondingly. 

  4. The function *compute_cov2* requires the temperature in its definition, which should be chosen to be the 
  average temperature of the simulation.

  5. Each function is called with the one of the two possible values for the flag "cov_flag", which controls how the
  covariance matrices are computed. The standard definition of the covariance matrix subtracts the mean value of each 
  coordinate (so one actually looks at the covariance of the fluctuations). However, both papers define their "covariances"
  without such a correction, so I'm not totally sure which option is the best and leave it to the user to decide. 


### Systems

  1. x0.xml - Si8 molecule simulated for 100 fs

  2. x0_h2.xml - H2 molecule simulated for 500 fs

  3. x0_co2.xml - CO2 molecule simulated for 1000 fs





## Example 4

### Instructions

  1.  Just run the script ```python test.py ```

### Explanations 

  This exmple demonstrates how to get the coordinates, element names, and normal modes for a given system from
  QE (phonon) output file (.dyn). This example simply shows how to invoke the processing of two files - the first one
  produces some verbose output on how thigs are going under the hood - good for demonstrations. The other file
  deals with a much larger system, so we don't really want to produce a lot of output for it - just get the objects with
  the data we need. For this, we use the default value of verbosity flag (0 - no output). 

  The objects generated can later be used for visualization of the normal modes, with for instance the
  normal_modes.get_xyz2 function
