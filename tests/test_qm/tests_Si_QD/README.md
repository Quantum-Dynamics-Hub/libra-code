# Electronic structure calculations on Si QDs

## Intro

  This example demonstrates how run electronic structure calculations (at the EHT
  level) on various Si QDs. The postprocessing includes computing partial densities of 
  states and molecular orbitals (for visualization with VMD)


## Files

  * control_parameters_eht.dat  - defines the simulation parameters for EHT calculations (if EHT is chosen)
  * control_parameters_indo.dat  - defines the simulation parameters for INDO calculations (if INDO is chosen)
  * elements.dat - file that defines the properties of elements in the periodic table
  * muller_params - parameters for the EHT Hamiltonian
  * params_indo - parameters for the INDO Hamiltonian
  * plot_dos.plt - the file to plot pDOS of the system
  * QD_6_h.inp.xyz,  ... QD_12_h.inp.xyz - are sample input files 
  * test_qd.py - the main file to run the calculations: define which Hamiltonan to use, which system to compute,
  which projections to compute, which orbitals to plot, and all other relevant parameters here


## Procedures

  * Edit the file test_qd.py as needed
  * Create the output directory, e.g. qd6, if absent
  * Run the calculations:  `python test_qd.py`
  * The first line of the generated file dos_proj.txt will contain 
    the value of the Fermi energy - use it to adjust the x-axis of the 
    plotting script 'plot_dos.plt'
  * Run the gnuplot to compute the pDOS:
    `gnuplot plot_dos.plt`
  * Open VMD and its Tcl shell, navigate in the Tcl shell to 
    the results directory (e.g. qd6)
  * Load the .cube files and plot visualize the structure and orbitals

    
