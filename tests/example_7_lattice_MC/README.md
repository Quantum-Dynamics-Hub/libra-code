# Lattice Monte Carlo (MC) Example

## Intro

  This example demonstrates how to run MC-base annealing of 
  a 2D lattice. The lattice is prepared with the dopants randomly distribute.
  You define the model Hamiltonian and run the Metropolis-based MC simulations
  to get the structure with the lowest possible energy.


## Files

  * ar.xyz  - defines the lattice unit cell, the atom type doesn't matter
  * elements.dat - file that defines the properties of elements in the periodic table
  * make_pics.tcl - file to generate a set of pictures for all structures along the simulation
  * model.py - the main file to run the calculations: define the lattice, Hamiltonian and 
  the simulation parameters here.
  * plot.plt - the file to plot the energy of the system as the function of annealing time

## Procedures

  * Edit the file model.py as needed
  * Create the output directory, e.g. 2D, if absent
  * Run the calculations:  `python model.py`
  * Copy plotting files to the results folder:
    `cp plot.plt 2D` and `cp make_pics.tcl 2D`
  * In the results directory, run the gnuplot:
    `gnuplot plot.plt`
  * Open VMD and its Tcl shell, navigate in the Tcl shell to 
    the results directory (2D)
  * Run the make_pics.tcl script: `source make_pics.tcl`

    
