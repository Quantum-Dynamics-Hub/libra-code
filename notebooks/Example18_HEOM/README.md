# 1. Reference data

   Calculations done with Amber Jain's & Joe Subotnik's code
   [ original code ](https://github.com/subotnikgroup/HEOM_Amber)
   [ clone on the quantum dynamics hub ](https://github.com/Quantum-Dynamics-Hub/HEOM_Amber)

   *spectra_heom.inp* - is the input file
   *rho.out* - is the output file, containing the populations of 2 states vs. time

   !! Reproducing figure 5 of Strumpfer Schulten, JCTC 8, 2808 (2012)
   !! Chen, Zheng, Shi, Tan, JCP 131, 094502 (2009)


# 2. Methodology development approach

   *test_heom.py* - the Python code that shows the steps of the calculations
    
   it defines its own functions that can be used in the execution, but currently
   the C++ functions exposed to Python from Libra core are used

   This script isn't the best example, since some of the functions use the global-scope variables,
   but it works and illustrates nicely the idea (i hope). 

   *pops.txt* - is the output file produced by *test_heom.py*

   Once produced, you can plot the curves pops.txt vs. rho.out using the *HEOM_plotting.ipynb* notebook


# 3. User-oriented  approach 

   To use the HEOM code in a most convenient way, simply use the *HEOM_via_libra_py.ipynb* notebook

   All you need to worry about there - just some properties describing the bath and the system.

   The code generates the */out* directory with the file that stores the parameters used in the calculations and 
   the HDF5 file with the simulation results.

   The notebook also uses the created data to plot the results (aginst the reference data obtained with the original code).


