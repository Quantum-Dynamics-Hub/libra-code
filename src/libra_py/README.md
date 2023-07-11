# Principles of function's arguments ordering:

 Function arguments list:

 1. dynamical variables or DOFs properties
 2. parameters controlling dynamics 
 3. function to compute Hamiltonians
 4. parameters of the Hamiltonian function
 5. parameters controlling the initialization
 6. random object
 7. additional variables that can be setup to default values


 For the dynamical variables:

 1. dynamical variables 
 2. DOF properties

 Within each group, when applies:

 1. nuclear DOFs
 2. electronic DOFs


 Again, whithin each group when applies:

 1. continuous variables
 2. discrete variables






# Description of the modules

## acf_matrix.py  - Calculation of ACF of a given data set with data elements of the type MATRIX.
   * average(data) - compute the average of the data [REDUNDANCY 1]
   * center_data(data) - - center data around the time-series average [REDUNDANCY 2]
   * acf(data)
     Computes the autocovariance function using the method with the least bias
     [Ref: https://www.itl.nist.gov/div898/handbook/eda/section3/eda331.htm]
   * acf2(data)
     Computes the autocovariance function using the method with more bias, but is less
     prone to noise for large values of "k"
     [Ref: https://www.itl.nist.gov/div898/handbook/eda/section3/eda331.htm]
   * ft(acf_data, wspan, dw, dt) - compute the Fourier Transform of given data (which is often the ACF) [REDUNDANCY 3]
   * recipe1(data, dt, wspan, dw, acf_filename="acf.txt", spectrum_filename="spectrum.txt", do_center=1)
   * recipe2(data, dt, wspan, dw, acf_filename="acf.txt", spectrum_filename="spectrum.txt", do_center=1)

## acf_vector.py -  Calculation of ACF of a given data set with data elements of the type VECTOR.
   * average(data) - compute the average of the data [REDUNDANCY 1]
   * center_data(data) - center data around the time-series average [REDUNDANCY 2]
   * acf(data,dt) - compute the autocorrelation function of the data
   * ft(acf_data, wspan, dw, dt) - compute the Fourier Transform of given data (which is often the ACF) [REDUNDANCY 3]
   * recipe1(data, dt, wspan, dw, acf_filename="acf.txt", spectrum_filename="spectrum.txt", do_center=1)


## autoconnect.py - Compute the connectivities between atoms
   * autoconnect(R, MaxCoord, Rcut, opt=0, verbosity=0) - computes the connectivities
   * autoconnect_pbc(R, MaxCoord, Rcut, tv1, tv2, tv3, pbc_opt, opt=0, verbosity=0) - as above + potentially PBC
   * example_1() - example

## build.py  - Creation of molecular/atomistic models
   * read_xyz(filename, inp_units=1, out_units=0) - read the xyz file 
   * make_xyz(L, R, inp_units=0, out_units=1) - returns an xyz-formatted string representing the molecular system
   * read_xyz_crystal(filename, a,b,c, inp_units=0, out_units=0) - same + periodic translation

   * generate_replicas_xyz(tv1,tv2,tv3, rep1, rep2, rep3 , filename, outfile, inp_units=0, out_units=0) - another version
   * generate_replicas_xyz2(L, R, tv1, tv2, tv3, Nx, Ny, Nz, inp_units=0, out_units=0) -  multiply the atoms into periodic images

   * crop_sphere_xyz(infile, outfile, Rcut) - another version
   * crop_sphere_xyz2(L, R, Rcut) - remove atoms outside a sphere
   * crop_sphere_xyz3(L, R, Rcut, pairs, new_L) - alternative version

   * add_atoms_to_system(syst, L, R, scl, mass, univ_file) - add atoms to a system
   * add_atom_to_system(syst, coords, MaxCoords, Nx,Ny,Nz, a,b,c, shift, elt, mass, scl, max_coord, rnd) - add an atom and its periodic images to a system

## datautils.py - some utility functions for data prosessing

   * show_matrix_splot(X, filename, set_diag_to_zero=0) - output a matrix in a gnuplot format 
   * find_maxima(s, verbose=0, logfile="run.log") - find maxima in the data and sort them out
   * scalar_stat(X) - get the descriptive statistics of the scalar data
   * matrix_stat(X) - get the descriptive statistics of the matrix data
   * matrix_freqs(X, a, b, dt, prefix, Nfreqs, verbose = [1,1,1], dw = 1.0, wspan = 3000.0, logfile=None) - get the frequencies of the data


## fit.py - fitting and regression module
   * Regression(X,Y) - linear regression
   * fit_exp(X,Y, x0, verbose=1) - fit data to an exponential function
   * fit_gau(X,Y, x0, verbose=1) - fit data to an gaussian function
   * get_data_from_file(filename, xindx, yindx, xminval=None, xmaxval=None, yminval=None, ymaxval=None) - read data from file

## hpc_utils.py - unitilities for calculations on HPC clusters (job management, etc.)
   * make_submit(Nstart,Nend,job_dir,submit_templ) - make submit files from a template
   * job(Nstart,Nend,job_dir,prefixes) - prepares a job 
   * distribute(Nmin,Nmax,max_steps,submit_templ,exp_files,prefixes,do_submit) - prepare submit files, organize filesystem, submit jobs


## hungarian.py - solves an assignment problem using the hungarian method
   * minimize(_X, verbosity=0) - minimizes the "cost" function
   * maximize(_X, verbosity=0) - maximizes the "cost" function


## LAMMSP_methods.py  - Utility functions to use with LAMMPS  
   * compute_dynmat(lmp, filename, atoms, dr)



#=========== Add more modules here =============


## namd.py  - functions for NA-MD calculations
   * compute_Hvib(ham_old, ham_cur, orb, dt) - vibronic Hamiltonian from two Hamiltonian objects in the MO-LCAO basis
   * compute_Hvib_sd(ham_old, ham_cur, orb, SD_basis, dt) -vibronic Hamiltonian from two Hamiltonian objects in the Slater Determinant basis

   
## normal_modes.py  - Normal modes analysis 
   * covariance_matrix(X, M) 
   * visualize_modes(E, R, U, w, params)
   * compute_cov(R, V, A, M, E, params)
     [Ref: Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics 
      or Monte Carlo Simulations. J. Chem. Phys. 2003, 120, 1–4]
   * compute_cov1(R, V, M, E, params) - same as compute_cov, but without acceleration matrix
     [Ref: Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics 
      or Monte Carlo Simulations. J. Chem. Phys. 2003, 120, 1–4]
   * compute_cov2(R, A, E, params)
     [Ref: Pereverzev, A.; Sewell, T. D. Obtaining the Hessian from the Force Covariance Matrix:
      Application to Crystalline Explosives PETN and RDX. J. Chem. Phys. 2015, 142, 134110]
   * compute_dynmat(R, D, E, params)
   * get_xyz(E, R, M, U, mode) - the auxiliary function to generate an xyz in the format good for py3Dmol visualization
   * get_xyz2(E, R, U, mode) - similar to get_xyz, but assumes the normal modes are already include the mass-factors

## nve_md.py - functions for doing NVE MD in a maximally-Pythonic way 
   * nve_md_init(syst, mol, el, ham) - initialize an NVE MD run 
   * nve_md_step(syst, mol, el, ham, dt=20.0, integrator="DLML") - do a single MD step


## parse_gamess.py - functions for parsing the GAMESS input/output/data 
  * name2no(name) - conversion of atomic labels to numbers (from PyQuante)
  * parse_gamess_basis(filename) - parse the file with the GAMESS-type format of the atomic basis data


## pdos.py - handling projected densities of states
  * def convolve(X0, Y0, dx0, dx, var) - convolve data with Gaussian distribution
  * QE_pdos(prefix, emin, emax, de, projections, Ef, outfile, do_convolve, de_new, var) - compute projected density of states (for QE)


## probabilities.py - the module that implements probabilities for certain events in statistical mechanics
  * Boltz_prob_up(E, T) - probability to have kinetic energy more or equal to E at temperature T
  * HO_prob(E, qn, T) - probability to find a multi-dimensional HO in one of its states
  * HO_prob_up(E, qn, T) - probability that each of the HO occupies states above given quantum number


## QE_methods.py - functions for helping with QE calculations
  * cryst2cart(a1,a2,a3,r) - conversion from crystal to Cartesian coordinates
  * read_qe_index(filename, orb_list, verbose=0) - read the QE index file (with all the essential info about the run)
  * read_qe_wfc_info(filename, verbose=0) - read the metadata about the wavefunction from the wavefunction file
  * read_qe_wfc_grid(filename, verbose=0) - read the wavefunction grid (may be too large and hence very slow!)
  * read_qe_wfc(filename, orb_list, verbose=0) - read the wavefunctions in their planewave representation
  * read_md_data(filename) - read the info about the MD run, saved in the xml file 
  * read_md_data_xyz(filename, PT, dt) - read the info about the MD run stored in xyz format
  * read_md_data_cell(filename) - read the cell dimensions at all times
  * out2inp(out_filename,templ_filename,wd,prefix,t0,tmax,dt) - convert the output file with MD trajectory to a set of QE input files
  * out2pdb(out_filename,T,dt,pdb_prefix) - convert the output file with MD trajectory to a set of PDB files
  * out2xyz(out_filename,T,dt,xyz_filename) - convert the output file with MD trajectory to an XYZ trajectory file
  * xyz2inp(out_filename,templ_filename,wd,prefix,t0,tmax,dt) - convert the XYZ MD trajectory file to a set of QE input files
  * get_QE_normal_modes(filename, verbosity=0) - read the normal modes information from QE calculations (.dyn files)
 
## regexlib.py - definitions of some common regular expressions


