# Description of the modules

## common_utils.py  - Auxiliary general-purpose functions
   * check_input(params, default_params, critical_params) - checks the presenece of required parameters
   * get_matrix(nrows, ncols, filename_re, filename_im, act_sp) - reads the files into a matrix
   * orbs2spinorbs(s) - conversion from the older PYXAID style to the new style
   * find_maxima(s, logname) - finds the maxima in the data series
   * flt_stat(X) - computes statistical properties on the scalar data series
   * mat_stat(X) - computes statistical properties on the MATRIX data series
   * cmat_stat(X) - computes statistical properties on the CMATRIX data series
   * printout(t, pops, Hvib, outfile) - one of the versions to print data out
   * show_matrix_splot(X, filename) - print out of a matrix in a gnuplot format for surface plotting
   * add_printout(i, pop, filename) - one of the versions for printing out the TSH information

## compute_hprime.py - module to compute the transition dipole moment matrices
   * compute_hprime_dia(es, info, filename) 
   * hprime_py(es, info, filename)

## compute_properties.py - auxiliary module to compute some properties
   * compute_properties_dia_gamma(params, es_curr, es_next, curr_index):

## decoherence_times.py - to compute decoherence times from the data series
   * energy_gaps(Hvib) - compute the energy gaps along a single trajectory
   * energy_gaps_ave(Hvib, itimes, nsteps) - compute the energy gaps averaged over several trajectories
   * decoherence_times(Hvib, verbosity=0) - compute decoherence times based on single trajectory data
   * decoherence_times(Hvib, itimes, nsteps, verbosity=0) - compute decoherence times based on several trajectories data

## excitation_spectrum.py - to generate Hvib matrices averaged over trajectories
   * calculate(energy_prefix,energy_suffix,dip_prefix,dip_suffix,isnap,fsnap,opt,scl1,scl2,outfile,HOMO,minE,maxE,dE)
   * ham_map(prefix, isnap, fsnap, suffix, opt, scl, outfile)

## influence_spectrum.py - computes the ACF of Hamiltonian matrix elements and its FT 
   * compute(X, a, b, dt, Nfreqs, filename, logname, dw=1.0, wspan=3000.0)

## lz.py - Landau-Zener theory as applied to NA-MD (experimental module)
   * get_data(params) - read in the vibronic Hamiltonian files 
   * Belyaev_Lebedev(Hvib, dt) - compute probabilities of NA transitions along a trajectory
   * run(params) - execute the NA-MD calculations based on the LZ hopping probabilities

## mapping.py - to map the 1-electron properties onto N-electron ones
   * energy_arb(SD, e) - compute energy of a Slater determinant (simplified)
   * energy_mat_arb(SD, e, dE) - compute matrix of the energies of the SD basis, with corrections
   * sd2indx(inp,nbasis) - mapping of spin-orbitals on the orbitals
   * ovlp_arb(SD1, SD2, S) - compute the overlap of two SDs
   * ovlp_mat_arb(SD1, SD2, S) - compute a matrix of overlaps of the SD basis


## qsh.py - generates the QSH files for longer trajectory
   * compute_freqs(nstates, H_vib_re, H_vib_im, dt, Nfreqs, filename, logname, dw, wspan) - compute
     the frequencies of the matrix elements in a time-series
   * compute_qs_Hvib(Nfreqs, freqs, t, nstates, 
                 H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                 H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                 dev) - compute the quasi-stochastic Hamiltonian
   * run(params) - execute the QSH calculations


## step2.py - run QE calculations and computes the time-overlaps of the orbitals
   * run_qe(params, t, dirname0, dirname1) - execute the QE calculations and do the file management
   * read_info(params) - read the info about the QE calculation
   * read_all(params) - reade the QE calculation results - energies, orbitlas, etc.
   * read_wfc_grid(params, info0, info1) - read the wavefunction grid
   * run(params) - execute the calculations for a "step2"

## step3.py - to generate the Hvib matrices for given basis of many-electron states with needed corrections
   * get_data(params) - read in the information about the 1-electron orbitals (time-overlaps, overlaps, and energies)
   * apply_state_reordering(St, E, params) - state reordering
   * apply_phase_correction(St, params) - phase correction
   * sac_matrices(coeff) - an auxiliary function to convert the user input about SAC to the internal format
   * scale_H_vib(hvib, en_gap, dNAC, sc_nac_method) - scissor operator to the energy gap and the NACs
   * compute_Hvib(basis, St_ks, E_ks, dE, dt) - compute Hvib in the Slater determinants basis
   * run(params) - using the input of the 1-el KS  orbitals, compute the Hvib files in the specified SAC basis

## step4.py - to compute the NBRA-based NA-MD
   * get_Hvib(params) - read in the vibronic Hamiltonians along the trajectories from the files (different data sets)
   * traj_statistics(i, Coeff, istate, Hvib, itimes) - compute the averaged populations and energies from TSH data
   * printout(t, res, outfile) - print out the data to a file
   * run(params) - exectue the NA-MD for multiple trajectories 


## utils.py - various auxiliary/utility functions
   * get_value(params,key,default,typ) - check if the dictionary has the required keys/set default values
   * split_orbitals_energies(C, E) - needed for processing SOC output
   * merge_orbitals(Ca, Cb) - change the format of the PW orbitals storage
   * post_process(coeff, ene, issoc) - transform the format of orbitals
   * orthogonalize_orbitals(C) - orthogonalize orbitals, if needed
   * orthogonalize_orbitals2(Ca,Cb) - orthogonalization of the spinor wavefunction




# Roadmap of workflows

## PYXAID (old version)
   * run adiabatic MD - one time
   * step2.run()  - execute one time
   * step3.run()  - without corrections on, no spin-adaptation
   * step4.run()

## PYXAID2 
   * run adiabatic MD - potentially several trajectories
   * step2.run()  - execute several times, with available trajectories
   * step3.run()
   * step4.run()

## Quasistochastic Hamiltonian
   * run adiabatic MD - potentially several trajectories
   * step2.run()  - execute one/several times, with available trajectories
   * step3.run()  - execute or not
   * qsh.run()
   * step4.run()


