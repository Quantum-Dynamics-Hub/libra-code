/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Control_Parameters.h
  \brief The file describes the class that stores the parameters controlling the calculations. Also the 
  auxiliary classes and functions are defined
    
*/

#ifndef CONTROL_PARAMETERS_H
#define CONTROL_PARAMETERS_H

#include "../math_linalg/liblinalg.h"
#include "../common_types/libcommon_types.h"
#include "../Units.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libcommon_types;


/// libcontrol_parameters namespace
namespace libcontrol_parameters{





class Control_Parameters{
/**
  The class that stores the parameters controlling calculations. Note: not all prameters and options are presently used.
*/

public:
//---------- Members --------

  //----------------- All simulation parameters and flags (set to default values) -------------------
  // <calculations>
  std::string runtype;           ///< Calculation type. 
                                 ///< Possible options: "scf", "scan", "dos", "md", "opt", "nac"
                                 ///< Default: "scf"
  std::string hamiltonian;       ///< Hamiltonian type. 
                                 ///< Possible options: "eht", "indo", "cndo2", etc.
                                 ///< Default: "eht"
  std::string spin_method;       ///< Way the spin is treated. 
                                 ///< Possible options:  "restricted", "unrestricted"
                                 ///< Default: "unrestricted"
  int DF;                        ///< Debug flag: if set to 1 - will print a lot of information. Be carefull!!!
                                 ///< Default: 0
  // </calculations>

  // <guess>
  std::string guess_type;        ///< Define how to make guess orbitals
                                 ///< Possible options: "sad" (superposition of atomic densities), "core" (core Hamiltonian states)
                                 ///< Default: "sad"
  // </guess>

  // <scf_options>
  std::string scf_algo;          ///< Algorithm for SCF iterations. 
                                 ///< Possible options: "none", "oda"
                                 ///< Default: "none"
  int use_disk;                  ///< write temporary variables to disk instead of RAM - this can help reducing memory costs
                                 ///< Possible options: 0 - do not use  disk (faster);  1 - use disk (less memory required)
                                 ///< Default: 0
  int use_rosh;                  ///< use restricted open-shell
                                 ///< Possible options: 1 (use), 0 (do not use)
                                 ///< Default: 0
  int do_annihilate;             ///< Do spin annihilation at the last iteration
                                 ///< Possible options: 1 (do annihilation), 0 (don't do annihilation)
                                 ///< Default: 0
  int pop_opt;                   ///< Occupation scheme - How to populate energy levels: 
                                 ///< Possble options: 0 - integer occupations, 1 - fractional occupations based on Fermi distribution 
                                 ///< Default: 0
  int use_diis;                  ///< flag to turn on/off DIIS calculations (presently not affecting calculations)
                                 ///< Possible options: 0 - do not use DIIS; 1 - use DIIS
                                 ///< Default: 0
  int diis_max;                  ///< Dimension of DIIS matrix, if used
                                 ///< Possible values: 1, 2, 3, ...
                                 ///< Default: 3
  int diis_start_iter;           ///< Iteration after which DIIS will start
                                 ///< Possible values: 0, 1, 2, ...
                                 ///< Default: 0
  int use_level_shift;           ///< Flag to turn on/off level shifting (LS) (not yet implemented)
                                 ///< Possible options: 0 - do not use LS; 1 - use LS
                                 ///< Default: use_level_shift = 0
  double shift_magnitude;        ///< The magnitude of the energy level shifts, if used
                                 ///< Possible options: any numerical (real) value, a.u.
                                 ///< Default: 2.5
  int use_damping;               ///< Flag to turn on/off damping        
                                 ///< Possible options: 0 - do not use damping (if ODA is used, then electronic optimization step
                                 ///< will be varying in the magnitude); 1 - use damping (if ODA is used - the electronic step
                                 ///< magnitude will be fixed, leading to the standard density matrix mixing scheme)
                                 ///< If "scf_algo" is = "none", it will not affect the calculations.
                                 ///< Default: 0
  int damping_start;             ///< The number of (standard) iterations before damping is in effect
                                 ///< Possible options: 0, 1, 2, ...
                                 ///< Default: 3
  double damping_const;          ///< Parameter for damping, if used - this is the magnitude of electronic iteration step
                                 ///< the smaller the constant, more likely the SCF will convege, but it may be slower than for a larger constant
                                 ///< Possible opions: any numerical (real) value in the [0.0, 1.0] interval
                                 ///< Default: 0.05
  double etol;                   ///< Energy convergence criterium, [Ha] 
                                 ///< Possible options: anything > 0.0
                                 ///< Default: 1e-6
  double den_tol;                ///< Density convergence criterium 
                                 ///< Possible options: anything > 0.0
                                 ///< Default: 1e-4
  int Niter;                     ///< The maximal number of SCF iterations before SCF is considered not converged
                                 ///< Possible options: > 1
                                 ///< Default: 300
  double degen_tol;              ///< The amount of population difference between two levels, when one can say the two 
                                 ///< levels are degenerate
                                 ///< Possible options: anything in the interval [0.0, 1.0]
                                 ///< Default: 0.2
  // </scf_options>

  // <hamiltonian_options>
  std::string parameters;        ///< Name of the file that contains parameters for Hamiltonian
                                 ///< Default: "none"
  std::string eht_params_format; ///< Format of the file that contains parameters for EHT Hamiltonian
                                 ///< Possible options:
                                 ///< "eht+0"   :  minimal EHT format
                                 ///< "eht+n"   :  add n more parameters at the end of each line, the meaning of the parameters
                                 ///<              may be different, depending on which method is used
                                 ///< "eht+n+K" :  same as eht+n, but also use K_ij as parameters - so far this is most general
                                 ///< Default: "eht+0"
  int eht_formula;               ///< formula for EHT Hamiltonian:
                                 ///< Possible options:
                                 ///< 0 - unweighted 
                                 ///< 1 - weighted
                                 ///< 2 - Calzaferi correction
                                 ///< 3 - my method (testing!)
                                 ///< Default: 1

  int eht_sce_formula;           ///< how to treat self-consistent electrostatics for EHT,
                                 ///< Posible options:
                                 ///< 0 - no self-consistent electrostatics
                                 ///< 1 - total-charge dependent IP
                                 ///< 2 - orbital-reolved corrections
                                 ///< 3 - my addition of perametric exchange 
                                 ///< Default: 0

  int eht_fock_opt;              ///< how to treat EHT Hamiltonian (H_eht) - self-consistency correction:
                                 ///< Opssible options:
                                 ///< 0 - as a Fock matrix (F_eht = H_eht) - no correction to self-consistency is needed, but the actual energy functional is different
                                 ///< 1 - as a model Hamiltonian (F_eht = 2*H_eht - H_eht0) - correction for self-consistency is needed
                                 ///< Default: 1         

  int eht_electrostatics;        ///< how to describe additional electrostatic interactions
                                 ///< Possible options:
                                 ///< 0 - no additional field effect
                                 ///< 1 - include pairwise Coulombic effects via Mulliken charges
                                 ///< Default: 0
  // </hamiltonian_options>

  // <properties>
  int compute_vertical_ip;       
  int compute_vertical_ea; 
  // </properties>

  // <md_options>
  double md_dt;                  ///< integration time step for MD [input in fs, internally in a.u.]
  int md_nsteps;                 ///< number of time steps for MD
  // </md_options>

  // <opt_options>
  double opt_dt;                 ///< integration time step for optimization [input in fs, internally in a.u.]
  int opt_nsteps;                ///< number of time steps for optimization
  // </opt_options>

  // <multipole_options>
  int compute_dipole;            ///< flag to turn dipole moment calculations
  // </multipole_options>

  // <dos_options>
  int compute_dos;               ///< flag to turn DOS calculations on/off
  std::string dos_opt;           ///< option for DOS comutation: "dens" - based on density matrix, "wfc" - based on wavefunction
  std::string dos_prefix;        ///< Prefix for the files in which atomic-projected DOS will be written  
  // </dos_options>

  // <charge_density_options>
  int compute_charge_density;         ///< flag to turn computation of charge density on
  int nx_grid, ny_grid, nz_grid;      ///< Number of voxels along each direction
  std::string charge_density_prefix;  ///< Prefix for the files in which CUBE orbitals will be written  
  vector<int> orbs;
  // </charge_density_options>


  // <nac_options>
  std::string nac_md_trajectory_filename;  ///< Name of the file that contains coordinates of the MD trajectory (in xyz format)
  std::string nac_prefix;        ///< Prefix of all files into wich NACs for different time steps will be written
  int nac_min_frame;             ///< index of minimal MD snapshot to include in NAC
  int nac_max_frame;             ///< index of maximal MD snapshot to include in NAC
  vector<int> nac_min_orbs;      ///< indexes of minimal orbitals to include in NAC matrix - for each fragment
  vector<int> nac_max_orbs;      ///< indexes of maximal orbitals to include in NAC matrix - for each fragment
  double nac_dt;                 ///< time step separating MD point in precomputed trajectory [fs]
  int nac_opt;                   ///< how to compute NACs between non-orthonormal orbitals
  // </nac_options>

  // <scan_options>
  int scan_mov_at;               ///< index of the atom that will be moved
  int scan_ref_at;               ///< index of the atom that serves as the reference atom
  VECTOR scan_dir;               ///< direction of the scan
  double scan_dxmin;             ///< initial displacement along the scan direction w.r.t position of the reference atoms
  double scan_dxmax;             ///< final displacement along the scan direction w.r.t position of the reference atoms
  double scan_dx;                ///< increment of displacement
  // </scan_options>

  // <excitations>
  int compute_excitations;       ///< flag turning on/off the actual computation of excitations
  int num_excitations;           ///< number of excitations to consider
  std::string excitations_opt;   ///< option for how to compute excitation energies
  double spectral_width;         ///< parameter for spectra calculation - width of the smearing
  vector<excitation> excitations; 
  // </excitations>

  // <unit_cell>
  VECTOR t1,t2,t3;               ///< vectors of periodic translations in 3 directions
  int x_period;                  ///< if periodic along t1
  int y_period;                  ///< if periodic along t2
  int z_period;                  ///< if periodic along t3
  // </unit_cell>

  // <coordinates>
  int Natoms;                    ///< Number of atoms
  double charge;                 ///< total charge of the system (-number of excess electrons)
  int spin;                      ///< 1 = singlet, 2 = doublet, 3 = triplet, etc.
  std::string coordinates;       ///< Direct or Cartesian
                                 ///< actual coordinates are stored separately
  // </coordinates>

  // <fragments>
  vector< vector<int> > fragments;  ///< list of atomic indices for atoms that are included in given fragment
  vector< int >         frag_size;  ///< size of the fragments
  vector< std::string>  frag_name;  ///< names of the fragments
  vector< double >      frag_charge;///< charges of the fragments, sum must be equal to charge - see <coordinates> sections
  // </fragments>
  


//--------- Methods ----------
    Control_Parameters();
   ~Control_Parameters(){ ;; }
    Control_Parameters(const Control_Parameters& x){ *this = x; }

};


void get_parameters_from_file(std::string, Control_Parameters&);


}// namespace libhcontrol_parameters
}// namespace liblibra






#endif // CONTROL_PARAMETERS_H
