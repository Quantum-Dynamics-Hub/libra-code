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

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <sstream>
#endif

#include "ctx_Control_Parameters.h"
#include "../Units.h"

/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace libio;
using namespace liblinalg;



namespace libcontext{

//-------------- Class methods implementation ------------------------


ctx_Control_Parameters::ctx_Control_Parameters(){

  //----------------- All simulation parameters and flags (set to default values) -------------------
  // Convert everything to internal units (mostly atomic)

  // <calculations>
  std::string hamiltonian = "eht";          add("calculations.hamiltonian", hamiltonian);
  std::string spin_method = "unrestricted"; add("calculations.spin_method", spin_method);
  int DF = 0;                               add("calculations.DF", DF);
  // </calculations>

  // <guess_options>
  std::string guess_type = "sad";           add("guess_options.guess_type", guess_type);
  // </guess_options>

  // <scf_options>
  std::string scf_algo = "oda";             add("scf_options.scf_algo", scf_algo);
  int use_disk = 0;                         add("scf_options.use_disk", use_disk);
  int use_rosh = 0;                         add("scf_options.use_rosh", use_rosh);
  int do_annihilate = 0;                    add("scf_options.do_annihilate", do_annihilate);
  int pop_opt = 0;                          add("scf_options.pop_opt", pop_opt);
  int use_diis = 0;                         add("scf_options.use_diis", use_diis);
  int diis_max = 3;                         add("scf_options.diis_max", diis_max);
  int diis_start_iter = 0;                  add("scf_options.diis_start_iter", diis_start_iter);
  int use_level_shift = 0;                  add("scf_options.use_level_shift", use_level_shift);
  int shift_magnitude = 2.5;                add("scf_options.shift_magnitude", shift_magnitude);
  int use_damping = 0;                      add("scf_options.use_damping", use_damping);
  int damping_start = 3;                    add("scf_options.damping_start", damping_start);
  double damping_const = 0.05;              add("scf_options.damping_const", damping_const);
  double etol = 1e-6;                       add("scf_options.etol", etol);
  double den_tol = 1e-4;                    add("scf_options.den_tol", den_tol);
  int Niter = 300;                          add("scf_options.Niter", Niter);
  double degen_tol = 0.2;                   add("scf_options.degen_tol", degen_tol);
  // </scf_options>
  
  // <hamiltonian_options>
  std::string parameters = "none";          add("hamiltonian_options.parameters", parameters);
  std::string eht_params_format = "eht+0";  add("hamiltonian_options.eht_params_format", eht_params_format);
  int eht_formula = 1;                      add("hamiltonian_options.eht_formula", eht_formula);
  int eht_sce_formula = 0;                  add("hamiltonian_options.eht_sce_formula", eht_sce_formula);
  int eht_fock_opt = 1;                     add("hamiltonian_options.eht_fock_opt", eht_fock_opt);
  int eht_electrostatics = 0;               add("hamiltonian_options.eht_electrostatics", eht_electrostatics);
  // </hamiltonian_options>

  // <md_options>
  double md_dt = 1.0 * FS;                  add("md_options.md_dt", md_dt);
  int md_nsteps = 10;                       add("md_options.md_nsteps", md_nsteps);
  // </md_options>

  // <opt_options>
  double opt_dt = 1.0 * FS;                 add("opt_options.opt_dt", opt_dt);
  int opt_nsteps = 10;                      add("opt_options.opt_nsteps", opt_nsteps);
  // </opt_options>

  // <multipole_options>
  int compute_dipole = 1;                   add("multipole_options.compute_dipole", compute_dipole);
  // </multipole_options>

  // <dos_options>
  int compute_dos = 0;                      add("dos_options.compute_dos", compute_dos);
  std::string dos_opt = "dens";             add("dos_options.dos_opt", dos_opt);
  std::string dos_prefix = "dos/";          add("dos_options.dos_prefix", dos_prefix);
  // </dos_options>

  // <charge_density_options>
  int compute_charge_density = 0;           add("charge_density_options.compute_charge_density", compute_charge_density);
  int nx_grid = 40;                         add("charge_density_options.nx_grid", nx_grid);
  int ny_grid = 40;                         add("charge_density_options.ny_grid", ny_grid);
  int nz_grid = 40;                         add("charge_density_options.nz_grid", nz_grid);
  std::string charge_density_prefix = "char_dens/";     add("charge_density_options.charge_density_prefix", charge_density_prefix);
  vector<int> orbs(1,0);                    add("charge_density_options.orbs", orbs);
  // </charge_density_options>


  // <nac_options>
  std::string nac_md_trajectory_filename = "md_trajectory.xyz";  add("nac_options.nac_md_trajectory_filename", nac_md_trajectory_filename);
  std::string nac_prefix = "/res/Ham_";     add("nac_options.nac_prefix", nac_prefix);
  int nac_min_frame = 1;                    add("nac_options.nac_min_frame", nac_min_frame);
  int nac_max_frame = 5;                    add("nac_options.nac_max_frame", nac_max_frame);
  vector<int> nac_min_orbs(1,0);            add("nac_options.nac_min_orbs", nac_min_orbs);
  vector<int> nac_max_orbs(1,1);            add("nac_options.nac_max_orbs", nac_max_orbs);
  double nac_dt = 1.0;                      add("nac_options.nac_dt", nac_dt);
  int nac_opt = 0;                          add("nac_options.nac_opt", nac_opt);
  // </nac_options>

  // <scan_options>
  int scan_mov_at = 1;                      add("scan_options.scan_mov_at", scan_mov_at);
  int scan_ref_at = 0;                      add("scan_options.scan_ref_at", scan_ref_at);
  VECTOR scan_dir(1.0,0.0,0.0);             add("scan_options.scan_dir", scan_dir);
  double scan_dxmin = 0.0;                  add("scan_options.scan_dxmin", scan_dxmin);
  double scan_dxmax = 2.0;                  add("scan_options.scan_dxmax", scan_dxmax);
  double scan_dx = 0.25;                    add("scan_options.scan_dx", scan_dx);
  // </scan_options>

  // <excitation_options>
  int compute_excitations = 0;              add("excitation_options.compute_excitations", compute_excitations);  
  int num_excitations = 1;                  add("excitation_options.num_excitations", num_excitations);  
  std::string excitations_opt = "scf";      add("excitation_options.excitations_opt", excitations_opt);  
  double spectral_width = 0.1;              add("excitation_options.spectral_width", spectral_width);  
//  if(excitations.size()>0){ excitations.clear(); }
//  excitations.push_back(excitation(0,1,0,1));  // ground state
  // </excitation_options>

  // <unit_cell>
  VECTOR t1;                                add("unit_cell.t1", t1);  
  VECTOR t2;                                add("unit_cell.t2", t2);  
  VECTOR t3;                                add("unit_cell.t3", t3);  
  int x_period = 0;                         add("unit_cell.x_period", x_period);  
  int y_period = 0;                         add("unit_cell.y_period", y_period);  
  int z_period = 0;                         add("unit_cell.z_period", z_period);  
  // </unit_cell>

  // <coordinates>
  int Natoms = 0;                           add("coordinates.Natoms", Natoms);  
  double charge = 0.0;                      add("coordinates.charge", charge);  
  int spin = 1;                             add("coordinates.spin", spin);  
  std::string coordinates = "Cartesian";    add("coordinates.coordinates", coordinates);  
  // </coordinates>
 
  // <fragments>
  int num_frags = 1;                        add("fragments.num_frags", num_frags);

  for(int i=0;i<num_frags;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string out;    
    (ss << i);  ss >> out;

    int frag_size = 1;                        add("fragments.fragment"+out+".frag_size", frag_size);   
    std::string frag_name = "frag"+out;       add("fragments.fragment"+out+".frag_name", frag_name);       
    double frag_charge = 0.0;                 add("fragments.fragment"+out+".frag_charge", frag_charge); 
    vector<int> frag_orbitals(1,0);           add("fragments.fragment"+out+".frag_orbitals", frag_orbitals); 
    vector<int> frag_atoms(1,0);              add("fragments.fragment"+out+".frag_atoms", frag_atoms); 

  }// for i
  // </fragments>


}


//------------------ Export -------------------------

void export_ctx_Control_Parameters_objects(){

  class_<ctx_Control_Parameters, bases<Context> >("ctx_Control_Parameters",init<>())
  ;

}



}// namespace libcontext
}// liblibra


