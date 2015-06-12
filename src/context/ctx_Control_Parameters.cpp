#include "ctx_Control_Parameters.h"



namespace libcontext{

//-------------- Class methods implementation ------------------------


ctx_Control_Parameters::ctx_Control_Parameters(){

  //----------------- All simulation parameters and flags (set to default values) -------------------
  // Convert everything to internal units (mostly atomic)

  // <calculations>
  std::string runtype = "scf";       add(".runtype", runtype);
  std::string hamiltonian = "eht";   add(".hamiltonian", hamiltonian);
  std::string spin_method = "unrestricted"; add(".spin_method", spin_method);
  int DF = 0;               add(".DF", DF);
  // </calculations>

  // <guess_options>
  std::string guess_type = "sad";    add(".guess_type", guess_type);
  // </guess_options>

  vector<int> nac_min_orbs(3); nac_min_orbs[0] = 1; nac_min_orbs[1] = 2; nac_min_orbs[2] = 3;  add(".nac_min_orbs",nac_min_orbs);







  // <scf_options>
  scf_algo = "oda";      // This is the most robust option
  use_disk = 0;          // 0 - do not use  disk (faster);  1 - use disk (less memory required)
  use_rosh = 0;          // use open shell (not to use restricted open-shell)
  do_annihilate = 0;     // do not do spin annihilation by default
  pop_opt = 0;           // if = 0 - integer occupations,  if = 1 - fractional occupations based on Fermi distribution

  use_diis = 0;
  diis_max = 3;
  diis_start_iter = 0;

  use_level_shift = 0;
  shift_magnitude = 2.5;

  use_damping = 0;       // if = 0 - ODA is used (default, robust), if = 1 - damping with constant mixing value is used
  damping_start = 3;     // 3-rd iteration will start dampbing
  damping_const = 0.05;  // the smaller the constant, more likely the SCF will convege, but it may be slower than for larger constant

  etol = 1e-6;       
  den_tol = 1e-4;   
  Niter = 300;       

  degen_tol = 0.2;
  // </scf_options>
  
  // <hamiltonian_options>
  parameters = "none";
  // For EHT
  eht_params_format = "eht+0";     // default format for EHT parameters
  eht_formula     = 1;             // weighted formula
  eht_sce_formula = 0;             // no self-consisten electrostatics by default
  eht_fock_opt    = 1;             // need self-consistency correction, if SC-EHT is used
  eht_electrostatics = 0;          // no additional electrostatic effects
  // </hamiltonian_options>


  // <md_options>
  md_dt = 1.0 * FS;             
  md_nsteps = 10;          
  // </md_options>

  // <opt_options>
  opt_dt = 1.0 * FS;  
  opt_nsteps = 10;         
  // </opt_options>

  // <multipole_options>
  compute_dipole = 1;             // do compute dipole moment
  // </multipole_options>

  // <dos_options>
  compute_dos = 0;                // do not compute DOS by default
  dos_opt = "dens";               // DOS computations based on density matrix
  dos_prefix = "dos/";  
  // </dos_options>

  // <charge_density_options>
  nx_grid = ny_grid = nz_grid = 40;
  charge_density_prefix = "char_dens/";
  orbs = vector<int>(1,0);
  // </charge_density_options>


  // <nac_options>
  nac_md_trajectory_filename = "md_trajectory.xyz";
  nac_prefix = "/res/Ham_"; 
  nac_min_frame = 1;
  nac_max_frame = 5;        
  nac_min_orbs = vector<int>(1,0); // one fragment, 0-eth orbital        
  nac_max_orbs = vector<int>(1,1); // one fragment, 1-st orbital
  nac_dt = 1.0;             // convention is to compute NACs in units of [Ha/fs], so don't transform nac_dt to a.u. of time
  nac_opt = 0;              // Tully-Hamess-Schiffer, 1 = add non-orthogonality correction - non-Hermitian
  // </nac_options>

  // <scan_options>
  scan_mov_at = 1;               
  scan_ref_at = 0;               
  scan_dir = VECTOR(1.0,0.0,0.0);
  scan_dxmin = 0.0;       
  scan_dxmax = 2.0;       
  scan_dx = 0.25;         
  // </scan_options>

  // <excitations>
  // Default - is just a ground state configuration
  compute_excitations = 0; // 0 = "no" by default, 1 = "yes"
  num_excitations = 1;
  excitations_opt = "scf";
  spectral_width = 0.1; // eV 
//  if(excitations.size()>0){ excitations.clear(); }
//  excitations.push_back(excitation(0,1,0,1));  // ground state
  // </excitations>


  // <unit_cell>
  t1 = t2 = t3 = 0.0;
  x_period = 0;  
  y_period = 0;  
  z_period = 0;  
  // </unit_cell>

  // <coordinates>
  Natoms = 0;               
  charge = 0.0;
  spin = 1;
  coordinates = "Cartesian";
  // </coordinates>
 


//  save(ctx_pt,"ctx_Control_Parameters");

}


/*
void ctx_Control_Parameters::save(boost::property_tree::ptree& pt,std::string path){
  
  int st;

  libio::save(pt,path+".runtype",runtype);
  libio::save(pt,path+".hamiltonian",hamiltonian);


}

void ctx_Control_Parameters::load(boost::property_tree::ptree& pt,std::string path,int& status){

  int st;

  libio::load(pt,path+".runtype",runtype, st);
  libio::load(pt,path+".hamiltonian",hamiltonian, st);

}
*/

//------------------ Export -------------------------

void export_ctx_Control_Parameters_objects(){

  void (Context::*expt_add_v1)(std::string varname, int varval) = &Context::add;
  void (Context::*expt_add_v2)(std::string varname, vector<int> varval) = &Context::add;
  void (Context::*expt_add_v3)(std::string varname, std::string varval) = &Context::add;
  void (Context::*expt_add_v4)(std::string varname, vector<std::string> varval) = &Context::add;
  void (Context::*expt_add_v5)(std::string varname, double varval) = &Context::add;
  void (Context::*expt_add_v6)(std::string varname, vector<double> varval) = &Context::add;


  class_<ctx_Control_Parameters>("ctx_Control_Parameters",init<>())
//      .def(init<const ctx_Control_Parameters&>())
//      .def("__copy__", &generic__copy__<ctx_Control_Parameters>)
//      .def("__deepcopy__", &generic__deepcopy__<ctx_Control_Parameters>)

       // Members
//      .def_readwrite("runtype",&ctx_Control_Parameters::runtype)
//      .def_readwrite("hamiltonian",&ctx_Control_Parameters::hamiltonian)

      .def("add",expt_add_v1)
      .def("add",expt_add_v2)
      .def("add",expt_add_v3)
      .def("add",expt_add_v4)
      .def("add",expt_add_v5)
      .def("add",expt_add_v6)


      .def("save_xml",&ctx_Control_Parameters::save_xml)
      .def("load_xml",&ctx_Control_Parameters::load_xml)
  ;


}



}// namespace libcontext
