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
  \file Control_Parameters.cpp
  \brief The file implement the class that stores the parameters controlling the calculations.     
*/

#include "Control_Parameters.h"

/// liblibra namespace
namespace liblibra{

/// libcontrol_parameters namespace
namespace libcontrol_parameters{



Control_Parameters::Control_Parameters(){
/**  
  Default constructor - sets the parameters to the default values:
*/

  //----------------- All simulation parameters and flags (set to default values) -------------------
  // Convert everything to internal units (mostly atomic)

  // <calculations>
  runtype = "scf";       /// runtype = "scf"
  hamiltonian = "eht";    /// hamiltonian = "eht"
  spin_method = "unrestricted";  /// spin_method = "unrestricted"
  DF = 0;                /// DF = 0  - no extra output by default, use it only for small systems and for benchmarking purposes
  // </calculations>

  // <guess_options>
  guess_type = "sad";  /// guess_type = "sad"
  // </guess_options>

  // <scf_options>
  scf_algo = "none";     /// scf_algo = "none" - This is the most robust option
  use_disk = 0;          /// use_disk = 0 
  use_rosh = 0;          /// use_rosh = 0 
  do_annihilate = 0;     /// do_annihilate = 0 -  do not do spin annihilation by default
  pop_opt = 0;           /// pop_opt = 0 - integer occupations

  use_diis = 0;          /// use_diis = 0
  diis_max = 3;          /// diis_max = 3
  diis_start_iter = 0;   /// diis_start_iter = 0

  use_level_shift = 0;   /// use_level_shift = 0
  shift_magnitude = 2.5; /// shift_magnitude = 2.5

  use_damping = 0;       /// use_damping = 0 
  damping_start = 3;     /// damping_start = 3 -  3-rd iteration will start damping
  damping_const = 0.05;  /// damping_const = 0.05 

  etol = 1e-6;           /// etol = 1e-6
  den_tol = 1e-4;        /// den_tol = 1e-4
  Niter = 300;           /// Niter = 300

  degen_tol = 0.2;       /// degen_tol = 0.2
  // </scf_options>
  
  // <hamiltonian_options>
  parameters = "none";   ///parameters = "none"
  // For EHT
  eht_params_format = "eht+0";  /// eht_params_format = "eht+0" default format for EHT parameters
  eht_formula     = 1;          /// eht_formula     = 1 - weighted formula
  eht_sce_formula = 0;          /// eht_sce_formula = 0 - no self-consistent electrostatics by default
  eht_fock_opt    = 1;          /// eht_fock_opt    = 1 - need self-consistency correction, if SC-EHT is used
  eht_electrostatics = 0;       /// eht_electrostatics = 0 -  no additional electrostatic effects
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
  compute_charge_density = 0;   
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
  if(excitations.size()>0){ excitations.clear(); }
  excitations.push_back(excitation(0,1,0,1));  // ground state
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
  

}


void get_parameters_from_file(std::string filename, Control_Parameters& prms){
/**
  Read the control parameters into the Control_Parameters object from an input file

  \param[in] filename The name of the input file
  \param[in,out] prms The object with control parameters
*/

  std::string st;
  vector< vector<std::string> > file;


  //---------------------------- Reading input file --------------------------------------

  cout<<"Reading input file = "<<filename<<endl;
  ifstream in(filename.c_str(), ios::in);
  if(in.is_open()){
    while(!in.eof()){
      getline(in,st); 

      vector<std::string> line;
      stringstream ss(st,stringstream::in|stringstream::out);
      while(ss>>st){ line.push_back(st);}

      file.push_back(line);

    }// while
  }else{ cout<<"Error: Can not open file\n";}
  in.close();


  cout<<"Echo tokenized file content\n";
  for(int i=0;i<file.size();i++){
    for(int j=0;j<file[i].size();j++){
      cout<<file[i][j]<<"  ";
    }
    cout<<endl;
  }// i


  //------------------ Now read in all parameters from the input files -----------------------------
  int f_sz = file.size();

  for(int i=0;i<f_sz;i++){   // line
    if(file[i].size()>0){
 
      if(file[i][0]=="<calculation>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</calculation>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){  
            if(file[i1][0]=="runtype"){  prms.runtype = file[i1][2]; } 
            else if(file[i1][0]=="hamiltonian"){  prms.hamiltonian = file[i1][2];  } 
            else if(file[i1][0]=="spin_method"){  prms.spin_method = file[i1][2];  }
            else if(file[i1][0]=="DF"){  prms.DF = atoi(file[i1][2].c_str());  }
          }
        }// for i1


        i = end_i + 1;
      }// <calculation>

      else if(file[i][0]=="<guess_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</guess_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){ 
            if(file[i1][0]=="guess_type"){  prms.guess_type = file[i1][2];  } 
          }
        }// for i1

        i = end_i + 1;
      }// <guess_options>


      else if(file[i][0]=="<hamiltonian>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</hamiltonian>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){ 
            // General
            if(file[i1][0]=="parameters"){  prms.parameters = file[i1][2];  } 

            // EHT-specific
            else if(file[i1][0]=="eht_params_format"){  prms.eht_params_format = file[i1][2];  } 
            else if(file[i1][0]=="eht_formula"){  prms.eht_formula = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="eht_sce_formula"){  prms.eht_sce_formula = atoi(file[i1][2].c_str());  }            
            else if(file[i1][0]=="eht_fock_opt"){  prms.eht_fock_opt = atoi(file[i1][2].c_str());  }            
            else if(file[i1][0]=="eht_electrostatics"){  prms.eht_electrostatics = atoi(file[i1][2].c_str());  }            
          }
        }// for i1

        i = end_i + 1;
      }// <hamiltonian>



      else if(file[i][0]=="<scf_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</scf_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){  
            if(file[i1][0]=="scf_algo"){  prms.scf_algo = file[i1][2].c_str();   } 
            else if(file[i1][0]=="use_disk"){ prms.use_disk = atoi(file[i1][2].c_str());   } 
            else if(file[i1][0]=="use_rosh"){ prms.use_rosh = atoi(file[i1][2].c_str());   } 
            else if(file[i1][0]=="do_annihilate"){ prms.do_annihilate = atoi(file[i1][2].c_str());   } 
            else if(file[i1][0]=="pop_opt"){  prms.pop_opt = atoi(file[i1][2].c_str());   } 
            else if(file[i1][0]=="use_diis"){  prms.use_diis = atoi(file[i1][2].c_str());   } 
            else if(file[i1][0]=="diis_max"){  prms.diis_max = atoi(file[i1][2].c_str());   } 
            else if(file[i1][0]=="diis_start_iter"){  prms.diis_start_iter = atoi(file[i1][2].c_str());   } 
            else if(file[i1][0]=="use_level_shift"){  prms.use_level_shift = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="shift_magnitude"){  prms.shift_magnitude = atof(file[i1][2].c_str());  } 
            else if(file[i1][0]=="use_damping"){  prms.use_damping = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="damping_start"){  prms.damping_start = atoi(file[i1][2].c_str());  }
            else if(file[i1][0]=="damping_const"){  prms.damping_const = atof(file[i1][2].c_str());  }
            else if(file[i1][0]=="etol"){  prms.etol = atof(file[i1][2].c_str());  }
            else if(file[i1][0]=="den_tol"){  prms.den_tol = atof(file[i1][2].c_str());  }
            else if(file[i1][0]=="Niter"){  prms.Niter = atoi(file[i1][2].c_str());  }
            else if(file[i1][0]=="degen_tol"){  prms.degen_tol = atof(file[i1][2].c_str());  }
          }
        }// for i1

        i = end_i + 1;
      }// <scf_options>


      else if(file[i][0]=="<md_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</md_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){ 
            if(file[i1][0]=="dt"){  prms.md_dt = atof(file[i1][2].c_str()) * FS;  } 
            else if(file[i1][0]=="nsteps"){  prms.md_nsteps = atoi(file[i1][2].c_str());  }
          }
        }// for i1

        i = end_i + 1;
      }// <md_options>



      else if(file[i][0]=="<opt_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</opt_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){  
            if(file[i1][0]=="dt"){  prms.opt_dt = atof(file[i1][2].c_str()) * FS;  } 
            else if(file[i1][0]=="nsteps"){  prms.opt_nsteps = atoi(file[i1][2].c_str());  } 
          }
        }// for i1

        i = end_i + 1;
      }// <opt_options>


      else if(file[i][0]=="<multipole_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</multipole_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){ 
            if(file[i1][0]=="compute_dipole"){  prms.compute_dipole = atoi(file[i1][2].c_str());  } 

          }
        }// for i1

        i = end_i + 1;
      }// <multipole_options>


      else if(file[i][0]=="<dos_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</dos_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){ 
            if(file[i1][0]=="compute_dos"){  prms.compute_dos = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="dos_opt"){  prms.dos_opt = file[i1][2];  } 
            else if(file[i1][0]=="dos_prefix"){  prms.dos_prefix = file[i1][2];  } 

          }
        }// for i1

        i = end_i + 1;
      }// <dos_options>


      else if(file[i][0]=="<charge_density_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</charge_density_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){  
            if(file[i1][0]=="compute_charge_density"){  prms.compute_charge_density = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="charge_density_prefix"){  prms.charge_density_prefix = file[i1][2];  } 
            else if(file[i1][0]=="nx_grid"){  prms.nx_grid = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="ny_grid"){  prms.ny_grid = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="nz_grid"){  prms.nz_grid = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="orbs"){  
              prms.orbs = vector<int>(file[i1].size()-2,0); 
              for(int i2=2;i2<file[i1].size();i2++){ prms.orbs[i2-2] = atoi(file[i1][i2].c_str()); }
            } 
          }
        }// for i1

        i = end_i + 1;
      }// <charge_density_options>



      else if(file[i][0]=="<nac_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</nac_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){  
            if(file[i1][0]=="nac_prefix"){  prms.nac_prefix = file[i1][2];  } 
            else if(file[i1][0]=="nac_md_trajectory_filename"){ prms.nac_md_trajectory_filename = file[i1][2]; }
            else if(file[i1][0]=="nac_min_frame"){  prms.nac_min_frame = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="nac_max_frame"){  prms.nac_max_frame = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="nac_min_orbs"){ 
              int nfr = atoi(file[i1][1].c_str());
              prms.nac_min_orbs = vector<int>(nfr,0); 
              for(int i2=0;i2<nfr;i2++){ prms.nac_min_orbs[i2] = atoi(file[i1][2+i2].c_str()); }
            } 
            else if(file[i1][0]=="nac_max_orbs"){ 
              int nfr = atoi(file[i1][1].c_str());
              prms.nac_max_orbs = vector<int>(nfr,0); 
              for(int i2=0;i2<nfr;i2++){ prms.nac_max_orbs[i2] = atoi(file[i1][2+i2].c_str()); }
            } 
            else if(file[i1][0]=="nac_dt"){  prms.nac_dt = atof(file[i1][2].c_str());  }
            else if(file[i1][0]=="nac_opt"){  prms.nac_opt = atoi(file[i1][2].c_str());  }
          }
        }// for i1

        i = end_i + 1;
      }// <nac_options>



      else if(file[i][0]=="<scan_options>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</scan_options>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){  
            if(file[i1][0]=="scan_mov_at"){  prms.scan_mov_at = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="scan_ref_at"){  prms.scan_ref_at = atoi(file[i1][2].c_str());  }
            else if(file[i1][0]=="scan_dxmin"){  prms.scan_dxmin = atof(file[i1][2].c_str());  } 
            else if(file[i1][0]=="scan_dxmax"){  prms.scan_dxmax = atof(file[i1][2].c_str());  } 
            else if(file[i1][0]=="scan_dx"){  prms.scan_dx = atof(file[i1][2].c_str());  } 

          }// size > 2

          if(file[i1].size()>4){  
            if(file[i1][0]=="scan_dir"){  
              prms.scan_dir.x  = atof(file[i1][1].c_str()); 
              prms.scan_dir.y  = atof(file[i1][2].c_str()); 
              prms.scan_dir.z  = atof(file[i1][3].c_str()); 
            }
          }// size > 4

        }// for i1

        i = end_i + 1;
      }// <scan_options>


      else if(file[i][0]=="<excitations>"){
        cout<<"Internalizing <excitations> info:\n";
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</excitations>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        for(int i1=i+1;i1<end_i;i1++){
          if(file[i1].size()>2){  

            if(file[i1][0]=="compute_excitations"){  prms.compute_excitations = atoi(file[i1][2].c_str());  } 
            else if(file[i1][0]=="excitations_opt"){  prms.excitations_opt = file[i1][2];  } 
            else if(file[i1][0]=="spectral_width"){  prms.spectral_width = atof(file[i1][2].c_str());  } 
            else if(file[i1][0]=="num_excitations"){  
              // Number of excitations
              int tmp_sz = atoi(file[i1][2].c_str());

              prms.num_excitations = 0;
              if(prms.excitations.size()>0){ prms.excitations.clear(); }


              // Actual excitations
              for(int n=0;n<tmp_sz;n++){
                int ex_size = atoi(file[i1+1+n][0].c_str());
            
                if(ex_size!=1){  cout<<"Warning: Only single excitations are currently implemented. Skipping this excitation\n"; }
                else{

                  cout<<"Excitation #"<<prms.num_excitations;
            
                  prms.num_excitations++;
                  //------------ From ---------------------
                  std::string _from = file[i1+1+n][1];
                  int len = _from.size();
                  std::string x; x="";
                  for(int l=0;l<(len-1);l++){ x = x + _from[l]; } 
                  
                  int _f_o = atoi(x.c_str());
                  cout<<_from[len-1]<<endl;
                  int _f_s = (_from[len-1]=='A')?1:-1;
                  cout<<"_f_o = "<<_f_o<<endl;
                  cout<<"_f_s = "<<_f_s<<endl;
            
                  //------------ To ---------------------
                  std::string _to = file[i1+1+n][3];
                  len = _to.size();
                  x="";
                  for(int l=0;l<(len-1);l++){ x = x + _to[l]; } 
                  
                  int _t_o = atoi(x.c_str());
                  int _t_s = (_to[len-1]=='A')?1:-1;
                  cout<<"_t_o = "<<_t_o<<endl;
                  cout<<"_t_s = "<<_t_s<<endl;

            
                  prms.excitations.push_back( excitation(_f_o,_f_s,_t_o,_t_s) );
            
            
                }// else
            
              }// for n
            
              cout<<"Number of excitations = "<<prms.num_excitations<<endl;


            }// num_excitations

          }// > 2
        }// for i1

     
        i = end_i + 1;
      }// <configuration>



      else if(file[i][0]=="<unit_cell>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</unit_cell>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // now analyze all lines in between
        prms.t1 = VECTOR( atof(file[i+2][0].c_str()),atof(file[i+2][1].c_str()),atof(file[i+2][2].c_str()) );
        prms.t2 = VECTOR( atof(file[i+3][0].c_str()),atof(file[i+3][1].c_str()),atof(file[i+3][2].c_str()) );
        prms.t3 = VECTOR( atof(file[i+4][0].c_str()),atof(file[i+4][1].c_str()),atof(file[i+4][2].c_str()) );
        prms.t1 *= Angst; prms.t2 *= Angst; prms.t3 *= Angst;
//        cell = Cell(t1,t2,t3);

        // Periodicity  - flags to check if the system is periodic in given direction
        prms.x_period = atoi(file[i+2][3].c_str());
        prms.y_period = atoi(file[i+3][3].c_str());
        prms.z_period = atoi(file[i+4][3].c_str());

        i = end_i + 1;
      }// <unit_cell>


      else if(file[i][0]=="<fragments>"){
        // search for end of this group
        int end_i = i;
        for(int i1=i+1;i1<f_sz;i1++){
          if(file[i1].size()>0){  if(file[i1][0]=="</fragments>"){ end_i = i1; break; }    }// non-empty line
        }// for i1


        // Number of atoms and coordinate system
        int nfrags = atoi(file[i+1][0].c_str());
  
        // Actual information about each fragment 
        for(int n=0;n<nfrags;n++){

          prms.frag_name.push_back(file[i+3+n][1]);
          prms.frag_charge.push_back(atof(file[i+3+n][2].c_str()));
          int sz = atoi(file[i+3+n][3].c_str());  prms.frag_size.push_back(sz);

          vector<int> frag;
          for(int i1=0;i1<sz;i1++){
            int at_indx = atoi(file[i+3+n][4+i1].c_str()) - 1; // !!! Because input contains atomic numbers (so min is 1, not 0)
            frag.push_back(at_indx);
          }
          prms.fragments.push_back(frag);           

        }// for n

        i = end_i + 1;
      }// <fragments>


    }// non-empty line
  }// for i - lines in the file


  //===================== Some analysis of what is read from input/used as default ===================
  if(prms.scf_algo=="oda"||prms.scf_algo=="diis_fock"||prms.scf_algo=="diis_dm"||prms.scf_algo=="none"){
  }else{
    cout<<"Error: prms.scf_algo = "<<prms.scf_algo<<" is unknown or not registered\n";
    cout<<" possible values are:\n";
    cout<<" none      - for standard SCF algorithm (default)\n";
    cout<<" oda       - for optimal damping algorithm\n";
    cout<<" diis_fock - for DIIS with Fock matrix mixing/extrapolation\n";
    cout<<" diis_dm   - for DIIS with density matrix mixing/extrapolation\n";
    exit(0);
  }

  if(prms.compute_excitations==1 && prms.compute_dipole!=1){
    cout<<"To compute excitations (spectra) dipole moments must be computed. Set compute_dipole to value 1\n";
    exit(0);
  }

  if(prms.hamiltonian=="eht"){
    if(prms.eht_sce_formula==1){
      if(prms.eht_params_format=="eht+0"){
        cout<<"Error: eht_sce_formula = 1 requires eht_params_format to be more general than eht+0\n"; exit(0);
      }       
    }// SC-EHT

/*
    if(prms.eht_sce_formula==0 && prms.guess_type!="core"){
      cout<<"Non self-consistent EHT must be used with guess_type=\"core\"\n"; exit(0);
    }

    if(prms.eht_sce_formula>0 && prms.guess_type!="sad"){
      cout<<"Self-consistent EHT must be used with guess_type=\"sad\"\n"; exit(0);
    }
*/

  }// eht



}


}// namespace libhcontrol_parameters
}// namespace liblibra






