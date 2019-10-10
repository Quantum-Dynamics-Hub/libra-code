/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_control_params.cpp
  \brief The file implements the methods to setup control parameters for dynamics
*/

#include "dyn_control_params.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;

dyn_control_params::dyn_control_params(){

/**

  This function initializes the default values of control parameters

*/

  rep_tdse = 1;
  rep_ham = 0;
  rep_sh = 1;
  rep_lz = 0;
  tsh_method = 0;
  force_method = 1;
  nac_update_method = 1;
  rep_force = 1;
  use_boltz_factor = 0;
  Temperature = 300.0;
  do_reverse = 1;
  vel_rescale_opt = 0;
  dt = 41.0;
  do_phase_correction = 1;
  state_tracking_algo = 2;
  MK_alpha = 0.0;
  MK_verbosity = 0;

  entanglement_opt = 0;
  ETHD3_alpha = 1.0;
  ETHD3_beta = 1.0;

  decoherence_algo = -1; 
  decoh_rates = NULL;

}



void dyn_control_params::sanity_check(){

  if(state_tracking_algo==0 || state_tracking_algo==1 ||
     state_tracking_algo==2 || state_tracking_algo==3){ ; ; }
  else{
    std::cout<<"Error in dyn_control_params::sanity_check: state_tracking_algo = "
        <<state_tracking_algo<<" is not allowed. Exiting...\n";
    exit(0);
  }

}



void dyn_control_params::set_parameters(bp::dict params){
/**
  Extract the parameters from the input dictionary
*/

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);


    if(key=="rep_tdse") { rep_tdse = bp::extract<int>(params.values()[i]); }
    else if(key=="rep_ham") { rep_ham = bp::extract<int>(params.values()[i]);   }
    else if(key=="rep_sh") { rep_sh = bp::extract<int>(params.values()[i]);  }
    else if(key=="rep_lz") { rep_lz = bp::extract<int>(params.values()[i]);  }
    else if(key=="tsh_method") { tsh_method = bp::extract<int>(params.values()[i]);  }
    else if(key=="force_method") { force_method = bp::extract<int>(params.values()[i]);  }
    else if(key=="nac_update_method") { nac_update_method = bp::extract<int>(params.values()[i]);  }
    else if(key=="rep_force") { rep_force = bp::extract<int>(params.values()[i]);  }
    else if(key=="use_boltz_factor") { use_boltz_factor = bp::extract<int>(params.values()[i]);  }
    else if(key=="Temperature") { Temperature = bp::extract<double>(params.values()[i]);  }
    else if(key=="do_reverse") { do_reverse = bp::extract<int>(params.values()[i]);  }
    else if(key=="vel_rescale_opt") { vel_rescale_opt = bp::extract<int>(params.values()[i]);  }

    else if(key=="dt") { dt = bp::extract<double>(params.values()[i]);  }

    // Phase correction
    else if(key=="do_phase_correction") { do_phase_correction = bp::extract<int>(params.values()[i]);  }

    // State tracking options
    else if(key=="state_tracking_algo"){  state_tracking_algo = bp::extract<int>(params.values()[i]);  }
    else if(key=="MK_alpha") { MK_alpha = bp::extract<double>(params.values()[i]);  }
    else if(key=="MK_verbosity") { MK_verbosity = bp::extract<int>(params.values()[i]);  }

    // Trajectory coupling
    else if(key=="entanglement_opt"){ entanglement_opt = bp::extract<int>(params.values()[i]); }
    else if(key=="ETHD3_alpha") { ETHD3_alpha = bp::extract<double>(params.values()[i]);   }
    else if(key=="ETHD3_beta") { ETHD3_beta = bp::extract<double>(params.values()[i]);   }


    else if(key=="decoherence_algo"){ decoherence_algo = bp::extract<int>(params.values()[i]); }
    else if(key=="decoh_rates"){ 

      MATRIX x( bp::extract<MATRIX>(params.values()[i]) );


      decoh_rates = new MATRIX(x.n_rows, x.n_cols);      
      for(int a=0;a<x.n_rows;a++){
        for(int b=0;b<x.n_cols;b++){
          decoh_rates->set(a, b, x.get(a,b));
        }
      } 

    }



  }

  sanity_check();

}


}// namespace libdyn
}// liblibra

