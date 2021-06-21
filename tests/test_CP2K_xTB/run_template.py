import os
import sys
from libra_py.workflows.nbra import step2


params = {}
params['nprocs'] = 36
params['istep'] = 
params['fstep'] = 
params['init_state'] = 2200
params['final_state'] = 2600
params['is_spherical'] =  True
params['remove_molden'] = True
params['res_dir'] = '/home/97425008/Libra/tests/test_CP2K_xTB/res'
params['all_pdosfiles'] = '/home/97425008/Libra/tests/test_CP2K_xTB/all_pdosfiles'
params['all_logfiles'] = '/home/97425008/Libra/tests/test_CP2K_xTB/all_logfiles'
params['cp2k_exe'] = '/home/97425008/cp2k-v7/cp2k/exe/local/cp2k.popt'
params['cp2k_ot_input_template'] = '/home/97425008/Libra/tests/test_CP2K_xTB/es_ot_temp.inp'
params['cp2k_diag_input_template'] = '/home/97425008/Libra/tests/test_CP2K_xTB/es_diag_temp.inp'
params['trajectory_xyz_filename'] = '/home/97425008/Libra/tests/test_CP2K_xTB/Si_QD_xtb_3.5nm_0.5fs-pos-1.xyz'

step2.run_cp2k_xtb_step2(params)

