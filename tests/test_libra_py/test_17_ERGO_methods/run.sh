./ergo <<EOINPUT > /dev/null
spin_polarization = 0
charge = 0
molecule_inline
F    0  0  0.0
H    0  0  1.0
EOF
basis = "6-31G"
use_simple_starting_guess=1
tmpdir = "./tmp"
set_nthreads("detect")
scf.create_mtx_files_D = 1
scf.create_mtx_files_F = 1
scf.create_mtx_file_S = 1
scf.output_homo_and_lumo_eigenvectors = 1
scf.number_of_occupied_eigenvectors = 3
scf.number_of_unoccupied_eigenvectors = 2
scf.eigenvectors_method = "projection"
XC.sparse_mode = 1
scf.force_unrestricted = 1
scf.starting_guess_disturbance = 0.01

run "HF"
