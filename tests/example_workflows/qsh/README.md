# Script to run quasi-stochastic Hamiltonian (QSH) NA-MD simulation

  * gen_data.py - script to generate the short time direct Hamiltonian files used for QSH simulation for tutorial purpose
                  or user can use the their own direct Hamiltonian files calculated from ab initio MD trajectory

  * run.py - python script to run the QSH simulation

  * Run the calculations by:   ```python run.py```

  * Need to make ```out``` and ```res``` directories before running the script. 
    the former contains the output files
    the latter contains the direct Hamiltonian files 


# output files descriptions

  * populations.txt - populations of each states

  * out/H_vib_re_E_acf_xx.txt - acf of each pair states

  * out/H_vib_re_E_spectrum_xx.txt - FT spectrum of each pair states

  * out/qsh_Ham_xx_re - QSH energies print at every time_intervals steps
 
  * out/qsh_Ham_xx_im - QSH couplings print at every time_intervals steps

