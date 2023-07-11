# QSH workflow

## step2  - generate the "short trajectory" data

  * ```cd step2```
  * ```mkdir res```
  * ```python gen_data.py```

### Description:
  * gen_data.py - script to generate the short time direct Hamiltonian files used for QSH simulation for tutorial purpose
                  or user can use the their own direct Hamiltonian files calculated from ab initio MD trajectory
  * the ```res``` directory will contain the model Hamiltonian files

## step3  - extend the trajectory data to longer times, using QSH

  * ```cd ../step3```
  * ```mkdir res_qsh```
  * ```python run_step3.py```

### Description:
  * The file ```res_qsh``` will contain the QSH files generated from the original input files (direct/model simulations)

## step4  - run the TSH calculations with original data and QSH data

  * ```cd ../step4```
  * Edit the file ```run_step4.py```
  * ```python run_step4.py```

### Description:
  * The files from the ```res_qsh``` or from ```res``` directories can be used and the input to this function
  In the former case, one obtains QSH dynamics, in the latter - the direct dynamics

  This script actually runs 3 variants of the TSH calculations:
  a) the standard one
  b) file-based QSH 
  c) file-less (on-the-fly) QSH





