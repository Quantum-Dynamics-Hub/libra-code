# Computing the molecular orbital overlaps for xTB calculations

Here are the files for computing the MO overlaps for extended tight-binding calculations using CP2K. 
We first need to run an OT calculations and then use the `wfn` file for diagonalization so that we can print out the 
MOs in `molden` or `MOLog` file formats. 

Then we will use `libint2` to compute the AO overlap matrix to compute the MO overlaps.

We first need to specify the correct variables in the `es_ot_template.inp` and `es_diag_template.inp` and 
then specify the full paths to these files:

* `res_dir`
* `all_pdosfiles`
* `all_logfiles`
* `cp2k_exe`
* `cp2k_ot_input_template`
* `cp2k_diag_input_template`
* `trajectory_xyz_filename`

It is required to specify the 

* `nprocs` 
* `init_state`
* `final_state` in the `run_template.py` 

but the `istep`, `fstep`, and `njobs` is defined in `distribute_jobs.py`. 
If you want to submit the jobs (and not running as bash), set `run_slurm` to `True`.

After specifying all the variables in both files you can run the calculations using `python distribute_jobs.py`

The MO overlaps are saved as `.npz` format. 
These files are binary files by `scipy.sparse` library. In order to load the saved files you can do this:

```
import numpy as np
import scipy.sparse

# Load the sparse matrix
S = scipy.sparse.load_npz('./res/S_ks_0_re.npz')
# The dense matrix
S_dense = S.todense()

# Show the diagonal elements of the dense matrix
np.diag(S_dense)
```


