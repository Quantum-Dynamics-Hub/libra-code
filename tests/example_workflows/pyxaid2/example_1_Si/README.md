# Step 1

  * Go to *step1* directory:  ```cd step1```

  * Edit the submit.slm file to change your email (or comment it out like this:)
   ```###SBATCH --mail-user=alexeyak@buffalo.edu```

  * Run the calculations by:   ```sbatch submit.slm```


# Step 2

  * Go to *step2* directoy:  ```cd ../step2```

  * Copy the file **x0.md.out** produced in the previous step to here: ```cp ../step1/x0.md.out .```

  * Edit the **submit_templ.slm** file to define the absolute path of the directory that will collect
    the results of this step:  e.g. ```res=/gpfs/scratch/alexeyak/example_1_Si/step2/res``` 
    Here, the *res* directory is (or will be created in) the current *step2* folder, for instance. 

  * Edit the **run_step2.py** file to define how many steps to prcess: e.g.
    ```nsteps_per_job = 20
       tot_nsteps = 100```
   
  * Run the calculations:  ```python run_step2.py```


# Step 3

  * Go to the *step3* directory:  ```cd ../step3```

  * In the **run_step3.py** edit the variable ```setup``` - select from options 1 or 2

  * Depending on the setup, create a folder *res_setup1* or *res_setup2*, e.g. ```mkdir res_setup1``` 
    the name can be different, but should be consistent with the definition in ```params["output_set_paths"]```

  * Edit other variables, if needed

  * Run the calculations:  ```python run_step3.py```


# Step 4

  * Go to the *step4* directory:  ```cd ../step4```

  * In the **run_step4.py** edit the desired variables. In particular, note the
    line containing ```params["data_set_paths"].append(...)```. If uncommented, it mimics the presence 
    of another directory (potentially for a data corresponding to another MD trajectory)

    Note the difference btween the ```params["nfiles"]``` and ```params["nsteps"]``` variables

  * Run the calculations:  ```python run_step4.py```
