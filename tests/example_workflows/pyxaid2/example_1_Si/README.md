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

  * Edit the **step2_noSOC.py** file to define how many steps to prcess: e.g.
    ```nsteps_per_job = 20
       tot_nsteps = 100```
   
  * Run the calculations:  ```python step2_noSOC.py```



