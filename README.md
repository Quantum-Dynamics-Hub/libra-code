# Libra

[![Build Status](https://travis-ci.org/Quantum-Dynamics-Hub/libra-code.svg?branch=master)](https://travis-ci.org/Quantum-Dynamics-Hub/libra-code)

This is the main page of the computational chemistry methodology discovery library, Libra

## About

* [extensive tutorials, example, and documentation](https://github.com/compchem-cybertraining/Tutorials_Libra).
This resource is actively developed, but the older tutorials may not fully reflect the current stage of the code itself, so the 
tutorials demonstrating the older features may be failing - please let us know if you need to use any of those and run into problems.

* [old program website](https://quantum-dynamics-hub.github.io/libra/index.html)
Not maintained for long time, but may still contain some useful info on some topics

## Community and Support

* Click below to join our community:
👉 [Join Slack](https://join.slack.com/t/quantumdynamicshub/shared_invite/zt-mjbhjssx-GGhsbYHxeBMvhmumK_j7LA)

* Open an Issue - to ask questions or report a problem

* The following [public forum](https://groups.google.com/forum/#!forum/quantum-dynamics-hub) exists, but
  hasn't been used for a while - use the Slack workspace instead


## Installation Videotutorials (as of 5/16/2022)

* [Installing WSL2](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=02184b70-7745-4eb4-a776-ae92014c652a&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)

* [Installing WSL2: After reboot](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=972aef79-e235-4a90-9ce1-ae92014d34db&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)

* [Installing Ubuntu on Windows 11](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=31a63536-f333-4242-9b56-ae92015ece64&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)

* [Creating Conda environment for Libra](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=d6ada23e-7e16-4b7a-b290-ae920188627c&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)

* [Installing Libra](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=7f8dd8c4-9f58-4ca0-a8cb-ae930166b7ec&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)

## Installation (as of after 3/23/2023)

### 1. Install miniconda (for Python 3.9) and activate Conda

#### 1.1 Download and install

    mkdir Conda
    cd Conda/
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh .
    sh ./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -u -p <install_dir>

  Here, 

  * the `-b` option will accept the license agreement and will skip all the prompts
  * the `-u` option will tell the installer to do all the needed updates
  * the `-p` option followed by the installation directory path (will be created), tells
     the installed where to install the package.

  Test it is working by doing:

    which conda


#### 1.2 Update the conda

  Actually, let's not do this:

    conda update -n base -c defaults conda


#### 1.3 Activate the environment

  Add the following line to you `.bashrc` or `.bash_profile` scripts:

    eval "$(<path to bin/conda> shell.bash hook)"

  For instance, 

    eval "$(/projects/academic/cyberwksp21/SOFTWARE/Conda/bin/conda shell.bash hook)"

  Restart your terminal or reload the `.bashrc` script:

    source ~/.bashrc


  When you do this, your command line should show up the (base) in front, indicating that
  the base environment is ready

  Test it is working by doing:

    which conda



### 2. Create the environment equipped with all Libra needs 

#### 2.1 Create the `libra` environment 

  In fact, you can call it whatever you like:

    conda create -n libra python=3.7

#### 2.2 Activate this environment

    conda activate libra

This is very important step - when activated, all the installs will go into that folder. 

**In case you mess up with an environment**, you can remove it with:

    conda remove --name libra --all


#### 2.3 Now, equip your environment with the required packages

Do this **one by one, and in this order**, (should not matter too much, but who knows...)

    > To automate the below procedures, you can use `-y` option to accept prompts (sometimes this will override)
    > previous packages/conflicts, so be careful
    > 
    > You can also use `-q` to get rid of all the messages to the output, although i'd keep it to keep track of what's going on


First let's install the most general packages:

    conda install -y -c conda-forge numpy scipy matplotlib imageio


Next, all what we actually need:

    conda install -y conda-build make
    conda install -y anaconda::py-boost
    conda install -y -c conda-forge gcc_linux-64=12.2.0 gxx_linux-64=12.2.0 cmake=3.24.2 python-devtools llvm-openmp
    conda install -y -c conda-forge/label/gcc7 eigen mpfr
    conda install -y -c psi4/label/dev libint2=2.7.1
    conda install -y -c anaconda h5py gmp

Install Jupyter Lab or traditional Jupyter notebook as explainted [here](https://jupyter.org/install):

    pip install -U jupyterlab

    or 

    pip install -U notebook


Install py3Dmol for viewing molecular structures:

    pip install -U py3Dmol


Installation instruction of Scikit-learn from its official website:

    pip install -U scikit-learn

 
    >
    >  YES - IT GOT SMALLER AND MORE COMPACT !
    >

Install PyTorch - since the current versions of Libra have a growing number of 
functions/modules written with PyTorch. The general instructions can be found [here](https://pytorch.org/get-started/locally/)

As a simple (most common case, peraps), we install PyTorch for CPU on Linux with:

    pip install -U torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

or

    conda install pytorch cpuonly -c pytorch
    

### 3. Download and build Libra

#### 3.1 Get it from the GitHub and choose the right branch

  Clone the repo from the GitHub

    git clone https://github.com/Quantum-Dynamics-Hub/libra-code.git libra

  and switch to the correct branch or tag - usually, it would be `devel` branch

    cd libra
    git checkout devel

#### 3.2 Create the build directory and make the Makefiles

   Then in the `libra` directory, create the build directory:

    mkdir _build
    cd _build
    cmake ../

#### 3.3 Compile the package
    
    make -j4


### 4. Make it ready to use

   Add the following exports to your `.bash_profile` file

    export PYTHONPATH=<path to the ppackage>/libra/_build/src:$PYTHONPATH


   Restart the terminal or source the bash profile and activate libra conda environment

    source .bash_profile 
    conda activate libra
    
   And you should be ready to use Libra.


## Useful notes:

### 1. 2/15/2025 (from Liz Stippell):

Notes on making Libra if you have python v3.7+ installed anywhere on your system (Linux):
Although the libra environment is made with python 3.7, during the `cmake ../` step it will search for any python, 
including versions outside of the libra environment. (Ex: my system kept finding python v3.9 in my 
Miniconda here: `/path/to/Conda/Miniconda3/include/python3.9` instead of searching within the libra 
environment: `/path/to/Conda/Miniconda3/envs/libra` )

To avoid this issue, you can add the following lines in the CMakeLists.txt file in your libra source code directory around line 44:

"""
set(Python3_ROOT_DIR "/path/to/Conda/Miniconda3/envs/libra")
set(Python3_EXECUTABLE "/path/to/Conda/Miniconda3/envs/libra/bin/python3")
set(Python3_LIBRARY "/path/to/Conda/Miniconda3/envs/libra/lib/libpython3.7m.so")
set(Python3_INCLUDE_DIR "/path/to/Conda/Miniconda3/envs/libra/include/python3.7m")
FIND_PACKAGE(Python3 3.6 REQUIRED COMPONENTS Development)
"""

This will force the make to search within the libra environment in the specified locations where you know the files exist. 


### 2. 4/17/2025 (Alexey Akimov)

A good way to setup the conda environment to have Boost and Python version consistent is this:


    conda install -c conda-forge boost=1.82 python=3.10


### 3. 4/17/2025 (Alexey Akimov)

Another useful recipe for setting up jupyter notebook specific to a selected Conda environment:

Step 1: Activate the environment

    conda activate libra


Step 2: Install ipykernel and register the kernel


    conda install ipykernel
    python -m ipykernel install --user --name=libra --display-name "Python (libra)"


Now, in Jupyter, you'll see a new kernel called "Python (libra)". Select that in your notebook.

### 4. 7/31/2025 (from Daeho Han)

Here are some installationinstructions from **Daeho Han** that may be used as a good revised summary of the above
installation instructions:

    conda create -n libra
    conda install -c conda-forge boost=1.82 python=3.10
    conda install -y -c conda-forge numpy scipy matplotlib imageio

    conda install -y conda-build make
    conda install -y -c conda-forge gcc_linux-64=12.2.0 gxx_linux-64=12.2.0 cmake=3.24.2 python-devtools llvm-openmp
    conda install -y -c conda-forge/label/gcc7 eigen mpfr
    conda install -y -c psi4/label/dev libint2=2.7.1
    conda install -y -c anaconda h5py gmp

    pip install -U jupyterlab
    pip install -U py3Dmol
    pip install -U scikit-learn
    pip install -U torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

## Developers and Contributors

  * Dr. Alexey Akimov (University at Buffalo, [link](https://akimovlab.github.io/index.html) )  
      The main developer and maintainer of the code
  
  * Dr. Daeho Han (University at Buffalo)
      Implementation of the exact-factorization (XF) methods, quantum trajectories surface hopping (QTSH),
      validating and testing Ehrenfest dynamics and other internals, implementing some model Hamiltonians 

  * Mr. Brendan Smith (University at Buffalo) 
      Entangled trajectories Hamiltonian, NA-MD with spin-orbit coupling, NBRA workflows, BL-LZ NA-MD tutorials and examples, 
      Libra/DFTB+, Libra/QE, Libra/ErgoSCF, Libra/CP2K, and Libra/Gaussian interfaces
      
  * Mr. Mohammad Shakiba (Shahid Bahonar University of Kerman, Iran)
      Cube file processing scripts, Libra/CP2K and Libra/Gaussian, Libra/Libint2 interfaces

  * Mrs. Story Temen (University at Buffalo)
      Implementation and testing of the HEOM codes

  * Dr. Wei Li (Hunan Agricultural University)
      NA-MD with spin-orbit coupling

  * Dr. Kosuke Sato (Toyota Research Lab) 
      State reordering scripts, Libra/GAMESS interface (Libra-X)

  * Dr. Ekadashi Pradhan (York University) 
      Libra/QE interface, delta-SCF NA-M (Libra-X)

  * Dr. Amber Jain (Indian Institute of Technology Bombay, India)
      Implementation and testing of the HEOM codes

  * Dr. Xiang Sun (NYU Shanghai, China)
      Implementation and testing of the FGR codes

  * Dr. Sophya Garashchuk (University of South Carolina)
      QTAG theory development

  * Dr. Matthew Dutra (University of South Carolina)
      Implementation and testing of the QTAG codes

## References

This code is provided in the hope it will be useful. 

  If you use Libra in your research, please cite the following paper:

  ### The "main" Libra papers 
  * [More recent overview of Libra's capabilities](https://doi.org/10.1016/j.simpa.2022.100445)
  Shakiba, M.; Smith, B.; Li, W.; Dutra, M.; Jain, A.; Sun, X.; Garashchuk, S.; Akimov, A.V.* "Libra: A modular software library for quantum
  nonadiabatic dynamics" *Software Impacts* **2022**  14, 100445
   
  * [The initial implementation](http://onlinelibrary.wiley.com/doi/10.1002/jcc.24367/full)
  Akimov, A. V. "Libra: An open-Source 'methodology discovery' library for quantum and classical dynamics simulations" 
  *J. Comput. Chem.*  **2016**  37, 1626-1649

  If you use any of the Libra's methods or implementation, please cite as appropriate:

  ### Papers that describe various developments in Libra
  * [FISH: Fully-integrated surface hopping](https://doi.org/10.1021/acs.jpclett.5c01509)
  Han, D.; Shakiba, M.; Akimov, A. V. "Fully-Integrated Surface Hopping as Quantum Decoherence Correction in Nonadiabatic Dynamics"
  *J. Phys. Chem. Lett.* **2025** 16, 7168-7176

  * [Multiple-state QTSH](https://doi.org/10.1021/acs.jctc.4c01751)
  Han, D.; Martens, C. C.; Akimov, A. V. "Generalization of Quantum-Trajectory Surface Hopping to Multiple Quantum States"
  *J. Chem. Theory Comput.* **2025** 20, 5022-5042

  * [F-tracking of excited states](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c02909)
   Akimov, A. V. "State Tracking in Nonadiabatic Molecular Dynamics Using Only Forces and Energies"
   *J. Phys. Chem. Lett.* **2024** 15, 11944-11953
    
  * [FSSH-3, modified FSSH-2, implmenentation of GFSH](https://www.tandfonline.com/doi/full/10.1080/00268976.2024.2376893)
  Akimov, A. V. "The Fewest Switches Surface Hopping as An Optimization Problem" *Mol. Phys.* **2024** e2376893

  * [Exact Factorization methods: SHXF, MQCXF, MFXF](https://doi.org/10.1021/acs.jctc.4c00343)
   Han, D.; Akimov, A.V. "Nonadiabatic Dynamics with Exact Factorization: Implementation and Assessment."
   *J. Chem. Theory Comput.* **2024** 20, 5022–5042

  * [ML Hamiltonian mapping approach](https://doi.org/10.1021/acs.jctc.4c00008)
  Shakiba, M.; Akimov, A.V. "Machine-Learned Kohn–Sham Hamiltonian Mapping for Nonadiabatic Molecular Dynamics."
  *J. Chem. Theory Comput.* **2024** 20, 2992–3007

  * [TC-NBRA](https://doi.org/10.1021/acs.jpclett.3c03029)
  Akimov, A.V. "Energy-Conserving and Thermally Corrected Neglect of Back-Reaction Approximation Method for Nonadiabatic Molecular Dynamics"
  *J. Phys. Chem. Lett.* **2023** 14, 11673-11683

  * [QTAG](https://doi.org/10.1002/qua.27078)
  Dutra, M.; Garshchuk, S.; Akimov, A. "The Quantum Trajectory-guided Adaptive Gaussian Methodology in the Libra Software Package"
  *Int. J. Quntum Chem.* **2023.** 123, e27078

  * [(generalized) Local diabatization](https://doi.org/10.1007/s00214-023-03007-7)
  Shakiba, M.; Akimov, A.V. "Generalization of the Local Diabatization Approach for Propagating Electronic Degrees of Freedom in Nonadiabatic Dynamics"
  *Theor. Chem. Acc.* **2023** 142, 68

  * [xTB from cp2k/Libra interface for NA-MD in large-scale systems, k-point formulation, spin-adaptation](https://doi.org/10.1021/acs.jctc.2c00297)
  Shakiba, M.; Stippel, E.; Li, W.; Akimov, A. V. "Nonadiabatic Molecular Dynamics with Extended Density Functional Tight-Binding:
  Application to Nanocrystals and Periodic Solids" *J. Chem. Theory Comput.* **2022** 18, 5157-5180

 * [many-body effects, NA-MD with TD-DFT states](https://pubs.acs.org/doi/10.1021/acs.jctc.0c01009)
  Smith, B.; Shakiba, M.; Akimov, A. V. "Nonadiabatic Dynamics in Si and CdSe Nanoclusters: Many-Body vs. Single-Particle
  Treatment of Excited States" *J. Chem. Theory. Comput.* **2021** 17, 678-693

  * [HEOM implementation](https://pubs.acs.org/doi/10.1002/qua.26373)
  Temen, S.; Jain, A.; Akimov, A. V. "Hierarchical equations of motion in the Libra software package"
  *Int. J. Quant. Chem.* **2020** 120, e26373

  * [Belyaev-Lebedev-Landau-Zener Surface Hopping within the Neglect of Back-Reaction Approximation](https://pubs.acs.org/doi/10.1021/acs.jpclett.9b03687)
  Smith, B.; Akimov, A. V. "Hot Electron Cooling in Silicon Nanoclusters via Landau-Zener Non-Adiabatic Molecular Dynamics: 
  Size Dependence and Role of Surface Termination" *J. Phys. Chem. Lett.* **2020** 11, 1456-1465 

  * [Phase correction, Ehrenfest dynamics details, basis transformations (see the SI)](https://doi.org/10.1021/acs.jpclett.8b02826)
  Akimov, A. V.; "A Simple Phase Correction Makes a Big Difference in Nonadiabatic Molecular Dynamics"
  *J. Phys. Chem. Lett.*  **2018** 9, 6096-6102



  You may find the following papers useful examples
  ### Papers that utilize Libra

  * [Formulation of a fragment-based NA-MD](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00955)
  Akimov, A. V. "Nonadiabatic Molecular Dynamics with Tight-Binding Fragment Molecular Orbitals" 
  *J. Chem. Theory Comput.*  **2016** 12, 5719-5736

  * [Quasi-stochastic Hamiltonian for longer NA-MD](http://pubs.acs.org/doi/abs/10.1021/acs.jpclett.7b02185)
  Akimov, A. V.; "Stochastic and Quasi-Stochastic Hamiltonians for Long-Time Nonadiabatic Molecular Dynamics"
  *J. Phys. Chem. Lett.*  **2017** 8, 5190-5195

  * [Entrangled-trajectories Hamiltonian dynamics to capture quantum effects of nuclei](https://doi.org/10.1063/1.5022573)
  Smith, B. A.; Akimov, A. V. "Entangled trajectories Hamiltonian dynamics for treating quantum nuclear effects" 
  *J. Chem. Phys.* **2018** 148, 144106

  * [Inclusion of the Spin-orbit coupling in NA-MD](https://doi.org/10.1021/acsenergylett.8b01226)
  Li, W.; Zhou, L.; Prezhdo, O. V.; Akimov, A. V. "Spin-Orbit Interactions Greatly Accelerate Nonradiative Dynamics
  in Lead Halide Perovskites" *ACS Energy Lett.* **2018** 3, 2159-2166




