# Installation

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

Here are some installation instructions from **Daeho Han** that may be used as a good summary of the above content:

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
