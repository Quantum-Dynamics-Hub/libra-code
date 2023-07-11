1.  wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh .

2.  sh ./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -u -p /mnt/d/WORK/SOFTWARE/Conda

3.  edit .bashrc and restart Ubuntu

4.  conda update -n base -c defaults conda

5.  conda create -n libra

6.  conda activate libra

7.  Install packages

  conda install conda-build=3.22.0 

  conda install -c conda-forge gcc_linux-64=12.2.0 gxx_linux-64=12.2.0

  conda install make=4.2.1

  conda install -c conda-forge cmake=3.24.2

  conda install git
  
  conda install -c conda-forge boost-cpp=1.80.0

  conda install -c psi4 libint2=2.7.1

  conda install -c conda-forge eigen=3.4.0

  conda install -c conda-forge gmp=6.2.1

  conda install -c conda-forge mpfr=4.1.0

  conda install -c conda-forge numpy  (installs version 1.21.6, I think)

  conda install -c conda-forge scipy (installs version 1.7.3)

  conda install -c conda-forge matplotlib (installs version 3.5.2)

  conda install -c conda-forge jupyter_core (installs version 4.11.1)

  conda install -c conda-forge imageio=2.6.1

  conda install -c conda-forge llvm-openmp=14.0.4

  conda install -c conda-forge h5py=3.7.0

  conda install -c conda-forge py3dmol=1.8.1



  conda install -c conda-forge python-devtools  (for the headers)

  python3 -m pip install python-dev-tools --user --upgrade


  (one of the above operations downgrades Python to 3.7)


  conda install -c conda-forge boost=1.80.0


  conda install -c conda-forge ninja=1.11.0


