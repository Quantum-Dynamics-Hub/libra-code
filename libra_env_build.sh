conda create -n libra
conda activate libra
conda install -y conda-build
conda install -y gcc_linux-64
conda install -y gxx_linux-64
conda install -y cmake
conda install -y make
conda install -y git
conda install -y boost=1.73.0
conda install -y -c anaconda h5py
conda install -y -c conda-forge/label/gcc7 eigen
conda install -y -c psi4/label/dev libint2
conda install -y -c anaconda gmp
conda install -y -c conda-forge/label/gcc7 mpfr 
conda install -y python=3.6
conda install -y -c psi4 psi4
conda install -y -c conda-forge matplotlib
conda install -y -c rmg py3dmol
conda install -y -c anaconda numpy
conda install -y -c anaconda scipy 
conda install -y -c conda-forge/label/gcc7 imageio
conda install -y -c conda-forge llvm-openmp

