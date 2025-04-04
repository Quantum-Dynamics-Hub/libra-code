conda install conda-forge::gcc
conda install conda-forge::cmake
conda install conda-forge::make
conda install conda-forge::gxx


It is essential to use this version of LibTorch

Download here (cxx11 ABI):
https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.6.0%2Bcpu.zip


If need to reinstall pytorch:

    pip uninstall torch torchvision torchaudio
    pip cache purge


Installing:

    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu


How to install pytorch with cxx11 ABI:

conda install pytorch torchvision torchaudio -c pytorch  // but this is 2.3.1 PyTorch


Before the use, it is important to add the LibTorch libraries path:

export LD_LIBRARY_PATH=/home/alexvakimov/SOFTWARE/libpytorch/libtorch/lib:$LD_LIBRARY_PATH


