# Installing Libra with Conda

## 1. Prerequisites

Install either:

* Miniconda
* Anaconda
* Miniforge (recommended)


```bash
mkdir Conda
cd Conda/
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh .
sh ./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -u -p <install_dir>
```
Here,

  * the `-b` option will accept the license agreement and will skip all the prompts
  * the `-u` option will tell the installer to do all the needed updates
  * the `-p` option followed by the installation directory path (will be created), tells
     the installed where to install the package.

Verify that Conda is available:

```bash
conda --version
```

Equip your miniconda with the basic tools, needed for the next steps, such as `git`:

```bash
conda install -y -c conda-forge git
```

## 2. Clone the repository and choose the branch to build

```bash
git clone https://github.com/Quantum-Dynamics-Hub/libra-code.git
cd libra-code
git checkout devel
```
Most of the time, the `main` version is behind the current development (`devel`)) version
so we often want to switch to the correct branch.


## 3. Create and activate the Libra environment

Create it using the `environment.yml` file located in the rood directory of Libra

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate libra
```

## 4. Configure and build

```bash
mkdir _build
cd _build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} ..
make -j2
```

## 5. Setup the environment variables

Add the following line to you `.bashrc` or `.bash_profile` scripts:

```bash
eval "$(<path to bin/conda> shell.bash hook)"
```

For instance,
```bash
eval "$(/projects/academic/cyberwksp21/SOFTWARE/Conda/bin/conda shell.bash hook)"
```

Restart your terminal or reload the `.bashrc` script:
```bash
source ~/.bashrc
```

When you do this, your command line should show up the (base) in front, indicating that
the base environment is ready

## 6. Verify the installation

```bash
python -c "import liblibra_core"
```

If no errors are reported, Libra was built successfully.


# Recommended additional packages for your Conda environment

Make sure to have your desired Conda environment activated before you do the following:

## 1. Jupyter notebook or JupyterLab

Install Jupyter Lab or traditional Jupyter notebook as explainted [here](https://jupyter.org/install):
```bash 
pip install -U jupyterlab
```

or
```bash
pip install -U notebook
```

## 2. py3Dmol

Install py3Dmol for viewing molecular structures:

```bash
pip install -U py3Dmol
```

# How to add the environment to Jupyter notebook

Another useful recipe for setting up jupyter notebook specific to a selected Conda environment:

## 1. Activate the environment

```bash
conda activate libra
```

## 2. Install ipykernel and register the kernel
```bash
conda install ipykernel
python -m ipykernel install --user --name=libra --display-name "Python (libra)"
```

Now, in Jupyter, you'll see a new kernel called "Python (libra)". Select that in your notebook before you do the calculations that require Libra


# Installation of WSL

## Installation Videotutorials (as of 5/16/2022)

* [Installing WSL2](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=02184b70-7745-4eb4-a776-ae92014c652a&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)

* [Installing WSL2: After reboot](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=972aef79-e235-4a90-9ce1-ae92014d34db&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)

* [Installing Ubuntu on Windows 11](https://ub.hosted.panopto.com/Panopto/Pages/Embed.aspx?id=31a63536-f333-4242-9b56-ae92015ece64&autoplay=false&offerviewer=true&showtitle=true&showbrand=true&captions=false&interactivity=all)


