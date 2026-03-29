# *********************************************************************************
# *
# * Copyright (C) 2020 Mohammad Shakiba, Brendan Smith, Alexey V. Akimov
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
.. module:: cube_file_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing with cube files
.. moduleauthors::
       Mohammad Shakiba, Brendan Smith, Alexey V. Akimov

"""


import os
import sys
import math
import re
import numpy as np

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *


import util.libutil as comn
import libra_py.packages.cp2k.methods as CP2K_methods

def read_cube(filename: str):
    """
    Read scalar field data (e.g., wavefunction or electron density) from a Gaussian
    .cube file and return it as a flattened 1D NumPy array.

    This function parses the cube file header to determine where the volumetric
    data begins, then reads all grid values and stores them in a single array.
    The ordering of values follows the convention used in the cube file
    (typically x fastest, then y, then z).

    Parameters
    ----------
    filename : str
        Path to the .cube file.

    Returns
    -------
    isovalues : numpy.ndarray
        1D array containing the scalar field values on the grid.

    Notes
    -----
    - The number of atoms is read from the third line of the file. Its absolute
      value is used to handle Gaussian cube files where this number may be negative.
    - The function skips:
        * 2 comment lines
        * 1 line with number of atoms and origin
        * 3 lines defining the grid
        * `natoms` lines with atomic coordinates
    - Some cube file variants include an extra line before the data block; this
      is handled automatically.
    - The output is flattened; reshaping into a 3D grid must be done separately
      using grid dimensions from the header if needed.
    """

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Number of atoms (may be negative in some Gaussian cube files)
    natoms = abs(int(lines[2].split()[0]))

    # Index of the first line containing volumetric data
    nstart = natoms + 2 + 1 + 3

    # Handle cube files with an extra line before data (format variation)
    if len(lines[nstart].split()) < 6:
        nstart += 1

    isovalues = []
    for line in lines[nstart:]:
        for val in line.split():
            isovalues.append(float(val))

    return np.array(isovalues)
    

def grid_volume(filename: str):
    """
    Compute the volume element (voxel volume) of a grid cell from a Gaussian
    .cube file.

    The cube file defines the volumetric grid using three lattice vectors
    (one per axis), given in lines 4–6 of the file. Each vector corresponds
    to the spacing and direction of the grid along x, y, and z. The volume
    of a single grid cell is the absolute value of the determinant of these
    three vectors.

    Parameters
    ----------
    filename : str
        Path to the .cube file.

    Returns
    -------
    dv : float
        Volume of a single grid cell (voxel) in Bohr³.

    Notes
    -----
    - Lines 4–6 of the cube file contain:
        * number of grid points along each axis (first column)
        * corresponding lattice vector components (remaining columns)
    - Only the vector components are used here.
    - The total grid volume would be `dv * Nx * Ny * Nz`, where Nx, Ny, Nz
      are the number of grid points along each axis.
    """

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Extract lattice vectors (skip the first column: number of grid points)
    axis_1 = [float(x) for x in lines[3].split()[1:4]]
    axis_2 = [float(x) for x in lines[4].split()[1:4]]
    axis_3 = [float(x) for x in lines[5].split()[1:4]]

    # Form the voxel matrix and compute its volume
    voxel = np.array([axis_1, axis_2, axis_3])
    dv = abs(np.linalg.det(voxel))

    return dv


def read_volumetric_data(filename: str):
    """
    Read volumetric data from a Gaussian .cube file and return it in a
    structured 3D form suitable for visualization and analysis.

    This function parses the cube file header to extract grid dimensions,
    lattice vectors, and atomic coordinates, then reshapes the volumetric
    data into a 3D array. It also constructs real-space grid coordinates
    corresponding to each voxel.

    Parameters
    ----------
    filename : str
        Path to the .cube file.

    Returns
    -------
    coordinates : numpy.ndarray
        Array of shape (natoms, 5) containing atomic data as read from the
        cube file (atomic number, charge, x, y, z). Values are returned as strings.

    x_grid, y_grid, z_grid : numpy.ndarray
        3D arrays of shape (nx, ny, nz) containing the Cartesian coordinates
        of each grid point.

    wave_fun : numpy.ndarray
        3D array of shape (nx, ny, nz) containing the volumetric data
        (e.g., wavefunction or electron density).

    spacing_vector : numpy.ndarray
        Vector representing the effective grid spacing (sum of the three
        lattice vectors). Often used in visualization routines.

    Notes
    -----
    - The cube file structure is:
        * 2 comment lines
        * 1 line: number of atoms and origin
        * 3 lines: grid size and lattice vectors
        * natoms lines: atomic coordinates
        * remaining lines: volumetric data
    - Grid vectors are given in Bohr; no unit conversion is applied.
    - The volumetric data is assumed to be ordered with x varying fastest,
      followed by y, then z (standard cube convention).
    """

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Number of atoms
    natoms = int(lines[2].split()[0])

    # Grid dimensions
    nx = int(lines[3].split()[0])
    ny = int(lines[4].split()[0])
    nz = int(lines[5].split()[0])

    # Lattice vectors
    axis_1 = np.array([float(x) for x in lines[3].split()[1:4]])
    axis_2 = np.array([float(x) for x in lines[4].split()[1:4]])
    axis_3 = np.array([float(x) for x in lines[5].split()[1:4]])

    # Effective spacing vector (useful for visualization libraries)
    spacing_vector = axis_1 + axis_2 + axis_3

    # Read volumetric data (flattened)
    data_start = natoms + 6
    isovals = []
    for line in lines[data_start:]:
        isovals.extend(float(val) for val in line.split())

    # Reshape into 3D array (cube convention)
    wave_fun = np.array(isovals).reshape((nx, ny, nz))

    # Construct coordinate grids
    x_grid = np.zeros((nx, ny, nz))
    y_grid = np.zeros((nx, ny, nz))
    z_grid = np.zeros((nx, ny, nz))

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                r = i * axis_1 + j * axis_2 + k * axis_3
                x_grid[i, j, k] = r[0]
                y_grid[i, j, k] = r[1]
                z_grid[i, j, k] = r[2]

    # Atomic coordinates (raw format from file)
    coordinates = np.array([lines[i].split() for i in range(6, natoms + 6)])

    return coordinates, x_grid, y_grid, z_grid, wave_fun, spacing_vector
    

def integrate_cube(cube_A, cube_B, grid_volume):
    """
    Compute the numerical integral of the product of two scalar fields
    defined on the same volumetric grid (e.g., wavefunctions from .cube files).

    The integral is approximated as a discrete sum over all grid points:
        ∫ A(r) B(r) dV ≈ Σ_i A_i * B_i * dv
    where dv is the volume of a single grid cell (voxel).

    Parameters
    ----------
    cube_A, cube_B : numpy.ndarray
        1D arrays containing the volumetric data (e.g., wavefunctions or
        densities) sampled on the same grid. These are typically obtained
        from a cube file reader (e.g., `read_cube`).
        Both arrays must have the same shape and ordering.

    grid_volume : float
        Volume of a single grid cell (voxel), typically computed using
        `grid_volume`. Units are usually Bohr³.

    Returns
    -------
    integral : float
        Numerical approximation of the integral ∫ A(r) B(r) dV.

    Notes
    -----
    - This operation corresponds to an overlap integral if A and B are
      wavefunctions defined on the same grid.
    - No normalization or unit conversion is performed.
    - The accuracy depends on the grid resolution and spacing.
    """

    # Element-wise product
    product = cube_A * cube_B

    # Discrete integration
    integral = product.sum() * grid_volume

    return integral



def plot_cubes(params):
    """
    This function plots the cubes for selected energy levels using VMD.

    Args:

        params (dict):

            min_band (int): The minimum state number.

            states_to_be_plotted (list): The list containing the Kohn-Sham orbitals to be plotted by VMD. This list is defined in the submit file.

            path_to_tcl_file (str): The path to the tcl file which contains the input for plotting the cubes in VMD.

            MO_images_directory (str): The molecular orbitals images directory.

            isUKS (int): This parameter is set for spin restricted and unrestricted calculations. When it is
                         set to 1 it means that unrestricted calculations were set in the input file otherwise
                         it is restricted.

            curr_step (int): The current time step used to save the images of the MOs.

            phase_factor_visual (numpy array): The phase correction factor list for each MOs for the current step.

    Returns:

        None

    """

   # Critical parameters
    critical_params = [
        "min_band",
        "states_to_be_plotted",
        "path_to_tcl_file",
        "MO_images_directory",
        "curr_step",
        "phase_factor_visual"]
    # Default parameters
    default_params = {"isUKS": 0}
    # Check input
    comn.check_input(params, default_params, critical_params)

    # Unpack the Kohn-Sham orbital indicies and the Kohn-Sham orbitals to be plotted. Also unpack the
    # path to the directory where the molecular orbitals will be plotted
    states_to_be_plotted = params["states_to_be_plotted"]
    # For VMD
    path_to_tcl_file = params["path_to_tcl_file"]
    # The molecular orbital images directory
    MO_images_directory = params["MO_images_directory"]

    # isUKS flag for spin-polarized and spin-unpolarized
    isUKS = int(params["isUKS"])
    # The current step
    curr_step = int(params["curr_step"])
    # If the path does not exist create it.
    if not os.path.isdir(MO_images_directory):
        os.makedirs(MO_images_directory)

    # Extracting the phase factor from params
    phase_factor_visual = params["phase_factor_visual"]

    min_band = params["min_band"]
    # We plot the cubes for the previous time step
    cubefile_names_prev = CP2K_methods.cube_file_names_cp2k(params)
    # read the lines of the tcl file
    tcl_file = open(path_to_tcl_file, 'r')
    tcl_lines = tcl_file.readlines()
    tcl_file.close()

    if isUKS == 1:
        # The cube file names of the alpha spin is the even indices of the cubefile_names_prev
        alp_cubefile_names_prev = cubefile_names_prev[0::2]
        # The same but for the phase factor of the alpha spin
        phase_factor_alpha = phase_factor_visual[0::2]
        # The cube file names of the beta spin is the odd indices of the cubefile_names_prev
        bet_cubefile_names_prev = cubefile_names_prev[1::2]
        # The same but for the phase factor of the beta spin
        phase_factor_beta = phase_factor_visual[1::2]

        for state_to_be_plotted in states_to_be_plotted:
            # Subtracting the min_band to obtain the index of the cube file for the state to be plotted for alpha spin
            alpha_cube_name = alp_cubefile_names_prev[state_to_be_plotted - min_band]
            # Subtracting the min_band to obtain the index of the cube file for the state to be plotted for beta spin
            beta_cube_name = bet_cubefile_names_prev[state_to_be_plotted - min_band]
            # open a new tcl file for alpha cubes
            new_file_alpha = open("vmd_alpha_cube_plot_%d.tcl" % curr_step, 'w')
            # open a new tcl file for beta cubes
            new_file_beta = open("vmd_beta_cube_plot_%d.tcl" % curr_step, 'w')

            for j in range(len(tcl_lines)):
                if 'mol load cube' in tcl_lines[j]:
                    # Open the cube file in VMD for alpha cubes
                    new_file_alpha.write('mol load cube %s\n' % alpha_cube_name)
                    # Open the cube file in VMD for beta cubes
                    new_file_beta.write('mol load cube %s\n' % beta_cube_name)
                elif 'render TachyonInternal' in tcl_lines[j]:
                    # Render the images to the MO_images_directory alpha cubes
                    new_file_alpha.write(
                        'render TachyonInternal %s/%s.tga\n' %
                        (MO_images_directory,
                         alpha_cube_name.replace(
                             'cubefiles/',
                             '').replace(
                             '.cube',
                             '')))
                    # Render the images to the MO_images_directory beta cubes
                    new_file_beta.write(
                        'render TachyonInternal %s/%s.tga\n' %
                        (MO_images_directory,
                         beta_cube_name.replace(
                             'cubefiles/',
                             '').replace(
                             '.cube',
                             '')))
                elif 'isosurface' in tcl_lines[j].lower():
                    # Correct the isovalues by multiplying it by the phase factor fo alpha orbitals
                    tmp_elements_alpha = tcl_lines[j].split()
                    tmp_elements_alpha[5] = str(phase_factor_alpha[state_to_be_plotted -
                                                min_band] * float(tmp_elements_alpha[5]))
                    isosurface_line_alpha = ' '.join(tmp_elements_alpha)
                    new_file_alpha.write(isosurface_line_alpha + '\n')

                    # Correct the isovalues by multiplying it by the phase factor for beta orbitals
                    tmp_elements_beta = tcl_lines[j].split()
                    tmp_elements_beta[5] = str(phase_factor_beta[state_to_be_plotted -
                                               min_band] * float(tmp_elements_beta[5]))
                    isosurface_line_beta = ' '.join(tmp_elements_beta)
                    new_file_beta.write(isosurface_line_beta + '\n')

                else:
                    # The rest of the tcl file lines
                    new_file_alpha.write(tcl_lines[j])
                    new_file_beta.write(tcl_lines[j])

            new_file_alpha.close()
            new_file_beta.close()
            # Run the VMD by tcl file
            os.system('vmd < vmd_alpha_cube_plot_%d.tcl' % curr_step)
            os.system('vmd < vmd_beta_cube_plot_%d.tcl' % curr_step)
            # os.system('rm vmd_alpha_cube_plot_%d.tcl' % curr_step)
            # os.system('rm vmd_beta_cube_plot_%d.tcl' % curr_step)

    else:
        # The same as above but with restricted spin calculations. No beta orbitals is considered.
        for state_to_be_plotted in states_to_be_plotted:
            cube_name = cubefile_names_prev[state_to_be_plotted - min_band]

            new_file = open("vmd_cube_plot_%d.tcl" % curr_step, 'w')
            for j in range(len(tcl_lines)):
                if 'mol load cube' in tcl_lines[j]:
                    new_file.write('mol load cube %s\n' % cube_name)
                elif 'render TachyonInternal' in tcl_lines[j]:
                    new_file.write(
                        'render TachyonInternal %s/%s.tga\n' %
                        (MO_images_directory,
                         cube_name.replace(
                             'cubefiles/',
                             '').replace(
                             '.cube',
                             '')))
                elif 'isosurface' in tcl_lines[j].lower():
                    tmp_elements = tcl_lines[j].split()
                    tmp_elements[5] = str(phase_factor_visual[state_to_be_plotted - min_band] * float(tmp_elements[5]))
                    isosurface_line = ' '.join(tmp_elements)
                    new_file.write(isosurface_line + '\n')
                else:
                    new_file.write(tcl_lines[j])

            new_file.close()

            os.system('vmd < vmd_cube_plot_%d.tcl' % curr_step)
            # os.system('rm vmd_cube_plot_%d.tcl' % curr_step)


def plot_cube_v2(params, cube_file_name, phase_factor):
    """
    This function plots the cube files using VMD

    Args:

        params (dictionary):

            vmd_input_template (string): The name of the VMD input template for visualizing the cube files.

            states_to_plot (list): The list of states that need to be plot.

            plot_phase_corrected (bool): Flag for plotting the molecular orbitals phase-corrected for the running job.

            vmd_exe (string): The VMD executable.

            tachyon_exe (string): The VMD Tachyon executable for rendering high-quality images.

            x_pixels (integer): Number of pixels in the X direction of the image.

            y_pixels (integer): Number of pixels in the Y direction of the image.

            remove_cube (bool): Flag for removing the cube files after plotting the molecular orbitals.

        cube_file_name (string): The name of the cube file

        phase_factor (integer): The phase factor for plotting the phase-corrected molecular orbital

    Returns:

        None
    """

    print(F'Plotting cube file {cube_file_name}...')
    file = open(params['vmd_input_template'], 'r')
    tcl_lines = file.readlines()
    file.close()
    vmd_exe = params['vmd_exe']
    tachyon_exe = params['tachyon_exe']
    x_pixels = params['x_pixels']
    y_pixels = params['y_pixels']
    image_format = params['image_format']
    together_mode = params['together_mode']
    if together_mode:
        new_tcl_name = 'vmd_tmode.tcl'
        file = open(new_tcl_name, 'w')
        cube_file_names = []
        for state in params['states_to_plot']:
            cube_file_names.append(cube_file_name.split('WFN')[0] + F'WFN_{str(state).zfill(5)}_1-1_0.cube')
        print(cube_file_names)
        state_counter = 0
        for i in range(len(tcl_lines)):
            if 'load cube' in tcl_lines[i]:
                tmp_name = cube_file_names[state_counter]
                file.write(F'mol load cube {tmp_name}\n')
                print(tmp_name)
                state_counter += 1
            elif 'render' in tcl_lines[i]:
                tmp_name = cube_file_name.split('WFN')[0]
                file.write(
                    F'render Tachyon {tmp_name} "{tachyon_exe} -aasamples 12 %s -format {image_format.upper()} -res {x_pixels} {y_pixels} -o %s.{image_format.lower()}"\n')
            else:
                file.write(tcl_lines[i])

        file.close()
    else:
        state_name = cube_file_name.replace('.cube', '')
        new_tcl_name = state_name + '.tcl'
        file = open(new_tcl_name, 'w')

        for i in range(len(tcl_lines)):
            if 'load cube' in tcl_lines[i]:
                file.write(F'mol load cube {cube_file_name}\n')
            elif 'Isosurface' in tcl_lines[i]:
                tmp = tcl_lines[i].split()
                for k in range(len(tmp)):
                    if tmp[k] == 'Isosurface':
                        break
                tmp[k + 1] = str(float(tmp[k + 1]) * phase_factor)
                tcl_lines[i] = ' '.join(tmp) + '\n'
                file.write(tcl_lines[i])
            elif 'render' in tcl_lines[i]:
                file.write(
                    F'render Tachyon {state_name} "{tachyon_exe} -aasamples 12 %s -format {image_format.upper()} -res {x_pixels} {y_pixels} -o %s.{image_format.lower()}"\n')
            else:
                file.write(tcl_lines[i])

        file.close()

    os.system(F'{vmd_exe} < {new_tcl_name}')
    if params['remove_cube']:
        if together_mode:
            for name in cube_file_names:
                os.system(F'rm {name}')
            os.system(F'rm {tmp_name}')
        else:
            os.system(F'rm {cube_file_name}')
            os.system(F'rm {state_name}')
