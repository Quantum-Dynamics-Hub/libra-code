#*********************************************************************************
#* Copyright (C) 2025 Kosar Yasin, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************

"""
.. module:: plotting
   :platform: Unix, Windows
   :synopsis: This module implements functions for plotting excess energy, density of states, and UV-Vis spectra from CP2K TD-DFT calculations.

.. moduleauthor:: Kosar Yasin, Alexey V. Akimov
"""
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit
import h5py
import warnings
from scipy.constants import physical_constants
from libra_py import units
from libra_py.units import au2ev, au2fs
from libra_py.workflows.nbra import lz, step4
from libra_py import data_outs, data_conv



def load_adiabatic_energies(start_step, end_step, path_template, scale=au2ev, time_scale=units.au2fs):
    """Load adiabatic energies from files.

    Args:
        start_step (int): Start step for loading.
        end_step (int): End step for loading.
        path_template (str): Path to the file.
        scale (float): Scaling factor for energies.
        time_scale (float): Time scaling factor.

    Returns:
        tuple: Arrays of energies and times.
    """
    energies = []
    time = None
    for step in range(start_step, end_step):
        path = path_template.format(step=step)
        energy = sp.load_npz(path).todense().real
        if time is None:
            time = np.arange(energy.shape[0]) * time_scale
        energies.append(np.diag(energy))
    return np.array(energies) * scale, time

def calculate_excess_energies(adiabatic_energies, methods, icond_range, method_dir_prefix, method_dir_suffix):
    """Calculate and fit excess energies for different methods.

    Args:
        params (dict): Dictionary containing parameters for calculation.

    Returns:
        tuple: Dictionaries with results of excess energy calculations.
    """
    excess_energy_results = {}
    average_excess_energies = {}
    time_results = {}
    
    lower_bound_A, lower_bound_tau, lower_bound_beta = params["lower_bounds"]
    upper_bound_A, upper_bound_tau, upper_bound_beta = params["upper_bounds"]
    
    bounds = ([lower_bound_A, lower_bound_tau, lower_bound_beta], 
              [upper_bound_A, upper_bound_tau, upper_bound_beta])
    
    if 'BLLZ' in methods:
        bllz_results = calculate_bl_results(Hvib_params, BLLZ_params)
        n =  Hvib_params["nstates"]*3+4
        energy_BLLZ = bllz_results[:, n] * au2ev 
        time_BLLZ =  bllz_results[:, 0] * au2fs

        popt, _ = curve_fit(energy_fit, time_BLLZ, energy_BLLZ, p0=[energy_BLLZ[0], energy_BLLZ[-1], 500])
        energy_fit_BLLZ =  energy_fit(time_BLLZ, *popt)
        tau_BLLZ = popt[2]
    
    for method in methods:
        try:
            method_results = []
            taus = []
            method_times = None
        
            for icond in icond_range:
                    if method == 'BLLZ':
                        method_results.append(energy_BLLZ)
                        taus.append(tau_BLLZ)
                        if method_times is None:
                            method_times = time_BLLZ
                    
                    else:
                        path = f'{method}_{method_dir_prefix}_{icond}{method_dir_suffix}.hdf'
                        with h5py.File(path, 'r') as f:
                            sh_pop = np.array(f['sh_pop_adi/data'])
                            time = np.array(f['time/data'])[0:params["icond_end"]] * au2fs  
                        
                            if method_times is None:
                                method_times = time                    

                            tmp = np.sum(sh_pop[:len(time), :] * np.roll(adiabatic_energies[:len(time), :], -icond, axis=0), axis=1)
                            method_results.append(tmp)
                        
                            popt, _ = curve_fit(exponential_fit, time, tmp, bounds=bounds)
                            taus.append(popt[1])
                            
            excess_energy_results[method] = (method_results, taus)
            average_excess_energies[method] = np.mean(method_results, axis=0)
            time_results[method] = method_times
                
        except Exception as e:
            print(f"Error with {method} at icond {icond}: {e}")

    return excess_energy_results, average_excess_energies, time_results

def calculate_bl_results(Hvib_params, BLLZ_params):
    """Calculate results for BLLZ method.

    Args:
        Hvib_params (dict): Parameters for vibrational Hamiltonian.
        BLLZ_params (dict): Parameters for BLLZ method.

    Returns:
        np.ndarray: Array of results from BLLZ calculation.
    """
    Hvibs = step4.get_Hvib_scipy(Hvib_params)
    res, _ = lz.run(Hvibs, BLLZ_params)
    return data_conv.MATRIX2nparray(res, _dtype=float)

def density_of_states(Hvib_params, BLLZ_params):
    """Calculate density of states.

    Args:
        Hvib_params (dict): Parameters for vibrational Hamiltonian.
        BLLZ_params (dict): Parameters for BLLZ method.

    Returns:
        np.ndarray: Energy levels for density of states.
    """
    res0 = calculate_bl_results(Hvib_params, BLLZ_params)
    energy_levels = []
    for i in range(Hvib_params["nstates"]):
        energy_levels.extend(res0[:, 3 * i + 1] * units.au2ev)
    return np.array(energy_levels)




def process_spectra(log_file_pattern, output_folder, fwhm=0.1, num_points=1000):
    """Process multiple log files to generate an average UV-VIS spectrum.

    Args:
        log_file_pattern (str): Pattern to match log files.
        output_folder (str): Directory to save results or plots.
        fwhm (float): Full width at half maximum for Gaussian broadening.
        num_points (int): Number of points in the spectrum.

    Returns:
        tuple: Average x and y values for the spectrum or None if no data.

    Raises:
        FileNotFoundError: If no files match the given pattern.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    log_files = glob.glob(log_file_pattern)
    if not log_files:
        raise FileNotFoundError(f"No files matched the pattern: {log_file_pattern}")

    x_values_list, y_values_list = [], []
    
    for log_file in log_files:
        print(f"Processing: {log_file}")
        energy_levels, intensities = parse_spectrum_data_from_log(log_file)
        x, y = gaussian_broadening(energy_levels, intensities, fwhm, num_points)
        x_values_list.append(x)
        y_values_list.append(y)

    if x_values_list and y_values_list:
        x_values_avg = np.mean(x_values_list, axis=0)
        y_values_avg = np.mean(y_values_list, axis=0)
        return x_values_avg, y_values_avg
    else:
        return None, None

def plot(params, Hvib_params, BLLZ_params, log_file_pattern, output_folder):
    """Create plots for excess energy, density of states, and UV-VIS spectrum.

    Args:
        params (dict): Parameters for plotting. Expected keys include:
            - "start_step": Starting step for adiabatic energies (int)
            - "end_step": Ending step for adiabatic energies (int)
            - "path_template": Path template for energy files (str)
            - "methods": List of methods to calculate excess energies (list)
            - "icond_start", "icond_end", "icond_step": Indices for the calculation (int)
            - "method_dir_prefix", "method_dir_suffix": Parts of file name for methods (str)
            - "colors": Colors for plotting (list of str)
        Hvib_params (dict): Parameters for vibrational Hamiltonian. Expected keys include:
            - "nstates": Number of vibrational states (int)
            - Other parameters specific to the vibrational model
        BLLZ_params (dict): Parameters for BLLZ method. Expected keys include:
            - Any relevant parameters for the BLLZ method
        log_file_pattern (str): Pattern for log files to process.
        output_folder (str): Directory to save plots.

    Returns:
        None
    """
    plt.figure(figsize=(30, 22))
    gs = plt.GridSpec(2, 3, width_ratios=[4, 1, 2], height_ratios=[1, 4], wspace=0.05, hspace=0.05)
    ax0 = plt.subplot(gs[1, 0])
    ax1 = plt.subplot(gs[1, 1], sharey=ax0)
    
    adiabatic_energies, time_adiabatic = load_adiabatic_energies(params["start_step"], params["end_step"], params["path_template"])
    excess_energy_results, average_excess_energies, time_results =  calculate_excess_energies(adiabatic_energies, params["methods"], range(params["icond_start"], params["icond_end"], params["icond_step"]), params["method_dir_prefix"], params["method_dir_suffix"])

    for method, color in zip(params["methods"], params["colors"]):
        ax0.plot(time_results[method], average_excess_energies[method], label=method, color=color, linewidth=9) 
        time = time_results[method]
    
    for i in range(adiabatic_energies.shape[1]):
        ax0.plot(time, adiabatic_energies[:, i], linestyle='--', color='gray', alpha=0.5, linewidth=1.5)

    if ax0.get_legend_handles_labels()[1]:
        ax0.legend(fontsize=28, loc='upper center', bbox_to_anchor=(0.5, 0.23), ncol=3, frameon=False)
    
    energy_levels = density_of_states(Hvib_params, BLLZ_params)

    hist_values, _, _ = ax1.hist(energy_levels, bins=50, density=True, orientation='horizontal', alpha=0.7, color='gray', edgecolor='black', label='DOS')
    ax1.set_xlabel('DOS', fontsize=54)
    ax1.tick_params(axis='x', labelsize=38)
    ax1.yaxis.set_visible(False)


    x_uvvis, y_uvvis = process_spectra(log_file_pattern, output_folder)

    if x_uvvis is not None and y_uvvis is not None:
        max_hist_value = hist_values.max()
        y_normalized = (y_uvvis / y_uvvis.max()) * max_hist_value
        ax1.plot(y_normalized, x_uvvis, linestyle='-', color='blue', label='UV-VIS', linewidth=9)
        
    ax0.set_title('System', fontsize=60)
    ax0.set_xlabel('Time (fs)', fontsize=54)
    ax0.set_ylabel('Energy (eV)', fontsize=54)
    ax0.tick_params(axis='both', which='major', labelsize=38)

    ax1.set_xlabel('UV-VIS', fontsize=54)
    ax1.tick_params(axis='x', which='major', labelsize=38)

    handles2, labels2 = ax1.get_legend_handles_labels()
    ax1.legend(handles2, labels2, fontsize=26, loc='lower right', bbox_to_anchor=(1.0, 0.0), frameon=False)
    plt.savefig('spectra_dynamics.png', dpi=500)
    plt.show()
   
    return None



