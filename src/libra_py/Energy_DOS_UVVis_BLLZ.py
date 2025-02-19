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
from libra_py.workflows.nbra import lz, step4
from libra_py import data_outs, data_conv

au2ev = physical_constants['hartree-electron volt relationship'][0]
au2fs = physical_constants['atomic unit of time'][0] * 1e15

def exponential_fit(t, A, tau, beta):
    """Fit an exponential decay function.

    Args:
        t (np.ndarray): Time vector.
        A (float): Amplitude of the decay.
        tau (float): Time constant for decay.
        beta (float): Exponential power.

    Returns:
        np.ndarray: Fitted values.
    """
    return A * np.exp(-np.power(t / tau, beta))

def energy_fit(t, E0, Einf, tau):
    """Fit an energy decay function.

    Args:
        t (np.ndarray): Time vector.
        E0 (float): Initial energy.
        Einf (float): Energy at infinite time.
        tau (float): Time constant for decay.

    Returns:
        np.ndarray: Fitted energy values.
    """
    return Einf + (E0 - Einf) * np.exp(-t / tau)

def load_adiabatic_energies(start_step, end_step, prefix, suffix, path, scale=au2ev, time_scale=units.au2fs):
    """Load adiabatic energies from files.

    Args:
        start_step (int): Start step for loading.
        end_step (int): End step for loading.
        prefix (str): Prefix for file naming.
        suffix (str): Suffix for file naming.
        path (str): Path to the file.
        scale (float): Scaling factor for energies.
        time_scale (float): Time scaling factor.

    Returns:
        tuple: Arrays of energies and times.
    """
    energies = []
    time = None
    for step in range(start_step, end_step):
        energy = sp.load_npz(path).todense().real
        if time is None:
            time = np.arange(energy.shape[0]) * time_scale
        energies.append(np.diag(energy))
    return np.array(energies) * scale, time

def calculate_excess_energies(params):
    """Calculate and fit excess energies for different methods.

    Args:
        params (dict): Dictionary containing parameters for calculation.

    Returns:
        tuple: Dictionaries with results of excess energy calculations.
    """
    excess_energy_results = {}
    average_excess_energies = {}
    time_results = {}
    
    bounds = ([lower_bound_A, lower_bound_tau, lower_bound_beta], 
              [upper_bound_A, upper_bound_tau, upper_bound_beta])

    if 'BLLZ' in params["methods"]:
        bllz_results = calculate_bl_results(params["Hvib_params"], params["BLLZ_params"])
        n = params["Hvib_params"]["nstates"] * 3 + 4
        energy_BLLZ = bllz_results[:, n] * au2ev 
        time_BLLZ = bllz_results[:, 0] * au2fs

        popt, _ = curve_fit(energy_fit, time_BLLZ, energy_BLLZ, p0=[energy_BLLZ[0], energy_BLLZ[-1], 500])
        energy_fit_BLLZ = energy_fit(time_BLLZ, *popt)
        tau_BLLZ = popt[2]

    for method in params["methods"]:
        try:
            method_results, taus, method_times = [], [], None
            for icond in range(params["icond_start"], params["icond_end"], params["icond_step"]):
                if method == 'BLLZ':
                    method_results.append(energy_BLLZ)
                    taus.append(tau_BLLZ)
                    if method_times is None:
                        method_times = time_BLLZ
                else:
                    path = f'{method}_{params["method_dir_prefix"]}_{icond}{params["method_dir_suffix"]}.hdf'
                    with h5py.File(path, 'r') as f:
                        sh_pop = np.array(f['sh_pop_adi/data'])
                        time = np.array(f['time/data'])[:params["icond_end"]] * au2fs  
                        
                        if method_times is None:
                            method_times = time                    

                        tmp = np.sum(sh_pop[:len(time), :] * np.roll(params["adiabatic_energies"][:len(time), :], -icond, axis=0), axis=1)
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

def parse_spectrum_data_from_log(file_path):
    """Parse energy levels and intensities from a CP2K log file.

    Args:
        file_path (str): Path to the log file.

    Returns:
        tuple: Arrays of energy levels and intensities.
    """
    energy_levels, intensities = [], []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(' TDDFPT|') and len(line.split()) >= 7:
                parts = line.split()
                try:
                    energy_levels.append(float(parts[2]))
                    intensities.append(float(parts[-1]))
                except ValueError:
                    continue
    return np.array(energy_levels), np.array(intensities)

def gaussian_broadening(energy_levels, intensities, fwhm, num_points=3500):
    """Apply Gaussian broadening to spectral data.

    Args:
        energy_levels (np.ndarray): Array of energy levels.
        intensities (np.ndarray): Array of corresponding intensities.
        fwhm (float): Full width at half maximum for Gaussian.
        num_points (int): Number of points for the broadened spectrum.

    Returns:
        tuple: Broadened energy and intensity arrays.

    Raises:
        ValueError: If input arrays are empty.
    """
    if len(energy_levels) == 0 or len(intensities) == 0:
        raise ValueError("Energy levels or intensities array is empty. Check input data.")
    x = np.linspace(np.min(energy_levels) - 5 * fwhm, np.max(energy_levels) + 5 * fwhm, num_points)
    y = np.sum([intensity * np.exp(-((x - energy) ** 2) / (2 * (fwhm / 2.35482) ** 2)) 
                for energy, intensity in zip(energy_levels, intensities)], axis=0)
    return x, y

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
        params (dict): Parameters for plotting.
        Hvib_params (dict): Parameters for vibrational Hamiltonian.
        BLLZ_params (dict): Parameters for BLLZ method.
        log_file_pattern (str): Pattern for log files to process.
        output_folder (str): Directory to save plots.

    Returns:
        None
    """
    plt.figure(figsize=(30, 22))
    gs = plt.GridSpec(2, 3, width_ratios=[4, 1, 2], height_ratios=[1, 4], wspace=0.05, hspace=0.05)
    ax0 = plt.subplot(gs[1, 0])
    ax1 = plt.subplot(gs[1, 1], sharey=ax0)
    
    adiabatic_energies, time_adiabatic = load_adiabatic_energies(params["start_step"], params["end_step"], params["adiabatic_dir_prefix"], params["adiabatic_dir_suffix"])
    excess_energy_results, average_excess_energies, time_results = calculate_excess_energies(params)

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

    return None




