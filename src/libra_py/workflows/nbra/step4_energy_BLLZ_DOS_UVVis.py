#*********************************************************************************
#* Copyright (C) 2024 Kosar Yasin, Alexey V. Akimov
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
   :synopsis: This module implements function for plotting the excess energy, density of state ,and UV-Vis spectrums of CP2K TD-DFT calculations
       Example:
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
from libra_py import units
from scipy.constants import physical_constants
from libra_py import units
from libra_py.workflows.nbra import lz, step4
from libra_py import data_outs, data_conv

au2ev = physical_constants['hartree-electron volt relationship'][0]
au2fs = physical_constants['atomic unit of time'][0] * 1e15  

# Exponential function for fitting
def exp_func_2(t, A, tau, beta):
    return A * np.exp(-np.power(t / tau, beta))

# Energy fitting function
def E_function(t, E0, Einf, tau):
    return Einf + (E0 - Einf) * np.exp(-t / tau)

# Loading adiabatic energies
def load_adiabatic_energies(start_step, end_step, adiabatic_dir_prefix, adiabatic_dir_suffix, scale=au2ev, time_scale=au2fs): 
    energies = []
    time = None
    for step in range(start_step, end_step):
        path = f'../{adiabatic_dir_prefix}_{step}_{adiabatic_dir_suffix}.npz'
        energy = sp.load_npz(path).todense().real
        if time is None:
            time_adiabatic = np.arange(energy.shape[0]) * time_scale
        energies.append(np.diag(energy))
    return np.array(energies) * scale, time_adiabatic

# Calculating excess energies
def excess_energies(adiabatic_energies, methods, icond_range, method_dir_prefix, method_dir_suffix):
    excess_energy_results = {}
    average_excess_energies = {}
    time_results = {}
    
    #bounds = ([0, 0, -np.inf], [np.inf, np.inf, np.inf])  # General case
    bounds = ([2.91 - 0.000001, 0.0, 0.5], [2.91 + 0.000001, np.inf, 4.0])
    
    if 'BLLZ' in methods:
        bllz_results = calculate_BL_results(Hvib_params, BLLZ_params)
        n =  Hvib_params["nstates"]*3+4
        energy_BLLZ = bllz_results[:, n] * au2ev 
        time_BLLZ =  bllz_results[:, 0] * au2fs

        popt, _ = curve_fit(E_function, time_BLLZ, energy_BLLZ, p0=[energy_BLLZ[0], energy_BLLZ[-1], 500])
        energy_fit_BLLZ = E_function(time_BLLZ, *popt)
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
                        
                            popt, _ = curve_fit(exp_func_2, time, tmp, bounds=bounds)
                            taus.append(popt[1])
                            
            excess_energy_results[method] = (method_results, taus)
            average_excess_energies[method] = np.mean(method_results, axis=0)
            time_results[method] = method_times
                
        except Exception as e:
            print(f"Error with {method} at icond {icond}: {e}")

    return excess_energy_results, average_excess_energies, time_results

# BLLZ
def calculate_BL_results(Hvib_params, BLLZ_params):
    Hvibs = step4.get_Hvib_scipy(Hvib_params)
    res, _ = lz.run(Hvibs, BLLZ_params)
    return data_conv.MATRIX2nparray(res, _dtype=float)

# DOS
def density_of_states(Hvib_params, BLLZ_params):
    res0 = calculate_BL_results(Hvib_params, BLLZ_params)
    energy_levels = []
    for i in range(Hvib_params["nstates"]):
        energy_levels.extend(res0[:, 3 * i + 1] * units.au2ev)
    energy_levels = np.array(energy_levels)
    return energy_levels

# Updated UV-VIS Spectrum Code
def parse_spectrum_data_from_log(file_path):
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
    if len(energy_levels) == 0 or len(intensities) == 0:
        raise ValueError("Energy levels or intensities array is empty. Check input data.")
    x = np.linspace(np.min(energy_levels) - 5 * fwhm, np.max(energy_levels) + 5 * fwhm, num_points)
    y = np.sum([intensity * np.exp(-((x - energy) ** 2) / (2 * (fwhm / 2.35482) ** 2)) 
                for energy, intensity in zip(energy_levels, intensities)], axis=0)
    return x, y

def process_spectra(log_file_pattern, output_folder, fwhm=0.1, num_points=1000):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    log_files = glob.glob(log_file_pattern)
    if not log_files:
        raise FileNotFoundError(f"No files matched the pattern: {log_file_pattern}")

    x_values_list = []
    y_values_list = []
    
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


# Plot
def plot(params, Hvib_params, BLLZ_params, log_file_pattern, output_folder):
    plt.figure(figsize=(30, 22))
    gs = plt.GridSpec(2, 3, width_ratios=[4, 1, 2], height_ratios=[1, 4], wspace=0.05, hspace=0.05)
    ax0 = plt.subplot(gs[1, 0])
    ax1 = plt.subplot(gs[1, 1], sharey=ax0)
    
    adiabatic_energies, time_adiabatic = load_adiabatic_energies(params["start_step"], params["end_step"], params["adiabatic_dir_prefix"], params["adiabatic_dir_suffix"])
    excess_energy_results, average_excess_energies, time_results = excess_energies(adiabatic_energies, params["methods"], range(params["icond_start"], params["icond_end"], params["icond_step"]), params["method_dir_prefix"], params["method_dir_suffix"])

    for method, color in zip(params["methods"], params["colors"]):
        ax0.plot(time_results[method], average_excess_energies[method], label=method, color=color, linewidth=9) 
        time = time_results[method]
    
    for i in range(adiabatic_energies.shape[1]):
        ax0.plot(time, adiabatic_energies[:, i], linestyle='--', color='gray', alpha=0.5, linewidth=1.5)

    if ax0.get_legend_handles_labels()[1]:
        ax0.legend(fontsize=28, loc='upper center', bbox_to_anchor=(0.5, 0.23), ncol=3, frameon=False)
    
    energy_levels = density_of_states(Hvib_params, BLLZ_params)

    hist_values, bin_edges, patches = ax1.hist(energy_levels, bins=50, density=True, orientation='horizontal', alpha=0.7, color='gray', edgecolor='black', label='DOS')
    ax1.set_xlabel('DOS', fontsize=54)
    ax1.tick_params(axis='x', labelsize=38)
    ax1.yaxis.set_visible(False)

    # UV-VIS Spectrum
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
    ax1.yaxis.set_visible(False)

    handles2, labels2 = ax1.get_legend_handles_labels()
    labels = []
    labels += labels2

    ax1.legend(handles2, labels, fontsize=26, loc='lower right', bbox_to_anchor=(1.0, 0.0), frameon=False)

    return None



# example of calling functions by user
params = {
    "methods": ['FSSH', 'IDA', 'MSDM', 'FSSH2', 'GFSH', 'MSDM_GFSH', 'BLLZ'],
    "colors": ['#3357ff', '#003f5c', '#FF5733', '#33A1FF', '#008080', '#FF33A6', 'red'],
    "start_step": 1000,
    "end_step": 4994,
    "icond_start": 1,
    "icond_end": 3994,
    "icond_step": 200,
    "adiabatic_dir_prefix": 'ustep3/res-mb-sd-DFT/Hvib_ci',
    "adiabatic_dir_suffix": 're',  # suffix of the filename, like "re", "im", or any suffix as needed
    "method_dir_prefix": 'latestnewNBRA_icond',
    "method_dir_suffix": '/mem_data'
}

Hvib_params = {
    "data_set_paths": ["/projects/academic/alexeyak/kosar/cp2k/fullerenes/c60-MOs60/ustep3/res-mb-sd-DFT/"],
    "Hvib_re_prefix": "Hvib_ci_", "Hvib_re_suffix": "_re", 
    "Hvib_im_prefix": "Hvib_ci_", "Hvib_im_suffix": "_im",
    "init_times": 1000, "nfiles": 3994, 
    "nstates": 75, "active_space": list(range(75))
}

BLLZ_params = { 
    "dt": 0.5*units.fs2au, "ntraj": 1, "nsteps": 3994, "istate": 54, 
    "Boltz_opt_BL": 1, "Boltz_opt": 1, "T": 300.0, "do_output": True, 
    "outfile": "BL.txt", "do_return": True, "evolve_Markov": True, 
    "evolve_TSH": False, "extend_md": False, "extend_md_time": 3994,
    "detect_SD_difference": False, "return_probabilities": False,
    "init_times": [0], "gap_min_exception": 0, "target_space": 0
}

log_file_pattern = '/projects/academic/alexeyak/kosar/cp2k/fullerenes/c60-MOs60/step2/all_logfiles/*.log'
output_folder = 'spectra_output'

plot(params, Hvib_params, BLLZ_params, log_file_pattern, output_folder)





