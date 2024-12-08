#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
import os
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

methods = ['FSSH', 'IDA', 'MSDM', 'FSSH2', 'GFSH', 'MSDM_GFSH']
bright_colors = ['#3357ff', '#003f5c', '#FF5733','#33A1FF', '#008080', '#FF33A6']

plt.figure(figsize=(30, 22))
gs = plt.GridSpec(2, 3, width_ratios=[4, 1, 2], height_ratios=[1, 4], wspace=0.05, hspace=0.05)

ax0 = plt.subplot(gs[1, 0])
ax1 = plt.subplot(gs[1, 1], sharey=ax0)

adiabatic_energies = []
time = None

# Adiabatic energies
for step in range(start_step, end_step):
    energy = sp.load_npz(f'/path/to/Hvib_ci_{step}_re.npz').todense().real
    if time is None:
        time = np.arange(energy.shape[0]) * au2fs
    adiabatic_energies.append(np.diag(energy))

adiabatic_energies = np.array(adiabatic_energies)
adiabatic_energies = au2ev * adiabatic_energies

# Excess energies
for method, color in zip(methods, bright_colors):
    all_excess_energies = []
    taus = []


    for icond in range(icond_start, icond_end, icond_step):
        F = h5py.File(f'{method}_File_name_icond_{icond}/mem_data.hdf')
        sh_pop = np.array(F['sh_pop_adi/data'])
        time = np.array(F['time/data'])[:Add_icond_end] * au2fs  # Adjust time as necessary

        tmp2 = np.roll(adiabatic_energies[:Add_icond_end, :], -icond, axis=0)
        tmp1 = np.multiply(sh_pop[:Add_icond_end, :], tmp2)
        tmp = np.sum(tmp1, axis=1)
        all_excess_energies.append(tmp)

# Fit tau 
        try:
            popt, _ = curve_fit(exp_func_2, time, data, p0=initial_params, bounds=bounds)
            _, tau_fit, _ = popt
            taus.append(tau_fit)
        except Exception as e:
            print(f"Error fitting data for method {method}: {e}")

    avg_excess_energy = np.mean(all_excess_energies, axis=0)
    ax0.plot(time, avg_excess_energy, label=method, color=color, linewidth=9)


for i in range(adiabatic_energies.shape[1]):
    ax0.plot(time, adiabatic_energies[:, i], linestyle='--', color='gray', alpha=0.5, linewidth=1.5)

# BLLZ data
params = {
    "data_set_paths": ["/path/to/res-mb-sd-DFT/"],
    "Hvib_re_prefix": "Hvib_ci_", "Hvib_re_suffix": "_re", 
    "Hvib_im_prefix": "Hvib_ci_", "Hvib_im_suffix": "_im",
    "init_times": params["the number of first step"], "nfiles": params["the number of steps involved"], 
    "nstates": params["the number of states involved"], "active_space": list(range(Add_nstate))
}

Hvibs = step4.get_Hvib_scipy(params)
params = { 
    "dt": 0.5*units.fs2au, "ntraj": 1, "nsteps": params["the number of first step"], "istate": params["the number of first state"], 
    "Boltz_opt_BL": 1, "Boltz_opt": 1, "T": 300.0, "do_output": True, 
    "outfile": "BL.txt", "do_return": True, "evolve_Markov": True, 
    "evolve_TSH": False, "extend_md": False, "extend_md_time": params["the number of steps involved"],
    "detect_SD_difference": False, "return_probabilities": False,
    "init_times": [0]
}

params["gap_min_exception"] = 0 
params["target_space"] = 0  #

res, _ = lz.run(Hvibs, params)
res0 = data_conv.MATRIX2nparray(res, _dtype=float)

time_BLLZ = res0[:, 0] * units.au2fs  
energy_BLLZ = res0[:, (params['nstates']*3) + 4] * units.au2ev  


time_BLLZ = time_BLLZ[time_BLLZ <= 2000]  
energy_BLLZ = energy_BLLZ[:len(time_BLLZ)]  

# Fit energy data
popt, _ = curve_fit(E_function, time_BLLZ, energy_BLLZ, p0=[energy_BLLZ[0], energy_BLLZ[-1], 500])

time_fit = np.linspace(time_BLLZ.min(), time_BLLZ.max(), 500)
energy_fit = E_function(time_fit, *popt)

ax0.plot(time_BLLZ, energy_BLLZ, color='red', label="BLLZ", linewidth=6)

# Adiabatic Energy Levels
for i in range(adiabatic_energies.shape[1]):
    ax0.plot(time, adiabatic_energies[:, i], linestyle='--', color='gray', alpha=0.5, label=f'Adiabatic Level {i+1}' if i == 0 else "", linewidth=1.5)

# DOS 
energy_levels = []
for i in range(Add_number_of_state):
    energy_levels.extend(res0[:, 3 * i + 1] * units.au2ev)
energy_levels = np.array(energy_levels)

hist_values, bin_edges, patches = ax1.hist(energy_levels, bins=50, density=True, orientation='horizontal', alpha=0.7, color='gray', edgecolor='black', label='DOS')
ax1.set_xlabel('DOS', fontsize=54)
ax1.tick_params(axis='x', labelsize=38)
ax1.yaxis.set_visible(False)

# UV-VIS 
file_path = '/file/path/to/UV-VIS/'
x, y = [], []
with open(file_path, 'r') as file:
    for line in file:
        if line.strip():
            parts = line.split()
            x.append(float(parts[0]))
            y.append(float(parts[1]))

x = np.array(x) 
y = np.array(y)
max_y = max(y)

# Normalize UV-VIS 
max_hist_value = hist_values.max()
y_normalized = (y / y.max()) * max_hist_value

ax1.plot(y_normalized, x, linestyle='-', color='b', label='UV-VIS', linewidth=9)
ax1.set_xlabel('UV-VIS', fontsize=54)
ax1.tick_params(axis='x', which='major', labelsize=38) 
ax1.yaxis.set_visible(False) 


ax0.set_title('Title', fontsize=60)
ax0.set_xlabel('Time (fs)', fontsize=54)
ax0.set_ylabel('Energy (eV)', fontsize=54)
ax0.tick_params(axis='both', which='major', labelsize=38)


if ax0.get_legend_handles_labels()[1]:
    ax0.legend(fontsize=28, loc='upper center', bbox_to_anchor=(0.5, 0.23), ncol=3, frameon=False)

handles2, labels2 = ax1.get_legend_handles_labels()
labels = []
labels += labels2

ax1.legend(handles2, labels, fontsize=26, loc='lower right', bbox_to_anchor=(1.0, 0.0), frameon=False)

plt.savefig('Energy+DOS+UV-VIS.png', dpi=600)
plt.show()

