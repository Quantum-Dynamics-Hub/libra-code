# *********************************************************************************
# * Copyright (C) 2024 Kosar Yasin, Mohammad Shakiba,Qingxin Zhang, Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
#
#
#
"""

.. module:: plotting
   :platform: Unix, Windows
   :synopsis: This module implements function for plotting the results of NBRA calculations
       Example:
.. moduleauthor:: Kosar Yasin, Mohammad Shakiba,Qingxin Zhang, Alexey V. Akimov

"""

import os
import warnings

import matplotlib.pyplot as plt   # plots
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit
import h5py
from scipy.constants import physical_constants
from scipy.special import stdtrit

# from liblibra_core import *
# import libra_py.common_utils as comn
import util.libutil as comn
from libra_py import units
# from libra_py.workflows.nbra import lz, step4
# from libra_py import data_outs, data_conv

# Fitting functions


def exp_func(x, a):
    return 1 - np.exp(-x / a)


def gau_func(x, a):
    return 1 - np.exp(-x**2 / a**2)


def run(_params):
    """
    Plotting all the data for TD-DFT/NBRA caclulations

    Args:
        _params ( dictionary ): parameters controlling the function execution

            Required parameter keys:

            * **_params["nstates"]** ( int ): how many lines/columns in the file [Required!]

    Returns:
        None:
            plots the figures
    """

    params = dict(_params)

    critical_params = []
    default_params = {"title_1": "System", "filename": "excess_energy.png", "do_show": 0}
    comn.check_input(params, default_params, critical_params)

    title_1 = params["title_1"]  # the title of the figure
    filename = params["filename"]
    do_show = params["do_show"]

    fig, ax = plt.subplots(figsize=(16, 12))  # Larger figure size for the combined plot

    ax.set_title(F"{title_1}", fontsize=60)
    ax.set_xlabel('Time, fs', fontsize=54)
    ax.set_ylabel('Excess Energy, eV', fontsize=54)
    ax.tick_params(axis='both', which='major', labelsize=40)
    ax.legend(fontsize=24, loc='upper right', ncol=3, frameon=False)
    plt.tight_layout()
    plt.savefig(F"{filename}", dpi=600, bbox_inches='tight')
    if do_show:
        plt.show()


def fit_namd_data(
        methods,
        iconds,
        fitting_for_method,
        figuretitle,
        figurename="figure.png",
        conf_level=0.95,
        r2_min=0.8):
    # Create a 1x1 grid of subplots with a more square-like figure size
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), sharex=True)

    low_r2_iconds = []  # To store icond values where r_squared < r2_min

    # Loop through each method
    for c, method in enumerate(methods):
        taus = []
        fit_function = fitting_for_method[method]

        for icond in iconds:
            try:
                with h5py.File(f'{method}_icond_{icond}/mem_data.hdf', 'r') as F:
                    sh_pop = np.array(F['sh_pop_adi/data'][:, 0])
                    md_time = np.array(F['time/data'][:]) * units.au2fs
            except FileNotFoundError:
                print(f" File not found: {method}_icond_{icond}/mem_data.hdf, skipping.")
                continue

            popt, pcov = curve_fit(fit_function, md_time, sh_pop, bounds=([0.0], [np.inf]))
            tau = popt
            residuals = sh_pop - fit_function(md_time, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sh_pop - np.mean(sh_pop))**2)
            r_squared = 1.0 - (ss_res / ss_tot)

            print(f" icond: {icond}, r_squared: {r_squared}, tau: {tau}")

            # Check r_squared value and plot accordingly
            if r_squared >= r2_min:
                ax.plot(md_time, sh_pop, color='blue', linewidth=2.0)
                taus.append(tau)
                print(f" Plotting valid data with tau: {tau}")
            else:
                ax.plot(md_time, sh_pop, color='grey', linewidth=2.0)  # Gray line for low R^2
                low_r2_iconds.append(icond)  # Store icond if r_squared < r2_min

        taus = np.array(taus)
        if len(taus) > 0:
            ave_tau = np.average(taus)
            s = np.std(taus)
            N = taus.shape[0]
            Z = stdtrit(N - 1, conf_level)  # Compute confidence interval
            print(f" Z = {Z}")
            error_bar = Z * s / np.sqrt(N)

            ave_tau_ps = ave_tau / 1000  # Convert to picoseconds
            error_bar_ps = error_bar / 1000

            print(f" Timescales for {method}: {ave_tau_ps} ± {error_bar_ps} ps, averaged over {len(taus)} samples")

            # Plotting the red fit line with error bars
            ax.plot(md_time, fit_function(md_time, ave_tau - error_bar), ls='--', color="red", linewidth=2.0)
            ax.plot(
                md_time,
                fit_function(
                    md_time,
                    ave_tau),
                ls='-',
                color="red",
                label=f"{method}: {ave_tau_ps:.1f} ± {error_bar_ps:.1f} ps\n$R^2$: {r2_min} confidence: {conf_level*100}%",
                linewidth=4)
            ax.plot(md_time, fit_function(md_time, ave_tau + error_bar), ls='--', color="red", linewidth=2.0)
            ax.fill_between(
                md_time,
                fit_function(
                    md_time,
                    ave_tau -
                    error_bar),
                fit_function(
                    md_time,
                    ave_tau +
                    error_bar),
                color='red',
                alpha=0.1,
                linewidth=2.0)

    # Set plot parameters as in the original
    ax.tick_params(axis='x', labelsize=24, width=2.5, length=10)
    ax.tick_params(axis='y', labelsize=24, width=2.5, length=10)
    ax.set_xlim(left=-500, right=10500)
    ax.set_ylim(bottom=-0.001, top=0.021)
    ax.set_facecolor('white')

    for spine in ax.spines.values():
        spine.set_linewidth(2.5)

    ax.set_title(f'{figuretitle}', fontsize=26, pad=20, x=0.45)
    ax.set_xlabel('Time, fs', fontsize=26)
    ax.set_ylabel('Population', fontsize=26)

    ax.legend(fontsize=20, loc='upper left')
    fig.patch.set_facecolor('white')
    ax.set_xticks(np.arange(0, 40001, 10000))
    ax.set_ylim(0, 0.007)

    plt.tight_layout()
    plt.savefig(f'{figurename}', dpi=600)
    plt.show()

    # Print out the icond values with r_squared < r2_min
    if low_r2_iconds:
        print(f" iconds with r_squared < {r2_min}: {low_r2_iconds}")
    else:
        print(f" No iconds with r_squared < {r2_min}")
