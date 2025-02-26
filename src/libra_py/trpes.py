# *********************************************************************************
# * Copyright (C) 2025 Qingxin Zhang and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: trpes
   :platform: Unix, Windows
   :synopsis: this module implements calculations of time-resolved NAMD spectra (potential energy surface)

.. moduleauthor:: Qingxin Zhang, Alexey V. Akimov

"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from matplotlib.ticker import FormatStrFormatter

from liblibra_core import *
from libra_py import units
import libra_py.packages.cp2k.methods as cp2k
import util.libutil as comn

#%matplotlib inline

def compute_trpes(_params):
    """
    Args:
        params (dict): control parameters. Can contain keywords:

        - istep ( int ): the first file to read, should be the same as used in NAMD [default: 1]
        - fstep ( int ): the last file to read, should be the same as used in NAMD [default: 10]
        - dt ( int ): time step in fs [ default: 1.0 ]
        - namd_nsteps ( int ): how many NAMD steps were computed in the dynamics. This may be 
          larger than the number of files read - in this case the (fsteps - istep) files will be
          repeated until the desired number of steps, `namd_nsteps` is met [ default: fstep - istep - 1]
        - input_file_type ( int ): the kind of files to get the energy from:
           - 0 : matrix-like Hvib files [ default ]
           - 1 : CP2K logfiles
           - 2 : .npz matrix-like Hvib files 
        - logfile_read_params ( dict ) : additional parameters to read the logfiles (only if input_file_type == 1)
           [default : {} ]
        - units (string): the units of the energy [ default: "Ha"]. Possible: "Ha", "eV"
        - eprefix (string): the common prefix of the energy files to read [default: "Hvib_ci_"]
        - esuffix (string): the common suffix of the energy files to read [default: "_re"]        
        - dprefix (string): the common prefix of all the dynamics simulation results (the .hdf file) [default: "icond_"]
        - dsuffix (string): the common prefix of all the dynamics simulation results (the .hdf file) [default: "/mem_data.hdf"]
        - iconds (list of ints): the beginning timesteps for all the separate simulations [default: [1] ]
        - de ( float ): spacing of the energy grid points [default: 0.01, units: eV]
        - emin ( float ): minmal energy of the range [default: 0.0, units: eV]
        - emax ( float ): maximal energy of the range [default: 6.0, units: eV]
        - sigma_e ( float ): smearing parameter for the energy grid [ default: 0.1, units: eV] 

    Returns:
        tuple: 
        - np.array(nsteps): time grid [units: no units]
        - np.array(n_e_grid_pts): energy grid [units: eV]
        - np.array(nsteps, n_e_grid_pts): TRPES - time-resolved potential energy surface
        
    """

    params = dict(_params)
    comn.check_input(params, {"istep": 1, "fstep": 10, "dt":1.0, "namd_nsteps":9,
                              "eprefix":"Hvib_ci_", "esuffix":"_re",
                              "input_file_type":0, "logfile_read_params": {},
                              "units":"Ha",
                              "iconds":[1], "dprefix":"icond_", "dsuffix":"/mem_data.hdf",
                              "de":0.01, "emin":0.0, "emax":6.0, "sigma_e":0.1,
                             }, [])
    
    # ============ Read the parameters ======================
    istep = params["istep"]
    fstep = params["fstep"]
    nsteps = fstep - istep
    namd_nsteps = params["namd_nsteps"]
    dt = params["dt"]
    eprefix = params["eprefix"]
    esuffix = params["esuffix"]
    input_file_type = params["input_file_type"]
    logfile_read_params = dict(params["logfile_read_params"])
    eunits = params["units"]
    iconds = params["iconds"]
    dprefix = params["dprefix"]
    dsuffix = params["dsuffix"]
    de = params["de"]
    emin = params["emin"]
    emax = params["emax"]
    sigma_e = params["sigma_e"]  

    print(params)
   
    # ================== Read energies =====================
    E = []
    for step in range(istep, fstep):
        if input_file_type == 0:
            energy_filename = F"{eprefix}{step}{esuffix}"
            energy_mat = np.loadtxt(energy_filename)
            E.append(np.array(np.diag(energy_mat)))
        elif input_file_type == 1:
            logfile_read_params.update({"logfile_name":F"{eprefix}{step}{esuffix}"})
            e, _, _, _ = cp2k.read_cp2k_tddfpt_log_file(logfile_read_params)
            x = [0.0]
            for e_val in e:
                x.append(e_val)
            E.append(np.array(x))
        elif input_file_type == 2:
            energy_filename = F"{eprefix}{step}{esuffix}"
            energy_mat = sp.load_npz(energy_filename)
            E.append( np.array( np.diag( energy_mat.todense() ) ) )

    E = np.array(E)
    if eunits=="Ha":
        pass
    elif eunits=="eV":
        E = np.multiply(E, units.ev2Ha)
       
    #nsteps = E.shape[0]
    nstates = E.shape[1]

    # ================== Determine the time length ===============
    nstps = 0
    with h5py.File(F'{dprefix}{iconds[0]}{dsuffix}', 'r') as F:
        sh_pop = np.array(F['sh_pop_adi/data'])
        nstps = sh_pop.shape[0] # how many time-points

    print(F"Trajectory data contains {nstps}, while we only want {namd_nsteps}")
    nsteps = min(nstps, namd_nsteps)
    print(F"The TRPES plot will be of the time legnth of {namd_nsteps} points")

    max_icond = max(iconds)
    nfiles = fstep - istep
    mult = int( (max_icond - istep + nsteps) / nfiles) + 1   # how many times we need to repeat the initially read files 
                                                             # to ensure we have the energies for at least as long as 
                                                             # the namd trajectory    
    print(F"The multiplication factor is {mult}")
    E_ext = np.array(E)
    for _ in range(mult):
        E_ext = np.concatenate( (E_ext, E))
    E_ext = E_ext[:nsteps]
    print(F"The shape of the effective E array is {E_ext.shape}")
    
    # ================ Initialize all the grids =============

    n_e_grid_pts = int((emax - emin)/de)+1

    t_grid = np.linspace(0, nsteps, nsteps)
    e_grid = np.linspace(emin, emax, n_e_grid_pts)

    TRPES = np.zeros((nsteps, n_e_grid_pts))
    TRPES_traj = np.zeros((nsteps, n_e_grid_pts))


    # =============== Average over all iconds ===============
    print(F"Initial iconds = {iconds}")
    print(F"Renormalizing iconds for istep = {istep}, fstep = {fstep} and nsteps = {nsteps}...")
    iconds_eff = []
    for icond in iconds:
        if icond < fstep:
            iconds_eff.append(icond)
    print(F"The effective iconds are = {iconds_eff} (dropping off meaningless iconds)")
 
    count = 0
    for icond_indx, icond in enumerate(iconds_eff):
        try:
            sh_pop = np.array([])
            print(F"Reading the file {dprefix}{icond}{dsuffix}")
            with h5py.File(F'{dprefix}{icond}{dsuffix}', 'r') as F:
                sh_pop = np.array(F['sh_pop_adi/data'])

            E_eff = np.roll(E_ext, -(icond-istep), 0)
            for i in range(nsteps):                
                for j, e_j in enumerate(e_grid):
                    for state in range(nstates):
                        e = E_eff[i, state] * units.au2ev
                        smear = np.exp(-(e_j - e) ** 2) / (2 * sigma_e ** 2)
                        TRPES_traj[i, j] += sh_pop[i, state] * de * smear
                        
            TRPES += TRPES_traj
            count += 1

        except Exception as e:
            print(f"Error processing MSDM, icond={icond}: {e}")
    if count > 0:
        TRPES /= count

    TRPES /= np.max(TRPES) 

    return t_grid * dt, e_grid, TRPES
    

def plot_trpes(t_grid, e_grid, trpes, _plt_params):
    """
    Args:
        t_grid ( np.array(nsteps)) : time grid [units: no units]
        e_grid ( np.array(n_e_grid_pts) ): energy grid [units: eV]
        trpes ( np.array(nsteps, n_e_grid_pts)) : TRPES - time-resolved potential energy surface
        plt_params (dict): control parameters. Can contain keywords:
        
        - fig_size ( [2-ints list]): aspect ratio and size of the figure [ default: [8,3] ] 
        - cmap_name ( string ): name of the color map [ default: "rainbow" ]
        - cfont_size ( int ): size of the fonts of the color bar [ default: 20 ]
        - ctick_size ( int ): size of the labels for the color bar tics [ default: 16 ]
        - xlabel_size ( int ): font size for the x label [ default: 20 ]
        - xtick_size ( int ): font size for the x ticks [ default: 16]
        - ylabel_size ( int ): font size for the y label [ default: 20 ]
        - ytick_size ( int ): font size for the x ticks [ default: 16]
        - figurename ( string ) : name of the figure to produce [ default: "trpes.png" ]
        - dpi ( int ): figure resolution [ default: 300 ]
        - do_show ( Bool ): whether to show the figure or not [ default: False - no ]
        - do_save ( Bool ): whether to save the figure or not [ default: True - yes ]
        - time_units ( string ): the units of the time scale [ default: "fs" ]
        - energy_units ( string ): the units of the energy scale [ default: "eV" ]
        

    Returns:
        None: but produces/plots the figure
    """

    plt_params = dict(_plt_params)
    comn.check_input(plt_params, {"fig_size":[8,3], "cmap_name":"rainbow",
                                  "cfont_size": 20, "ctick_size": 16,
                                  "xlabel_size" : 20, "xtick_size": 16,
                                  "ylabel_size" : 20, "ytick_size": 16,
                                  "figurename" : "trpes.png", "dpi": 300,
                                  "do_show":False, "do_save":True,
                                  "time_units": "fs", "energy_units": "eV"
                                 }, [])

    # ============ Read the parameters ======================
    _fig_size = plt_params["fig_size"]
    cmap_name = plt_params["cmap_name"]
    _cfont_size = plt_params["cfont_size"]
    _ctick_size = plt_params["ctick_size"]
    _xlabel_size = plt_params["xlabel_size"]
    _xtick_size = plt_params["xtick_size"]
    _ylabel_size = plt_params["ylabel_size"]
    _ytick_size = plt_params["ytick_size"]
    figurename = plt_params["figurename"]
    _dpi = plt_params["dpi"]
    do_show = plt_params["do_show"]
    do_save = plt_params["do_save"]
    time_units = plt_params["time_units"]
    energy_units = plt_params["energy_units"]
    
    
    
    fig, ax = plt.subplots(1, 1, figsize=_fig_size, sharex=True)
    cmap = cm.get_cmap(cmap_name)
    img = ax.imshow(trpes.T, extent=[t_grid.min(), t_grid.max(), e_grid.min(), e_grid.max()], 
                    origin='lower', aspect='auto', cmap=cmap, vmin=0, vmax=1)

    ax.set_xlim(t_grid.min(), t_grid.max())
    ax.set_ylim(e_grid.min(), e_grid.max())
    ax.set_xticks(np.arange(t_grid.min(), t_grid.max() + 1, 500))
    ax.tick_params(axis='x', labelsize=_xtick_size)
    ax.tick_params(axis='y', labelsize=_ytick_size)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    cbar = fig.colorbar(img, ax=ax, aspect=20, shrink=1.0)  
    cbar.set_label("Intensity", fontsize=_cfont_size)
    cbar.ax.tick_params(labelsize=_ctick_size)

    ax.set_ylabel(F'Energy, {energy_units}', fontsize=_ylabel_size)
    ax.set_xlabel(F'Time, {time_units}', fontsize=_xlabel_size)

    if do_save:
        plt.savefig(F'{figurename}', dpi=300)
    if do_show:
        plt.show()
    

