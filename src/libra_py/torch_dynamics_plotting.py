# *********************************************************************************
# * Copyright (C) 2025 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: torch_dynamics_plotting
   :platform: Unix, Windows
   :synopsis: This module implements the functions to read the HDF5 files with the results of the dynamical calculations
       List of functions:
           * plot_surfaces(_compute_model, _param_sets, states_of_interest, xmin, xmax, dx, plot_params)

.. moduleauthor:: Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"


import h5py
import numpy as np

import os
import sys
import math
import copy
import torch

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

import util.libutil as comn

import matplotlib.pyplot as plt   # plots
# from matplotlib.mlab import griddata



def plot_surfaces(_compute_model, _param_sets, states_of_interest, xmin, xmax, dx, _plot_params,
                  _ndof=1, _active_dof=0, _all_coordinates=[0.0]):
    """
    Args:
        _compute_model ( PyObject ): the function that returns the class with Hamiltonian properties
        _param_sets ( list of lists of dictionaries ): parameters of the model Hamiltonian, many sets are possible (hense the list)

            For each set, the following keywords are required:
            * **nstates** ( int ): the dimensionality of the Hamiltonian

        states_of_interest ( list of ints ): indices of the states we want to plot
        xmin ( double ): minimal value of the x axis used in the actual PES calculations  [a.u.]
        xmax ( double ): maximal value of the x axis used in the actual PES calculations  [a.u.]
        dx ( double ): step size of PES scan [a.u.]
        plot_params ( dictionary ): the parameters of plotting

            The dictionary can contain the following parameters:
            * **plot_params["colors"]** ( dictionary ) : defines the list of color definition (similar to the one found below)
            * **plot_params["clrs_index"]** ( list of strings ): the mapping of a color definition to its "name" (similar to the one found below)
            * **plot_params["xlim"]** ( list of 2 doubles ): the minimal and maximal values of the x axis in the plotted frame [ default: [-7.5, 15.0]]
            * **plot_params["ylim"]** ( list of 2 doubles ): the minimal and maximal values of the y axis in the plotted frame,
               for diabatic and adiabatic energies [default: [ -0.005, 0.025] ]
            * **plot_params["ylim2"]** ( list of 2 doubles ): the minimal and maximal values of the y axis in the plotted frame,
               for NACs [default: [-0.01, 0.01]]
            * **plot_params["do_show"]** (0 or 1) :
               - 0: - don't show the plot in the display (e.g. when there is no tunneling or proper terminal)
               - 1: - do show it. Be careful - sometimes showing is not recommended because you'll get a bunch of windows which
                 you'll need to close before the code would proceed to other parts. Turn it on if you use Jupyter notebooks. [ default: 1]
            * **plot_params["save_figures"]** (0 or 1): 0 - don't save the figures into separate files, 1 - do save them (e.g. when you can't or don't want
               to visualize them in terminal). The figures will be saved into a directory defined by the `prefix` variable. [default: 1]
            * **plot_params["prefix"]** (string): the name of the folder to which the figures will be saved (only if save_figures is set to 1). If the
               folder doesn't exist, it will be created [default: "out"]
            * **plot_params["dpi"]** (int): the quality of the figures to be saved [default: 300]
            * **plot_params["nac_idof"]** (int): when plotting derivative NACs, this variable defines the index of the nuclear DOF that defines such
               a NAC, e.g. <\\psi_i | d/dR_{nac_idof}| \\psi_j > [default: 0]
            * **plot_params["plotting_option"]** (0 or 1): how to plot the results:
               - 0 : Plot diabatic and adiabatic surfaces separately, plot projections too [ default ]
               - 1 : Plot diabatic and adiabatic surfaces in one picture, using dashed lines for diabatic. Plot NACs in a separate panel
            * **plot_params["figsize"]** (list of ints) : define the size of the figure [ default: [36, 18]]
            * **plot_params["titlesize"]** (int): size of the title [default: 38]
            * **plot_params["labelsize"]** (int): size of the label [default: 38]
            * **plot_params["fontsize"]** (int): size of the fonts [default: 36]
            * **plot_params["xticksize"]** (int): size of the x-ticks [default: 28]
            * **plot_params["yticksize"]** (int): size of the y-ticks [default: 28]
            * **plot_params["show_nac_abs"]** (int): when plotting derivative NACs, this will also add plotting of the modulus of NACs:
               - 0 : do not show [default]
               - 1 : do show it
            * **plot_params["margin_left"]** (double): left margin [ default: 0.2]
            * **plot_params["margin_right"]** (double): right margin [ default: 0.95]
            * **plot_params["margin_bottom"]** (double): bottom margin [ default: 0.13]
            * **plot_params["margin_top"]** (double): top margin [ default: 0.88]

        _ndof ( int ): the dimensionality of the PES [ default: 1 ]
        _active_dof ( int ): the index of the DOF used to construct the PES [ default: 0 ]
        _all_coodinates ( list of doubles ): values of all coordinates, the one at the position of `_active_dof` will be disregarded, while all
            other will be fixed at the values provided

    """

    # I know - this is a super-bad practice, but let's do it this way
    # import matplotlib.pyplot as plt

    colors_ = {}

    colors_.update({"11": "#8b1a0e"})  # red
    colors_.update({"12": "#FF4500"})  # orangered
    colors_.update({"13": "#B22222"})  # firebrick
    colors_.update({"14": "#DC143C"})  # crimson

    colors_.update({"21": "#5e9c36"})  # green
    colors_.update({"22": "#006400"})  # darkgreen
    colors_.update({"23": "#228B22"})  # forestgreen
    colors_.update({"24": "#808000"})  # olive

    colors_.update({"31": "#8A2BE2"})  # blueviolet
    colors_.update({"32": "#00008B"})  # darkblue

    colors_.update({"41": "#2F4F4F"})  # darkslategray

    clrs_index_ = ["11", "21", "31", "41", "12", "22", "32", "13", "23", "14", "24"]

    plot_params = dict(_plot_params)

    # Parameters and dimensions
    critical_params = []
    default_params = {"colors": colors_, "clrs_index": clrs_index_,
                      "xlim": [-7.5, 15], "ylim": [-0.005, 0.025], "ylim2": [-0.01, 0.01], "do_show": 1,
                      "save_figures": 1, "prefix": "out", "dpi": 300, "nac_idof": 0,
                      "plotting_option": 0, "figsize": [36, 18], "titlesize": 38, "labelsize": 38,
                      "fontsize": 36, "xticksize": 28, "yticksize": 28, "show_nac_abs": 0,
                      "margin_left": 0.2, "margin_right": 0.95, "margin_bottom": 0.13, "margin_top": 0.88,
                      "linewidth": 7
                      }
    comn.check_input(plot_params, default_params, critical_params)

    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]  # range of diabatic and adiabatic energies
    ylim2 = plot_params["ylim2"]  # range of NACs
    do_show = plot_params["do_show"]
    save_figures = plot_params["save_figures"]
    prefix = plot_params["prefix"]
    dpi_value = plot_params["dpi"]
    plotting_option = plot_params["plotting_option"]
    nac_idof = plot_params["nac_idof"]  # index of the DOF for which we plot the derivative coupling
    _figsize = plot_params["figsize"]
    _titlesize = plot_params["titlesize"]
    _labelsize = plot_params["labelsize"]
    _fontsize = plot_params["fontsize"]
    _xticksize = plot_params["xticksize"]
    _yticksize = plot_params["yticksize"]
    show_nac_abs = plot_params["show_nac_abs"]
    margin_left = plot_params["margin_left"]
    margin_right = plot_params["margin_right"]
    margin_top = plot_params["margin_top"]
    margin_bottom = plot_params["margin_bottom"]
    _linewidth = plot_params["linewidth"]

    fig_size = (_figsize[0], _figsize[1])

    X = []
    nsteps = int((xmax - xmin) / dx) + 1

    for i in range(nsteps):
        X.append(xmin + i * dx)

    plt.rc('axes', titlesize=_titlesize)      # fontsize of the axes title
    plt.rc('axes', labelsize=_labelsize)      # fontsize of the x and y labels
    plt.rc('legend', fontsize=_fontsize)      # legend fontsize
    plt.rc('xtick', labelsize=_xticksize)     # fontsize of the tick labels
    plt.rc('ytick', labelsize=_yticksize)     # fontsize of the tick labels

    plt.rc('figure.subplot', left=margin_left)
    plt.rc('figure.subplot', right=margin_right)
    plt.rc('figure.subplot', bottom=margin_bottom)
    plt.rc('figure.subplot', top=margin_top)

    sz = len(_param_sets)
    for iset in range(sz):

        # The "nstates" field must be specified
        model_params = dict(_param_sets[iset])
        comn.check_input(model_params, {}, ["nstates"])
        n = model_params["nstates"]
        nstates = n

        ham = sHamiltonian(1, nstates, _ndof)  # nbeads, nstates, nnucl
        
        """
        #ham.init_all(2)
        hdia, hadi, nac, nac_abs = [], [], [], []
        uij = []              # projecitions of the MOs onto elementary basis

        for k1 in range(nstates):
            hadi.append([])
            hdia.append([])
            nac_k1 = []
            nac_abs_k1 = []
            uij_k1 = []
            for k2 in range(nstates):
                uij_k1.append([])
                nac_k1.append([])
                nac_abs_k1.append([])
            uij.append(uij_k1)
            nac.append(nac_k1)
            nac_abs.append(nac_abs_k1)
        """
        
        hadi = torch.zeros(nsteps, nstates, nstates)
        hdia = torch.zeros(nsteps, nstates, nstates)
        nac = torch.zeros(nsteps, nstates, nstates)
        u = torch.zeros(nsteps, nstates, nstates)
        
        
        for i in range(nsteps):
            scan_coord = torch.zeros(1, _ndof, dtype=torch.float64) # MATRIX(_ndof, 1)
            
            for j in range(_ndof):
                scan_coord[0, j] = _all_coordinates[j]
            scan_coord[0, _active_dof] = X[i]

            # Diabatic properties
            ham.compute(_compute_model, scan_coord, model_params)

            # Adiabatic properties
            ham.dia2adi()
            ham.compute_nacs_and_grads()

           
            u[i,:,:] = cpp2py(ham.basis_transform).abs().pow(2) # as a tensor [1, nstates, nstates]
            # P = U * U.H()  # population matrix
            
            hadi[i,:,:] = cpp2py( ham.ham_adi )[0, :, :]
            hdia[i,:,:] = cpp2py( ham.ham_dia )[0, :, :]
            nac[i,:,:] = cpp2py( ham.dc1_adi )[0, _active_dof :, :]
            
            """
            for k1 in range(nstates):
                hadi[k1].append(ham.get_ham_adi().get(k1, k1).real)
                hdia[k1].append(ham.get_ham_dia().get(k1, k1).real)

                for k2 in range(nstates):
                    uij[k1][k2].append(U.get(k1, k2).real**2 + U.get(k1, k2).imag**2)
                    nac_k1_k2 = ham.get_dc1_adi(0).get(k1, k2).real
                    nac[k1][k2].append(nac_k1_k2)
                    nac_abs[k1][k2].append(abs(nac_k1_k2))
            """
        # print(hadi)
        nac_abs = torch.abs( nac )
        

        if plotting_option == 0:  # Plot diabatic and adiabatic surfaces separately, plot projections too

            findx = plt.gcf().number
            print("Current figre index is ", findx)

            # ==================== Diabatic and Adiabatic surfaces ==============
            # plt.figure(findx + 1, figsize=fig_size ) # dpi=300, frameon=False)
            plt.figure(figsize=fig_size)  # dpi=300, frameon=False)

            plt.subplot(1, 2, 1)
            plt.ylim(ylim[0], ylim[1])
            plt.xlim(xlim[0], xlim[1])

            # plt.title('Params set %i: Ham_dia' % (iset) )
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('Energy, a.u.')
            for k1 in states_of_interest:
                plt.plot(X, hdia[:, k1, k1], label='$H_{%i%i}$' %
                         (k1, k1), linewidth=_linewidth, color=colors[clrs_index[k1]])
            plt.legend()

            plt.subplot(1, 2, 2)
            plt.ylim(ylim[0], ylim[1])
            plt.xlim(xlim[0], xlim[1])
            # plt.title('Params set %i: Ham_adi' % (iset))
            plt.xlabel('Coordinate, a.u.')
            # plt.ylabel('Energy, a.u.')
            plt.yticks([])
            for k1 in states_of_interest:
                plt.plot(X, hadi[:,k1,k1], label='$E_{%i}$' % (k1), linewidth=_linewidth, color=colors[clrs_index[k1]])
            plt.legend()

            if save_figures:
                if not os.path.exists(F"{prefix}"):
                    os.system(F"mkdir {prefix}")
                plt.savefig(F"{prefix}/Ham_dia_E_adi_set_{iset}.png", dpi=dpi_value)

            if do_show:
                plt.show()

            # plt.clf()
            # plt.close()

            # =================== Projections ============================
            findx = plt.gcf().number
            print("Current figre index is ", findx)
            # plt.figure( findx + 1, figsize=fig_size ) # dpi=300, frameon=False)
            plt.figure(figsize=fig_size)  # dpi=300, frameon=False)
            sz1 = len(states_of_interest)

            for k2 in states_of_interest:
                indx = states_of_interest.index(k2)

                # plt.figure(2*iset+1+indx, figsize=(36, 18)) # dpi=300, frameon=False)
                # plt.figure(2*iset+1+indx, figsize=fig_size) # dpi=300, frameon=False)
                plt.subplot(sz1, 1, indx + 1)

                # plt.subplot(1, sz1, 1+indx)
                # plt.subplot(1, 1, 1+indx)
                # plt.title('Params set %i: Adi state %i' % (iset, k2) )
                plt.xlabel('Coordinate, a.u.')
                plt.ylabel('Projection')

                for k1 in range(nstates):
                    # plt.plot(X, uij[k1][k2], label='$dia_{%i} | adi_{%i}$' % (k1, k2), linewidth=7, color = colors[clrs_index[k1]])
                    plt.plot(X, u[:, k1,k2], label='$< \\psi^{dia}_{%i} | \\psi^{adi}_{%i} >$' % (
                        k1, k2), linewidth=_linewidth, color=colors[clrs_index[k1]])
            plt.legend()

            if save_figures:
                if not os.path.exists(F"{prefix}"):
                    os.system(F"mkdir {prefix}")
                plt.savefig(F"{prefix}/projections_set_{iset}.png", dpi=dpi_value)

            if do_show:
                plt.show()

            # plt.clf()
            # plt.close()

        elif plotting_option == 1:  # Plot diabatic and adiabatic surfaces in one picture, using dashed lines for diabatic
            # Plot NACs in a separate panel

            findx = plt.gcf().number
            plt.figure(figsize=fig_size)  # dpi=300, frameon=False)

            plt.subplot(1, 2, 1)
            plt.ylim(ylim[0], ylim[1])
            plt.xlim(xlim[0], xlim[1])

            # plt.title('Params set %i: Ham_dia' % (iset) )
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('Energy, a.u.')

            for k1 in states_of_interest:
                plt.plot(X, hdia[:, k1, k1], label='$H_{%i%i}$' %
                         (k1, k1), linewidth=_linewidth, ls="--", color=colors[clrs_index[k1]])
                plt.plot(X, hadi[:, k1, k1], label='$E_{%i}$' % (k1), linewidth=_linewidth, color=colors[clrs_index[k1]])
            plt.legend()

            plt.subplot(1, 2, 2)
            plt.ylim(ylim2[0], ylim2[1])
            plt.xlim(xlim[0], xlim[1])
            # plt.title('Params set %i: Ham_adi' % (iset))
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('NAC, a.u.')
            cnt = 0
            for k1 in states_of_interest:
                for k2 in states_of_interest:
                    if k2 > k1:
                        plt.plot(X, nac[:, k1, k2], label='$NAC_{%i%i}$' %
                                 (k1, k2), linewidth=_linewidth, color=colors[clrs_index[cnt]])
                        if show_nac_abs:
                            plt.plot(X, nac_abs[:,k1,k2], label='$NAC_{%i%i}$' % (
                                k1, k2), linewidth=_linewidth, ls="--", color=colors[clrs_index[cnt + 1]])
                        cnt = cnt + 1
            plt.legend()

            if save_figures:
                if not os.path.exists(F"{prefix}"):
                    os.system(F"mkdir {prefix}")
                plt.savefig(F"{prefix}/Ham_dia_E_adi_NAC_set_{iset}.png", dpi=dpi_value)

            if do_show:
                plt.show()
            # plt.clf()
            # plt.close()
    # plt.close()
