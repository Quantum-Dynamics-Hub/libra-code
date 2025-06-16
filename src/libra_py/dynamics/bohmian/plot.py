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
.. module:: plot
   :platform: Unix
   :synopsis: This module implements functions for plotting results of Bohmian calculations 
       List of functions:

.. moduleauthor:: Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2025 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://github.com/Quantum-Dynamics-Hub/libra-code"

import torch
import matplotlib.pyplot as plt
import libra_py.dynamics.tsh.plot as tsh_plot
import libra_py.units as units

def plot_scatter(_plt, t, q_traj, plot_params):

    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    which_energies = plot_params["which_energies"]

    if xlim is not None:
        _plt.xlim(xlim[0], xlim[1])
    if ylim is not None:
        _plt.ylim(ylim[0], ylim[1])

    _plt.title("Dynamics", fontsize=title_fontsize)
    _plt.xticks(fontsize=axes_fontsize[0])
    _plt.yticks(fontsize=axes_fontsize[1])
    _plt.xlabel('X coordinate, a.u.', fontsize=axes_label_fontsize[0])
    _plt.ylabel('Y coordinate, a.u.', fontsize=axes_label_fontsize[1])
    

    x = q_traj[:,:,0].clone().detach()
    y = q_traj[:,:,1].clone().detach()

    for indx, i in enumerate(plot_params["which_timesteps"]):
        _plt.scatter(x[i,:], y[i, :], label=F"{ round( t[i].item()/units.fs2au)}, fs", linewidth=Lw, color=colors[clrs_index[indx]])

    _plt.legend(fontsize=legend_fontsize)
    _plt.tight_layout()


def plot_pop(_plt, t, P, plot_params):
    #plot_params = tsh_plot.common_defaults(plot_params_)

    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    which_energies = plot_params["which_energies"]

    if xlim is not None:
        _plt.xlim(xlim[0], xlim[1])
    if ylim is not None:
        _plt.ylim(ylim[0], ylim[1])

    _plt.title("Transmission probability", fontsize=title_fontsize)
    _plt.xticks(fontsize=axes_fontsize[0])
    _plt.yticks(fontsize=axes_fontsize[1])
    _plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    _plt.ylabel('Probability', fontsize=axes_label_fontsize[1])

    _plt.plot(t[:]/units.fs2au, P.detach(), label=F"{round(P[-1].item(), 4)}", linewidth=Lw, color=colors[clrs_index[0]])

    _plt.legend(fontsize=legend_fontsize)
    _plt.tight_layout()


def plot_energies(_plt, t, E, plot_params):
    #plot_params = tsh_plot.common_defaults(plot_params_)

    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    which_energies = plot_params["which_energies"]

    if xlim is not None:
        _plt.xlim(xlim[0], xlim[1])
    if ylim is not None:
        _plt.ylim(ylim[0], ylim[1])

    _plt.title("Energies", fontsize=title_fontsize)
    _plt.xticks(fontsize=axes_fontsize[0])
    _plt.yticks(fontsize=axes_fontsize[1])
    _plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    _plt.ylabel('Energy, a.u.', fontsize=axes_label_fontsize[1])


    _plt.plot(t[:]/units.fs2au, E[:, 0].detach(), label="Kinetic energy", linewidth=Lw, color=colors[clrs_index[0]] )  # kinetic
    _plt.plot(t[:]/units.fs2au, E[:, 1].detach(), label="Potential energy", linewidth=Lw, color=colors[clrs_index[1]])  # potential
    _plt.plot(t[:]/units.fs2au, E[:, 2].detach(), label="Quantum potential", linewidth=Lw, color=colors[clrs_index[2]])  # quantum
    _plt.plot(t[:]/units.fs2au, E[:, 3].detach(), label="Total energy", linewidth=Lw, color=colors[clrs_index[3]]) # total

    _plt.legend(fontsize=legend_fontsize)
    _plt.tight_layout()


def plot(plot_params_):

    plot_params = tsh_plot.common_defaults(plot_params_)

    filename = plot_params["filename"]
    prefix = plot_params["prefix"]
    output_level = plot_params["output_level"]
    which_dofs = plot_params["which_dofs"]
    which_trajectories = plot_params["which_trajectories"]
    which_adi_states = plot_params["which_adi_states"]
    which_dia_states = plot_params["which_dia_states"]
    what_to_plot = plot_params["what_to_plot"]
    out_prefix = prefix

    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    do_show = plot_params["do_show"]


    #========== Read the data ======
    f = torch.load(F"{filename}")
    t = f["t"]
    q_traj = f["q_traj"]
    P = f["P"]
    E = f["E"]

    #========= Plot or show ========

    plt.figure( figsize=plot_params["figsize"], dpi=plot_params["dpi"],
                       edgecolor='black', frameon=plot_params["frameon"])
    #plt.subplot(1, 1, 1)
    plot_energies(plt, t, E, plot_params)
    if plot_params["save_figures"] == 1:
        plt.savefig(F"{prefix}-energies.png", dpi=plot_params["dpi"])
    if do_show:
        plt.show()
    plt.close()

    plt.figure( figsize=plot_params["figsize"], dpi=plot_params["dpi"],
                       edgecolor='black', frameon=plot_params["frameon"])
    #plt.subplot(1, 1, 1)
    plot_pop(plt, t, P, plot_params)
    if plot_params["save_figures"] == 1:
        plt.savefig(F"{prefix}-transmission.png", dpi=plot_params["dpi"])
    if do_show:
        plt.show()
    plt.close()

    plt.figure(figsize=plot_params["figsize"], dpi=plot_params["dpi"],
                       edgecolor='black', frameon=plot_params["frameon"])
    #plt.subplot(1, 1, 1)
    plot_scatter(plt, t, q_traj, plot_params)
    if plot_params["save_figures"] == 1:
        plt.savefig(F"{prefix}-scatter.png", dpi=plot_params["dpi"])
    if do_show:
        plt.show()

    plt.close()
    #plot_energies(t, E, prefix)
    #plot_pop(t, P, prefix)
    #plot_scatter(T, q_traj, prefix)





