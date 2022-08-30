#*********************************************************************************                     
#* Copyright (C) 2022 Matthew Dutra and Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
..module:: plot
  :platform: Unix, Windows
  :synopsis: This module contains functions for plotting various qtag outputs.

..moduleauthors :: Matthew Dutra
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from libra_py import data_read

def plot_wf_1D(plt_params, data, *ref_data):
    """
        plt_params (dict): Dictionary containing plot control parameters.

          * **plt_params[`which_states`]** (list of ints) : the list containing surfaces for which data is to be plotted

          * **plt_params[`xmin`]** (list of floats) : the list of domain minima for the plot

          * **plt_params[`xmax`]** (list of floats) : the list of domain maxima for the plot

        data (list of lists) : the data to be plotted, extracted from plots.wf_plot()
    """

    fig = plt.figure(figsize=plt_params['size'])
    nsnaps = len(data)
    states = sorted(plt_params['states'])
    which_states = plt_params['which_states']
    colors = plt_params['1Dcolors']
    m = int(np.ceil(nsnaps/3.0))
    ndof = 1

    xlimits = [plt_params['xmin'][0],plt_params['xmax'][0]]

    labels = []
    for state in which_states:
        labels.append("State "+str(state))

    for i in range(nsnaps):
        ax = fig.add_subplot(m,3,i+1)
        for j in range(len(which_states)):
            x, y = data[i][0], data[i][j+ndof]
            ax.plot(x, y, color=colors[j], label=labels[j])

            ax.set_xlabel(plt_params['xlabel'])
            ax.set_ylabel(plt_params['ylabel'])

            ax.axes.set_xlim(xlimits[0],xlimits[1])
            ax.legend(loc = plt_params['legend_loc'], prop={'size' : plt_params['legend_size']})

#            if ref_data:
#                ax1.plot(ref_data[0], ref_data[1+j], 'o', color=colors[j],markersize=2)
#            ax1.plot(data[0], data[1+j],label=F"{labels[j]}")

    plt.tight_layout()

#    plt.subplots()

def plot_wf_2D(plt_params, data, *ref_data):
    """
        plt_params (dict): Dictionary containing plot control parameters.

          * **plt_params[`which_states`]** (list of ints) : the list containing surfaces for which data is to be plotted

          * **plt_params[`npoints`]** (list of ints) : the grid used to compute the wavefunction in compute.wf_calc_nD()

          * **plt_params[`xmin`]** (list of floats) : the list of domain minima for the plot

          * **plt_params[`xmax`]** (list of floats) : the list of domain maxima for the plot

        data (list of lists) : the data to be plotted, extracted from plots.wf_plot()
    """

    fig = plt.figure(figsize=plt_params['size'])
    nsnaps = len(data)
    npoints = plt_params['npoints']
    states = sorted(plt_params['states'])
    which_states = plt_params['which_states']
    colors = plt_params['2Dcolors']
    m = int(np.ceil(nsnaps/3.0))
    ndof = 2

    xlimits = [plt_params['xmin'][0],plt_params['xmax'][0]]
    ylimits = [plt_params['xmin'][1],plt_params['xmax'][1]]

    for i in range(nsnaps):
        for j in range(len(which_states)):
            x, y, z = data[i][0], data[i][1], data[i][j+ndof]
            ax = fig.add_subplot(m, 3, i+1, projection = '3d')
            ax.contour3D(x, y, z, npoints[0], cmap=colors[j])

            ax.set_xlabel(plt_params['xlabel'])
            ax.set_ylabel(plt_params['ylabel'])
            ax.set_zlabel(plt_params['zlabel'])

            ax.axes.set_xlim3d(xlimits[0],xlimits[1])
            ax.axes.set_ylim3d(ylimits[0],ylimits[1])

def wf_plot(dyn_params, plt_params):
    """

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`states`]** (list of ints) : the states present in the output data [ default: 2 ]

          * **dyn_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]

        plt_params (dict): Dictionary containing plot control parameters.

          * **plt_params[`prefix`]** (str) : the name of the subdirectory where QTAG output is stored

          * **plt_params[`which_states`]** (list of ints) : the list containing surfaces for which data is to be plotted

          * **plt_params[`snaps`]** (list of ints) : the list of snapshots used in compute.wf_calc_nD()

          * **plt_params[`npoints`]** (list of ints) : the grid used to compute the wavefunction in compute.wf_calc_nD()
    """

    ndof = dyn_params['ndof']

    prefix = plt_params['prefix']
    states = sorted(dyn_params['states'])
    which_states = sorted(plt_params['which_states'])
    snaps = plt_params['snaps']
    npoints = plt_params['npoints']

    which_cols = [dof for dof in range(ndof)]
    for state in which_states:
        which_cols.append(ndof+states.index(state))

    filenames = []
    for snap in snaps:
        filenames.append(f"{prefix}/wfc/wfcr_snap_"+str(snap)+"_dens_rep_0")

    data = []
    for i in range(len(filenames)):
        data.append(data_read.get_data_from_file2(filenames[i], which_cols))
        if ndof == 2:
            for j in range(len(which_cols)):
                data[i][j] = np.reshape(data[i][j],tuple(npoints))

    plt_params['states'] = states
    print("Number of files = " + str(len(filenames)))
    if ndof == 1:
        plot_wf_1D(plt_params, data)
    elif ndof == 2:
        plot_wf_2D(plt_params, data)
    else:
        print("Plotting for ndof > 2 not yet supported!")

def energy_and_pops(dyn_params, plt_params):
    """

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`nsteps`]** (int) : the number of simulation steps

          * **dyn_params[`states`]** (list of ints) : the states present in the output data [ default: 2 ]

          * **dyn_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]

        plt_params (dict): Dictionary containing plot control parameters.

          * **plt_params[`prefix`]** (str) : the name of the subdirectory where QTAG output is stored

          * **plt_params[`xmin`]** (list of floats) : the list of domain minima for the plot

          * **plt_params[`xmax`]** (list of floats) : the list of domain maxima for the plot
    """

    prefix = plt_params['prefix']
    colors = plt_params['1Dcolors']
    states = sorted(dyn_params['states'])
    nstates = len(states)

    tfile = f"{prefix}/time.txt"
    popsfile = f"{prefix}/pops.txt"
    efile = f"{prefix}/Etot.txt"

    tdata = data_read.get_data_from_file2(tfile, [0])
    popdata = data_read.get_data_from_file2(popsfile, states)
    edata = data_read.get_data_from_file2(efile, [0])

    e0 = edata[0][0]
    nsteps = len(tdata[0])
    totalpop=np.zeros(nsteps)

    for n in range(nsteps):
        edata[0][n]=(edata[0][n]-e0)/e0*100
        for state in range(nstates):
            totalpop[n]+=popdata[state][n]

    fig = plt.figure(figsize=plt_params['size'])

    xlimits = [plt_params['xmin'][0],plt_params['xmax'][0]]

    labels = []
    for state in states:
        labels.append("State "+str(state))

    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    for j in range(nstates):
        ax1.plot(tdata[0],popdata[j],color=colors[j],label=labels[j])
    ax1.plot(tdata[0],totalpop,color='black',label='Total')

    ax2.plot(tdata[0],edata[0],color='green')

    ax1.set_xlabel("Time (a.u.)")
    ax1.set_ylabel("Populations")
#    ax1.tick_params(labelsize=15)

    ax2.set_xlabel("Time (a.u.)")
    ax2.set_ylabel("% Error, Energy")
#    ax2.tick_params(labelsize=15)

    ax1.axes.set_xlim(xlimits[0],xlimits[1])
    ax2.axes.set_xlim(xlimits[0],xlimits[1])

    ax1.legend(loc = plt_params['legend_loc'], prop={'size' : plt_params['legend_size']})

    plt.tight_layout()



def trajectories(dyn_params, plt_params):
    """

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`grid_dims`]** (list of ints) : the *grid_dims* list used to specify the initial basis passed to initialize.py

          * **dyn_params[`states`]** (list of ints) : the states present in the output data [ default: 2 ]

          * **dyn_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]

        plt_params (dict): Dictionary containing plot control parameters.

          * **plt_params[`prefix`]** (str) : the name of the subdirectory where QTAG output is stored

          * **plt_params[`which_states`]** (list of ints) : the list containing surfaces for which data is to be plotted

          * **plt_params[`xmin`]** (list of floats) : the list of domain minima for the plot

          * **plt_params[`xmax`]** (list of floats) : the list of domain maxima for the plot

    """

    prefix = plt_params['prefix']
    colors = plt_params['1Dcolors']
    states = sorted(dyn_params['states'])
    which_states = plt_params['which_states']
    nstates = len(states)
    ndof = dyn_params['ndof']

    which_states = plt_params['which_states']
    if type(which_states) == str:
        which_states == which_states.lower()
    elif type(which_states) == list:
        which_states = sorted(which_states)
    else:
        sys.exit('Error in which_states parameter in plt_params dict! Type should be list or "all"')
    if which_states == 'all':
        which_states = [state for state in range(nstates)]

    ntraj = 1
    for n in range(ndof):
        ntraj *= dyn_params["grid_dims"][n]

    tfile = f"{prefix}/time.txt"
    qfile = f"{prefix}/q.txt"

    tdata = data_read.get_data_from_file2(tfile, [0])
    qdata = data_read.get_data_from_file2(qfile, [traj for traj in range(ntraj*nstates)])

    nsteps = len(tdata[0])
    qnew = np.zeros((nstates,ntraj,nsteps),dtype=float)

    for state in range(nstates):
        for traj in range(ntraj):
            qnew[state][traj] = qdata[traj+ntraj*state]

    which_traj = plt_params['which_traj']
    if type(which_traj) == str:
        which_traj = which_traj.lower()
    elif type(which_traj) == list:
        which_traj = sorted(which_traj)
    else:
        sys.exit('Error in which_traj parameter in plt_params dict! Type should be list or "all"')
    if which_traj == 'all':
        which_traj = [traj for traj in range(ntraj)]

    labels = []
    for state in which_states:
        labels.append("State "+str(state))

    fig = plt.figure(figsize=plt_params['size'])
    xlimits = [plt_params['xmin'][0],plt_params['xmax'][0]]
    ylimits = [plt_params['xmin'][1],plt_params['xmax'][1]]
    m = int(np.ceil(len(which_states)/3.0))    
    for state in which_states:
        ax = fig.add_subplot(m,3,state+1)
        for traj in which_traj:
            x, y = tdata[0], qnew[state][traj]

            ax.plot(x, y, color=colors[state])

        ax.set_xlabel(plt_params['xlabel'])
        ax.set_ylabel(plt_params['ylabel'])

        ax.axes.set_xlim(xlimits[0],xlimits[1])
        ax.axes.set_ylim(ylimits[0],ylimits[1])
        ax.text(6,6.5,labels[state],fontsize=12)

    plt.tight_layout()
