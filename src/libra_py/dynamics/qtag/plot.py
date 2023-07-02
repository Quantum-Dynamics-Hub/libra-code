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

from liblibra_core import *
import util.libutil as comn

from libra_py import data_read


def plot_wf_1D(_plt_params, data):
    """Plots the 1D wavefunction

    Args:
        _plt_params (dict): Dictionary containing plot control parameters.

          * **_plt_params[`which_states`]** (list of ints) : the list containing surfaces for which data is to be plotted 
              [default : [0] ] 

          * **_plt_params[`xmin`]** ( list: [float] ) : the domain minimum for the plot [ default: [-4.0] ]
 
          * **_plt_params[`xmax`]** ( list: [float] ) : the domain maximum for the plot [ default: [4.0] ]
 
          * **_plt_params[`snaps`]** (list of ints) : for which timesteps to plot the wavefunction [ default: [0] ]

          * **_plt_params[`size`]** ( list: [int, int] ): size of the plot [ default: [16,16] ]

          * **_plt_params[`xlabel`]** ( string ): x label of the plot [ default: "Position (a.u.)"]

          * **_plt_params[`ylabel`]** ( string ): y label of the plot [ default: "Density"]

          * **_plt_params[`1Dcolors`]** ( list of strings ): color names for each surface, should be of the same shape as
            `which_states` [ default: ["Black"] ]

          * **_plt_params[`legend_loc`]** ( [float, float]): position of the legend [ default: [0.6, 0.8] ]

          * **_plt_params[`legend_size`]** ( int ) : font size of the legend [ default: 12 ] 

        data (list of lists) : the data to be plotted, extracted from plots.wf_plot()
    """

    plt_params = dict(_plt_params)
    critical_params = []
    default_params = {'which_states':[0], 'xmin':[-4.0],'xmax':[4.0],  'snaps':[0],
                      'size':[16, 16], 'xlabel':'Position (a.u.)', 'ylabel':'Density',
                      '1Dcolors':['Black'],'legend_loc':[0.6,0.8], 'legend_size':12
                     }
    comn.check_input(plt_params, default_params, critical_params)

    fig = plt.figure(figsize=(plt_params['size'][0], plt_params['size'][1]))
    snaps = plt_params["snaps"]
    nsnaps = len(snaps)
    colors = plt_params['1Dcolors']
    m = int(np.ceil(nsnaps/3.0))
    ndof = 1

    which_states = plt_params['which_states']
    labels = []
    for state in which_states:
        labels.append(F"State {state}")


    for i in range(len(snaps)):
        ax = fig.add_subplot(m,3,i+1)
        for j in range(len(which_states)):
            x, y = data[i][0], data[i][j+ndof]
            ax.plot(x, y, color=colors[j], label=labels[j])

            ax.set_xlabel(plt_params['xlabel'])
            ax.set_ylabel(plt_params['ylabel'])

            ax.axes.set_xlim(plt_params['xmin'][0], plt_params['xmax'][0] )
            ax.legend(loc = (plt_params['legend_loc'][0],plt_params['legend_loc'][1]), prop={'size' : plt_params['legend_size']})

    plt.tight_layout()


def plot_wf_2D(_plt_params, data):
    """Plots the 2D wavefunction

    Args:
        _plt_params (dict): Dictionary containing plot control parameters.

          * **_plt_params[`which_states`]** (list of ints) : the list containing surfaces for which data is to be plotted
              [default : [0] ]

          * **_plt_params[`xmin`]** ( list: [float, float] ) : the domain minimum for the plot [ default: [-4.0, -4.0] ]

          * **_plt_params[`xmax`]** ( list: [float, float] ) : the domain maximum for the plot [ default: [4.0, 4.0] ]

          * **_plt_params[`snaps`]** (list of ints) : for which timesteps to plot the wavefunction [ default: [0] ]

          * **_plt_params[`size`]** ( list: [int, int] ): size of the plot [ default: [16,16] ]

          * **_plt_params[`xlabel`]** ( string ): x label of the plot [ default: "X Coordinate (a.u.)"]

          * **_plt_params[`ylabel`]** ( string ): y label of the plot [ default: "Y Coordinate (a.u.)"]

          * **_plt_params[`zlabel`]** ( string ): z label of the plot [ default: "Density"]

          * **_plt_params[`2Dcolors`]** ( list of strings ): color names for each surface, should be of the same shape as
            `which_states` [ default: ["Reds"] ]

          * **_plt_params[`legend_loc`]** ( [float, float]): position of the legend [ default: [0.6, 0.8] ]

          * **_plt_params[`legend_size`]** ( int ) : font size of the legend [ default: 12 ]

        data (list of lists) : the data to be plotted, extracted from plots.wf_plot()

    """

    plt_params = dict(_plt_params)
    critical_params = []
    default_params = {'which_states':[0], 'xmin':[-4.0, -4.0],'xmax':[4.0, 4.0],  'snaps':[0],
                      'xlabel':'X Coordinate (a.u.)', 'ylabel':'Y Coordinate (a.u.)', 'zlabel':'Density',
                      '2Dcolors':['Reds'],'legend_loc':[0.6,0.8], 'legend_size':12
                     }
    comn.check_input(plt_params, default_params, critical_params)


    fig = plt.figure(figsize=plt_params['size'])
    snaps = plt_params["snaps"]
    nsnaps = len(snaps) # used to be len(data)
#    npoints = plt_params['npoints']
#    states = sorted(plt_params['states'])
    which_states = plt_params['which_states']
    colors = plt_params['2Dcolors']
    m = int(np.ceil(nsnaps/3.0))
    ndof = 2

    for i in range(len(snaps)):
        for j in range(len(which_states)):
            x, y, z = data[i][0], data[i][1], data[i][j+ndof]
            ax = fig.add_subplot(m, 3, i+1, projection = '3d')
            ax.contour3D(x, y, z, 33, cmap=colors[j])  # used to be npoints[0] instead of 33

            ax.set_xlabel(plt_params['xlabel'])
            ax.set_ylabel(plt_params['ylabel'])
            ax.set_zlabel(plt_params['zlabel'])

            ax.axes.set_xlim3d(plt_params['xmin'][0],plt_params['xmax'][0])
            ax.axes.set_ylim3d(plt_params['xmin'][1],plt_params['xmax'][1])



def wf_plot(_plt_params):
    """Plots 1D or 2D wavefunctions

    Args:
        _plt_params (dict): Dictionary containing plot control parameters.

          * **_plt_params[`ndof`]** (int) : the number of nuclear degrees of freedom in data [ default : 1 ]

          * **_plt_params[`states`]** (list of ints) : the states present in the output data [ default: [0] ]

          * **_plt_params[`prefix`]** (string) : the name of the subdirectory where QTAG output is stored [ default: "out" ]

          * **_plt_params[`npoints`]** ( list of ints) : how many points to use in the plotting [ default: [100] for 1D, and [25, 25] for 2D ]

    Also see: see docs for `plot_wf_1D` or `plot_wf_2D` for the following parameters:

          * **_plt_params[`which_states`]**

          * **_plt_params[`xmin`]**

          * **_plt_params[`xmax`]**

          * **_plt_params[`snaps`]**

          * **_plt_params[`size`]**

          * **_plt_params[`xlabel`]**

          * **_plt_params[`ylabel`]**

          * **_plt_params[`zlabel`]**

          * **_plt_params[`1Dcolors`]**

          * **_plt_params[`2Dcolors`]**

          * **_plt_params[`legend_loc`]**

          * **_plt_params[`legend_size`]**

        data (list of lists) : the data to be plotted, extracted from plots.wf_plot()

    """

    plt_params = dict(_plt_params)
    critical_params = []
    default_params = {'states':[0], "which_states":[0],
                      "ndof":1, "prefix":"out", "snaps":[0],
                     }
    comn.check_input(plt_params, default_params, critical_params)

    npts_1D = [100]
    npts_2D = [25,25]
    ndof = plt_params['ndof']
    if ndof==1:
      default_params.update({"npoints": npts_1D })
    elif ndof==2:
      default_params.update({"npoints": npts_2D })
    comn.check_input(plt_params, default_params, critical_params)


    prefix = plt_params['prefix']
    states = sorted(plt_params['states'])
    which_states = sorted(plt_params['which_states'])
    snaps = plt_params['snaps']
    npoints = plt_params['npoints']

    which_cols = [dof for dof in range(ndof)]
    for state in which_states:
        which_cols.append(ndof + states.index(state))

    filenames = []
    for snap in snaps:
        filenames.append(F"{prefix}/wfc/wfcr_snap_{snap}_dens_rep_0")

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



def energy_and_pops(_plt_params):
    """

    Args:
        _plt_params (dict): Dictionary containing simulation parameters.

          * **_plt_params[`nsteps`]** (int) : the number of simulation steps [ default: 1 ]

          * **_plt_params[`states`]** (list of ints) : the states present in the output data [ default: [0] ]

          * **_plt_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]

          * **_plt_params[`prefix`]** (string) : the name of the subdirectory where QTAG output is stored

          * **_plt_params[`xmin`]** (list of floats) : the list of domain minima for the plot [ default : [0.0] ]

          * **_plt_params[`xmax`]** (list of floats) : the list of domain maxima for the plot [ default: [1.0] ]

          * **_plt_params[`1Dcolors`]** ( list of strings ): color names for each surface, should be of the same shape as
            `which_states` [ default: ["Black"] ]

          * **_plt_params[`legend_loc`]** ( [float, float]): position of the legend [ default: [0.6, 0.8] ]

          * **_plt_params[`legend_size`]** ( int ) : font size of the legend [ default: 12 ]

    """

    plt_params = dict(_plt_params)
    critical_params = []
    default_params = {"nsteps":1, 'states':[0], "ndof":1, 
                      "prefix":"out", "xmin":[0.0], "xmax":[1.0], "1Dcolors":["Black"],
                      "legend_loc":[0.6, 0.8], "legend_size":12
                     }
    comn.check_input(plt_params, default_params, critical_params)


    prefix = plt_params['prefix']
    colors = plt_params['1Dcolors']
    states = sorted(plt_params['states'])
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



def trajectories(_plt_params):
    """

    Args:
        _plt_params (dict): Dictionary containing simulation parameters.

          * **_plt_params["data_type"]** (int): 0 - coordinates, 1 - momenta, 2 - width parameters [default: 0]

          * **_plt_params[`grid_dims`]** (list of ints) : the *grid_dims* list used to specify the initial basis passed to initialize.py

          * **_plt_params[`states`]** (list of ints) : the states present in the output data [ default: [0] ]

          * **_plt_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]

          * **_plt_params[`prefix`]** (str) : the name of the subdirectory where QTAG output is stored

          * **_plt_params[`which_states`]** (list of ints) : the list containing surfaces for which data is to be plotted

          * **_plt_params["which_traj"]** (list of ints or string) : the list of trajectories to print

          * **_plt_params[`xmin`]** (list of floats) : the list of domain minima for the plot

          * **_plt_params[`xmax`]** (list of floats) : the list of domain maxima for the plot

          * **_plt_params[`active_dof`]** (int) : index of the active DOF to plot. [defualt : 0 ]

             The format of the q, p, etc. files is: 
             q0 (t,    traj=0, state=0) q1(t,    traj=0, state=0) ... q0(t,    traj=ntraj-1, state0) .. then state 1  state 2 ,etc
             q0 (t+dt, traj=0, state=0) q1(t+dt, traj=0, state=0) ... q0(t+dt, traj=ntraj-1, state0) .. then state 1  state 2 ,etc

    """

    plt_params = dict(_plt_params)
    critical_params = [ "grid_dims" ]
    default_params = { "data_type":0,
                       'states':[0], "ndof":1, "which_states":[0],"active_dof":0,
                       "prefix":"out", "xmin":[0.0], "xmax":[1.0], "1Dcolors":["Black"],
                       "legend_loc":[0.6, 0.8], "legend_size":12
                     }
    comn.check_input(plt_params, default_params, critical_params)

    data_type = plt_params["data_type"]
    prefix = plt_params['prefix']
    colors = plt_params['1Dcolors']
    states = sorted(plt_params['states'])
    which_states = plt_params['which_states']
    nstates = len(states)
    ndof = plt_params['ndof']
    active_dof = plt_params["active_dof"]

    which_states = plt_params['which_states']
    if type(which_states) == str:
        which_states == which_states.lower()
    elif type(which_states) == list:
        which_states = sorted(which_states)
    else:
        sys.exit('Error in which_states parameter in plt_params dict! Type should be list or "all"')
    if which_states == 'all':
        which_states = list(range(nstates))

    ntraj = 1
    for n in range(ndof):
        ntraj *= plt_params["grid_dims"][n]

    tfile = f"{prefix}/time.txt"

    qfile = ""
    if data_type==0:
        qfile = F"{prefix}/q.txt"
    elif data_type==1:
        qfile = F"{prefix}/p.txt"
    elif data_type==2:
        qfile = F"{prefix}/a.txt"


    tdata = data_read.get_data_from_file2(tfile, [0])
    qdata = data_read.get_data_from_file2(qfile, list(range(ntraj*nstates*ndof)) )

    nsteps = len(tdata[0])
    qnew = np.zeros((nstates,ntraj, ndof, nsteps),dtype=float)

    for state in range(nstates):
        for traj in range(ntraj):
            for idof in range(ndof):
                qnew[state][traj][idof] = qdata[ idof + ndof*(traj+ntraj*state)]

    which_traj = plt_params['which_traj']
    if type(which_traj) == str:
        which_traj = which_traj.lower()
    elif type(which_traj) == list:
        which_traj = sorted(which_traj)
    else:
        sys.exit('Error in which_traj parameter in plt_params dict! Type should be list or "all"')
    if which_traj == 'all':
        which_traj = list( range(ntraj) )


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
            x, y = tdata[0], qnew[state][traj][active_dof]

            ax.plot(x, y, color=colors[state])

        ax.set_xlabel(plt_params['xlabel'])
        ax.set_ylabel(plt_params['ylabel'])

        ax.axes.set_xlim(xlimits[0],xlimits[1])
        ax.axes.set_ylim(ylimits[0],ylimits[1])
        ax.text(6,6.5,labels[state],fontsize=12)

    plt.tight_layout()


