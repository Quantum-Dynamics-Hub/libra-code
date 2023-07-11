#*********************************************************************************
#* Copyright (C) 2021 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: spectrum
   :platform: Unix, Windows
   :synopsis: 
       This module implements functions for computing absorption spectra

.. moduleauthor:: Alexey V. Akimov

"""

import math
import os
import sys
import time
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from . import units, data_conv, pdos
import libra_py.packages.cp2k.methods as CP2K_methods
import util.libutil as comn
import numpy as np


def cp2k_spectrum(_params):
    """Computes various absorption using outputs of the CP2K software

    Args:
        _params (dict) : control parameters of simulations, can contain the following keys:

          * **states** (list of ints): indices of the excited states we want to extract,  [default : [1] ]
              Remember, the index 1 corresponds to the 1-st excited state, but it will be the entry 0 of the
              returned values
              
          * **prefix** ( string ): a common prefix of the filenames for files containing the TDDFPT log files [Required]
          
          * **snapshots** ( list of ints ): integers enumerating the files over which to do the averaging [default: [0] ]

          * **emin** ( double ): the minimal energy of the spectral window [eV, default: 0.0 eV]

          * **emax** ( double ): maximal energy of the spectral window [eV, default: 5.0 eV]

          * **de** ( double ): the original grid spacing of the spectrum  [eV, default: 0.1 eV]

          * **outfile_prefix** ( string ): the prefix of the output file that will contain the final spectrum  [ default: "spectrum"]

          * **do_convolve** ( Bool ): the flag telling whether we want to convolve the original data with the
            Gaussian envelope. The convolution is done with :func:`convolve` [ default: True ] 

          * **de_new** ( double ): the new energy grid spacing, in effect only if do_convolve == True [eV, default: 0.025]

          * var ( double ): standard deviation of the Gaussian with which we do the  convolution
            in effect only if do_convolve == True [eV, default: 0.05] 
                    
    Returns:
        tuple: ( E, spectrum ), where:

            * E ( MATRIX(N, 1) ): new energy grid, N - the new number of energy grid points
            * spectrum ( MATRIX(N, nstates+1) ): spectrum at the grid point for each state and the total
 
    """

    params = dict(_params)

    critical_params = [ "prefix" ] 
    default_params = {"snapshots": [0], "states":[1],
                      "emin":0.0, "emax":5.0, "de":0.1,                        
                      "outfile_prefix":"spectrum", 
                      "do_convolve":True, "de_new":0.025, "var": 0.05
                     }
    comn.check_input(params, default_params, critical_params)

    states = params["states"]
    nstates = len(states)
    prefix = params["prefix"]
    snapshots = params["snapshots"]    
    emin, emax, de = params["emin"], params["emax"], params["de"]
    outfile_prefix = params["outfile_prefix"]
    do_convolve, de_new, var = params["do_convolve"], params["de_new"], params["var"]    


    #============= Dimensions  =================

    N = int(math.floor((emax - emin)/de))+1 # number of the gridpoints
    en   = MATRIX(N, 1)  # energy of the grid points
    spectr = MATRIX(N, nstates+1)  # oscillator strengths for each grid point, and the total one
    
    for i in range(0,N):
        en.set(i, 0, emin + i*de)   # this is a scale centered on Ef

    #============= Data gathering  =================

    nsnaps = float(len(snapshots))
    count = 0.0
    for snap in snapshots:
        
        filename = F"{prefix}{snap}.log"  # file 
        if not os.path.exists(filename):
            print(F"The file {filename} is not found.")             
            
        else:           
            count +=1
            
            E, F = CP2K_methods.cp2k_find_excitation_energies(filename)
        
            for istate, state in enumerate(states):
                
                e = E[state-1]
                f = F[state-1]
                
                state_indx = int(math.floor((e - emin)/de))
        
                spectr.add(state_indx, istate,  f )  # state-resolved spectrum
                spectr.add(state_indx, nstates, f )  # the last entry - the total spectrum

    spectr = spectr / count
                
    #============= Optional convolution =================

    E, SPECTR = MATRIX(en), MATRIX(spectr)

    if do_convolve==True:
        E, SPECTR = pdos.convolve(en, spectr, de, de_new, var)

    #============= Print out ==================

    f = open(outfile_prefix+"_spectrum.txt","w"); 
    
    N = E.num_of_rows
    for i in range(0,N):  # loop over grid points
        line = str(E.get(i,0))+"   "
        tot = 0.0
        for j in range(0,nstates+1):            
            line = line + str(SPECTR.get(i,j))+"   "
        line = line + str(tot)+"\n"        
        f.write(line)
        
    f.close()

    
    return E, SPECTR




def spectrum_plot(plt, E, spectr, _params):
    """
    This function plots spectrum into files and picture for needed projections

    Args:

        plt ( matplotlib instance ): Matplotlib object for plotting

        E ( MATRIX(N, 1) ): the energy grid, N - the new number of energy grid points

        spectr ( MATRIX(N, nstates+1) ): the state-resolved spectra, +1 total spectrum

        _params ( dict ): control parameters for plotting, includes the following keys:

          * **which_projections** ( list of ints): selects which of the available projections we want to plot [ default: [0] ]

          * **labels** ( list of stings ): selects the labels for each data line, should be of the same size as `which_projections` [ default: ["s"] ]
 
          * **colors** ( list of strings or hex codes): defines the colors of the lines to be plotted [ default: ["black"] ]

          * **title** ( string ) : The title of the produced figure [default: "No title"]

          * **output_prefix** ( string ): the name of the file where the spectra will be printed out as a png

    Returns:
        None: but produced .txt files with the needed spectra as well as the pictures of these spectra

    """

    params = dict(_params)

    critical_params = [ ] 
    default_params = {"which_projections": [0], "labels": ["s"],  "colors": ["black"], "title":"No title", "output_prefix":"spectrum_"}
    comn.check_input(params, default_params, critical_params)

    which_projections = params["which_projections"]
    labels = params["labels"]
    colors = params["colors"]
    title = params["title"]
    out_prefix = params["output_prefix"]


    e = data_conv.MATRIX2nparray(E)
    spec = data_conv.MATRIX2nparray(spectr)

    nplots = len(which_projections)

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
   
    # Plot the total density by black color
    for i in range(nplots):
        plt.plot(e[:, 0], spec[:, which_projections[i]], label=labels[i],color=colors[i], linewidth=2.0)

    plt.title(title,fontsize=12)
    plt.legend(fontsize=6.75, ncol=1, loc='upper center')
    plt.xlabel('Energy, eV',fontsize=12)
    plt.ylabel('Absorbance',fontsize=12)

    plt.tight_layout()
    plt.savefig(F'{out_prefix}.png', dpi=600)
    plt.show()



