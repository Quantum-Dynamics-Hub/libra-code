#*********************************************************************************
#* Copyright (C) 2018-2019 Wei Li and Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""

.. module:: qsh
   :platform: Unix, Windows
   :synopsis:   Implementation of the QuasiStochastic Hamiltonian method
       Akimov, J. Phys. Chem. Lett. 2017, 8, 5190

.. moduleauthor:: Wei Li and Alexey V. Akimov


"""


import cmath
import math
import os
import sys
import json

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


#import libra_py.common_utils as comn
import util.libutil as comn

import libra_py.units as units
import libra_py.influence_spectrum as influence_spectrum
import libra_py.data_stat as data_stat
import libra_py.data_conv as data_conv


def compute_freqs(H_vib, params):
    """Compute a matrix of frequencies for each matrix element

    Args:

        H_vib ( list of CMATRIX objects ): the vibronic Hamiltonian for all time-points
            H_vib[istep].get(i,j) - i,j matrix element for the step `istep`
    
        params ( dictionary ): the parameters that control the execution of this function

                * nsteps: the number of the QSH-MD steps to do

            SeeAlso: influence_spectrum.compute_mat_elt for the description of other parameters:

                * filename
                * logname
                * nfreqs

            SeeAlso: influence_spectrum.recipe1(data, params) for the description of other parameters: 
            
                * dt
                * wspan
                * dw
                * do_output
                * acf_filename
                * spectrum_filename
                * do_center
                * acf_type
                * data_type

    Returns:
        tuple: (freqs, dev), where:

            SeeAlso: influence_spectrum.compute_mat_elt for the description of the output `freqs`

            dev ( list of lists of double ): 
                dev[i][j] is the standard deviation the QSH terms (sum over frequencies)

    """

    # Set defaults and check critical parameters
    critical_params = [ ] 
    default_params = { "output_files_prefix":"_qsh", "nsteps":len(H_vib)     }
    comn.check_input(params, default_params, critical_params)


    output_files_prefix = params["output_files_prefix"]


    # Local variables and dimensions    
    nstates = H_vib[0].num_of_rows

    freqs, T,  norm_acf,  raw_acf,  W,  J, J2 = influence_spectrum.compute_all(H_vib, params)  # T in fs, W and freqs in cm^-1



    ## ============ Plot the spectra =============
    nsteps = params["nsteps"]
    nstates = len(freqs)

    for i in range(nstates):
        for j in range(i, nstates):

            figure = plt.figure(num=None, figsize=(3.21*2, 2.41), dpi=300, edgecolor='black', frameon=True)
            
            plt.subplot(1,2,1)
            plt.xlabel('Time, fs',fontsize=10)
            plt.ylabel('normalized ACF',fontsize=10)
            plt.plot(T[i][j], norm_acf[i][j], color="blue", label="", linewidth=2)
            
            plt.subplot(1,2,2)
            plt.xlabel('Wavenumber, $cm^{-1}$',fontsize=10)
            plt.ylabel('Influence spectrum',fontsize=10)
            plt.plot(W[i][j], J2[i][j], color="blue", label=F"", linewidth=2)
            
            plt.tight_layout()
            plt.savefig(F"{output_files_prefix}-{i}-{j}-acf-ifs.png", dpi=300)
            plt.show()




    dt = params["dt"]
    conv = units.wavn2au * units.fs2au

    
    # Ok, now we have the function - sum of sines, so let's compute the standard deviation
    # This is a silly method - just do it numerically
    dev = [ [ 0.0 for i in range(0,nstates)] for j in range(0,nstates)]
        
    for i in range(0,nstates):
        for j in range(0,nstates):

            # Adjust the number of frequencies
            nfreqs = params["nfreqs"]            
            if len(freqs[i][j]) < nfreqs:
                nfreqs = len(freqs[i][j])

            # Compute the variation of the QSH terms
            fu_ave, fu2_ave = 0.0, 0.0
            for r in range(0, nsteps):

                fu = 0.0
                for k in range(0,nfreqs):
                    #fu = fu + freqs[i][j][k][2] * math.sin(conv*freqs[i][j][k][0]*r*dt)
                    fu = fu + freqs[i][j][k][3] * math.sin(conv*freqs[i][j][k][0]*r*dt)

                fu_ave = fu_ave + fu
                fu2_ave = fu2_ave + fu*fu

            fu_ave = fu_ave / nsteps
            fu2_ave = fu2_ave / nsteps

            dev[i][j] = math.sqrt( (fu2_ave - fu_ave**2) ) 	
    
    return freqs, dev




def compute_qs_Hvib(Nfreqs, freqs, t, 
                 H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                 H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                 dev):
    """Compute the QSH Hamiltonians

    Args:
        Nfreqs ( int ): the number of frequencies we want to use in the QSH calculations (upper limit, the actual number could be smaller)
        freqs ( list of lists of doubles ): contains the spectral info for various matrix elements of the sampled Hvib
            SeeAlso: influence_spectrum.compute_mat_elt for the description of the output `freqs`

        t ( double ): time at which we want to reconstruct the QSH [ units: a.u. ]
        H_vib_re_ave ( MATRIX(nstates, nstates) ): average of energies for direct Hamiltonian
        H_vib_im_ave ( MATRIX(nstates, nstates) ): average of couplings for direct Hamiltonian
        H_vib_re_std ( MATRIX(nstates, nstates) ): std of energies for direct Hamiltonian
        H_vib_im_std ( MATRIX(nstates, nstates) ): std of couplings for direct Hamiltonian
        up_Hvib_re ( MATRIX(nstates, nstates) ): maximum value of energies for direct Hamiltonian
        up_Hvib_im ( MATRIX(nstates, nstates) ): maximum value of couplings for direct Hamiltonian
        dw_Hvib_re ( MATRIX(nstates, nstates) ): minimal value of energies for direct Hamiltonian
        dw_Hvib_im ( MATRIX(nstates, nstates) ): minimal value of couplings for direct Hamiltonian
        dev ( list of lists of doubles, nstates x nstates ): std for direct Hamiltonian
            SeeAlso: compute_freqs for the description of the output `dev`
    
    Returns: 
        CMATRIX(nstates, nstates): contains QSH vibronic Hamiltonian at a given time

    """

    nstates = H_vib_re_ave.num_of_cols
    conv = units.wavn2au
    
    Hvib_stoch_re = MATRIX(nstates,nstates)
    Hvib_stoch_im = MATRIX(nstates,nstates)

    fu = [ [ 0.0 for i in range(0,nstates)] for j in range(0,nstates)]

    for i in range(0,nstates):
        for j in range(0,nstates):
            # Adjust the number of frequencies
            nfreqs = Nfreqs
            if len(freqs[i][j]) < Nfreqs:
                nfreqs = len(freqs[i][j])

            for k in range(0,nfreqs):
                #fu[i][j] = fu[i][j] + freqs[i][j][k][2] * math.sin(conv*freqs[i][j][k][0]*t)
                fu[i][j] = fu[i][j] + freqs[i][j][k][3] * math.sin(conv*freqs[i][j][k][0]*t)


    for i in range(0,nstates):
        for j in range(0,nstates):
            if i==j:
                xab = 0.0
                if dev[i][j]>0.0:
                    xab = H_vib_re_ave.get(i,j) + H_vib_re_std.get(i,j) * (fu[i][j]/dev[i][j] ) 

                if xab < dw_Hvib_re.get(i,j):
                    xab = dw_Hvib_re.get(i,j)
                elif xab > up_Hvib_re.get(i,j):
                    xab = up_Hvib_re.get(i,j)
                Hvib_stoch_re.set(i,j,   xab )

            elif i<j:
                xab = 0.0
                if dev[i][j]>0.0:
                    xab = H_vib_im_ave.get(i,j) + H_vib_im_std.get(i,j) * (fu[i][j]/dev[i][j] ) 

                if xab < dw_Hvib_im.get(i,j):
                    xab = dw_Hvib_im.get(i,j)
                elif xab > up_Hvib_im.get(i,j):
                    xab = up_Hvib_im.get(i,j)
                Hvib_stoch_im.set(i,j,  -xab )
                Hvib_stoch_im.set(j,i,   xab )

    Hvib_stoch = CMATRIX(Hvib_stoch_re, Hvib_stoch_im)

    return Hvib_stoch




def run(H_vib, params):
    """

    The procedure to convert the results of QE/model Hvib calculations to longer 
    timescales using the QSH approach

    Args:
        H_vib ( list of lists of CMATRIX ): the vibronic Hamiltonian for all data sets and all time-points

            Such that H_vib[idata][istep].get(i,j) is the i,j matrix element for the data 
            set ```idata``` and step in that data set ```istep```

        params ( dictionary ): controls the present and all underlying level calculations

            * **params["dt"]** ( double ): nuclear dynamics integration time step
                this is also the spacing between initial data points [ units: a.u. of time, default: 41.0]

            * **params["nfreqs"]** ( int ): maximal number of frequencies to use to reconstruct all the matrix elements [ default: 1]

            * **params["nsteps"]** ( int ): how many time-steps of QSH to generate [ default: 10]

            * **params["do_QSH_output"]** ( Boolean ): the flag that determines whether to generate the 
                output files containing the QSH Hamiltonians. [ default : False] 

            * **params["output_set_paths"]** ( list of strings ): 
                Directories where the resulting QSH Hvib files will be printed out. These directories should already exist
                [default: QSH_#_of_data_set, e.g. QSH_0, QSH_1, etc.]

            * **params["qsh_Hvib_re_prefix"]** ( string ): prefixes of the output files with real part 
                of the QSH vibronic Hamiltonian at time t

            * **params["qsh_Hvib_re_suffix"]** ( string ): suffixes of the output files with real part 
                of the QSH vibronic Hamiltonian at time t

            * **params["qsh_Hvib_im_prefix"]** ( string ): prefixes of the output files with imaginary part
                of the QSH vibronic Hamiltonian at time t

            * **params["qsh_Hvib_im_suffix"]** ( string ): suffixes of the output files with imaginary part
                of the QSH vibronic Hamiltonian at time t

    Returns: 
        ( list of lists ): qsh_H_vib, such as
            qsh_H_vib[idata][istep] - is a CMATRIX(nstates, nstates) object representing a vibronic Hamiltonian
            predicted for the time step `istep` using the training dataset `idata`
     
    """

    # General dimensions
    ndata = len(H_vib)
    nstates = H_vib[0][0].num_of_cols
    output_set_paths = []
    for i in range(0,ndata): 
        output_set_paths.append( "QSH_%s" % (i) )

    # Check the input parameters and setup defaults
    critical_params = [ ] 
    default_params = { "dt":41.0, "nfreqs":1, "nsteps":10,
                       "do_QSH_output":False,
                       "output_set_paths":output_set_paths,
                       "qsh_Hvib_re_prefix":"qsh_Hvib_", "qsh_Hvib_im_prefix":"qsh_Hvib_",
                       "qsh_Hvib_re_suffix":"_re",   "qsh_Hvib_im_suffix":"_im",
                       "meta_info":"qsh-meta" }
    comn.check_input(params, default_params, critical_params)

    # Local parameters
    nfreqs = params["nfreqs"]
    nsteps = params["nsteps"]
    dt = params["dt"]
    meta_info_file = params["meta_info"]

    # Parameters for the undelying functions
    params1 = dict(params)
    params1["dt"] = params["dt"] * units.au2fs  # compute_freqs takes dt in fs


    qsh_H_vib = []

    for idata in range(0,ndata):   # over all MD trajectories (data sets)
        
        #======== Split complex-valued matrices into real and imaginary sets ============

        H_vib_re = []  # list of MATRIX
        H_vib_im = []  # list of MATRIX

        nsteps0 = len(H_vib[idata])   
        for i in range(0,nsteps0):
            H_vib_re.append(H_vib[idata][i].real())
            H_vib_im.append(H_vib[idata][i].imag())


        #======== Analyze the Hvib time-seris  ============
        H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re = data_stat.mat_stat(H_vib_re)
        H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im = data_stat.mat_stat(H_vib_im)    

        
        freqs, dev = compute_freqs(H_vib[idata], params1)   # freqs are in cm^-1



        ## ============= Now save meta-information and the data ==========
        meta_info = { "freqs": freqs, "dev":dev, 
                      "H_vib_re_ave":data_conv.matrix2list(H_vib_re_ave),
                      "H_vib_im_ave":data_conv.matrix2list(H_vib_im_ave),
                      "H_vib_re_std":data_conv.matrix2list(H_vib_re_std),
                      "H_vib_im_std":data_conv.matrix2list(H_vib_im_std),
                      "dw_Hvib_re":data_conv.matrix2list(dw_Hvib_re),
                      "dw_Hvib_im":data_conv.matrix2list(dw_Hvib_im),
                      "up_Hvib_re":data_conv.matrix2list(up_Hvib_re),
                      "up_Hvib_im":data_conv.matrix2list(up_Hvib_im)
                    }

        with open(F"{meta_info_file}.json", 'w') as outfile:
            json.dump(meta_info, outfile)



        
        #============= Output the resulting QSH Hamiltonians ===========================
        Hvib = []
        for i in range(0,nsteps):
            # compute QSH Hvib at time t_i = i * dt
            qs_Hvib = compute_qs_Hvib(nfreqs, freqs, i*dt, H_vib_re_ave, H_vib_re_std, dw_Hvib_re, 
                            up_Hvib_re, H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im,  dev)

            Hvib.append(CMATRIX(qs_Hvib))

            if params["do_QSH_output"]==True:
                #============= Output the resulting QSH Hamiltonians ===========================
                re_filename = params["output_set_paths"][idata] + params["qsh_Hvib_re_prefix"] + str(i) + params["qsh_Hvib_re_suffix"]
                im_filename = params["output_set_paths"][idata] + params["qsh_Hvib_im_prefix"] + str(i) + params["qsh_Hvib_im_suffix"]        
                qs_Hvib.real().show_matrix(re_filename)
                qs_Hvib.imag().show_matrix(im_filename)

        qsh_H_vib.append(Hvib)        
        

    return qsh_H_vib

