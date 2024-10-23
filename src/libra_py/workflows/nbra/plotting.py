#*********************************************************************************
#* Copyright (C) 2024 Kosar Yasin, Mohammad Shakiba, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#
#
"""

.. module:: plotting
   :platform: Unix, Windows
   :synopsis: This module implements function for plotting the results of NBRA calculations
       Example:
.. moduleauthor:: Kosar Yasin, Mohammad Shakiba, Alexey V. Akimov

"""

import os, warnings

import matplotlib.pyplot as plt   # plots
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit
import h5py
from scipy.constants import physical_constants


#from liblibra_core import *
#import libra_py.common_utils as comn
import util.libutil as comn
#from libra_py import units
#from libra_py.workflows.nbra import lz, step4
#from libra_py import data_outs, data_conv


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

    critical_params = [ ]
    default_params = { "title_1":"System", "filename":"excess_energy.png", "do_show":0 }
    comn.check_input(params, default_params, critical_params)


    title_1 = params["title_1"]  # the title of the figure
    filename = params["filename"]
    do_show = params["do_show"]


    fig, ax = plt.subplots(figsize=(16, 12))  # Larger figure size for the combined plot

    ax.set_title(F"{title_1}", fontsize=60)
    ax.set_xlabel('Time, fs', fontsize=54)
    ax.set_ylabel('Excess Energy, eV', fontsize=54)
    ax.tick_params(axis='both', which='major', labelsize=40)
    ax.legend(fontsize=24, loc='upper right', ncol=3, frameon=False  )
    plt.tight_layout()
    plt.savefig(F"{filename}", dpi=600, bbox_inches='tight')
    if do_show==True:
        plt.show()



