#*********************************************************************************
#* Copyright (C) 2017-2019 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
"""
.. module:: common_untils
   :platform: Unix, Windows
   :synopsis: This module implements various auxiliary functions data handling

.. moduleauthor:: Alexey V. Akimov

"""

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import units



def check_input(params, default_params, critical_params):
    """

    This function checks the input Python dictionary and verifies if it contains
    all the required keys. Exits if there are no. It also checks other variables 
    and sets the values of these keys to the default values, if they weren't defined

    Args:
        params ( Python dictionary ): the dictionary (usually containing parameters) to be checked
            This dictionary is passed by reference and can be changed inside this function.

        default_params ( Python dictionary ): the dictionary with the default key-value pairs to 
            be inserted into the `params` variable, if those pairs are not defined in it.

        critical_params ( list of strings ): the names of the keys that must be present in the
            `params` variable. The function exits and prints out an error, if any of the keys in 
            the `critical_params` variable are not found in `params`.

    Returns:
        None: 

            This function checks whether the required keys are defined in the input.
            If not - print out a worning and set the parameters to the default values provided.

            But the critical parameters should be present - otherwise will exit

            Revises the parameters dictionary - but only according to the subset of 
            default parameters - other ones are not affects

    Example:

        Assume there is a function called `my_function`, which prints a greeting message.
        So, it needs to know the name. There is too much of an ambiguity to set up a default name.
        Therefore, the "Name" key is required (critical) for the function not to crash. 
        Also, the function prints out the age of the person. If that info is provided, no problem.
        But if the user forgets about it, one can guess something more-or-less reasonable. So,
        in this example, the "Age" key is set up to a default value before we call the `my_function`.
          
        >>> def my_function(data, params):
        >>>     print "Hello ", params["Name"]
        >>>     print "You age is ", params["Age"]


        The following snippet 

        >>> params = {"Name":"Alexey", "Age":33 }
        >>> default_params = { "Department": "Chemistry", "Institution":"UB", "Age":1 }
        >>> critical_parames = ["Name"]
        >>> check_input(params, default_params, critical_params)
        >>> my_function(data, params)

        will produce:
        
        >>> Hello Alexey
        >>> Your age is 33


        The following snippet 

        >>> params = {"Name":"Alexey" }
        >>> default_params = { "Department": "Chemistry", "Institution":"UB", "Age":1 }
        >>> critical_parames = ["Name"]
        >>> check_input(params, default_params, critical_params)
        >>> my_function(data, params)

        will produce:
        
        >>> Hello Alexey
        >>> Your age is 1


        The following snippet 

        >>> params = { }
        >>> default_params = { "Department": "Chemistry", "Institution":"UB", "Age":1 }
        >>> critical_parames = ["Name"]
        >>> check_input(params, default_params, critical_params)
        >>> my_function(data, params)

        will crash:

        >>> ERROR: The critical parameter Name must be defined!
        >>> Exiting now...


    """

    for key in critical_params:
        if key not in params:
            print "ERROR: The critical parameter ", key, " must be defined!"
            print "Exiting now..."
            sys.exit(0)

    for key in default_params:
        if key not in params:
            print "WARNING: Parameter ", key, " is not defined! in the input parameters"
            print "Use the default value = ", default_params[key]
            params.update({key: default_params[key]})
        else:
            pass
            


def orbs2spinorbs(s):
    """
    This function converts the matrices in the orbital basis (e.g. old PYXAID style)
    to the spin-orbital basis.
    Essentially, it makes a block matrix of a double dimension: 
           ( s  0 )
    s -->  ( 0  s )

    This is meant to be used for backward compatibility with PYXIAD-generated data
    """

    sz = s.num_of_cols
    zero = CMATRIX(sz, sz)    
    act_sp1 = range(0, sz)
    act_sp2 = range(sz, 2*sz)
    
    S = CMATRIX(2*sz, 2*sz)

    push_submatrix(S, s, act_sp1, act_sp1)
    push_submatrix(S, zero, act_sp1, act_sp2)
    push_submatrix(S, zero, act_sp2, act_sp1)
    push_submatrix(S, s, act_sp2, act_sp2)

    return S






def printout(t, pops, Hvib, outfile):
    """
    t - time [a.u.] 
    pops - [MATRIX] - populations
    Hvib - [CMATRIX] - vibronic Hamiltonian
    outfile - filename where we'll print everything out
    """

    nstates = Hvib.num_of_cols


    line = "%8.5f " % (t)
    P, E = 0.0, 0.0
    for state in xrange(nstates):
        p, e = pops.get(state,0), Hvib.get(state, state).real
        P += p
        E += p*e
        line = line + " %8.5f  %8.5f " % (p, e)
    line = line + " %8.5f  %8.5f \n" % (P, E)

    f = open(outfile, "a") 
    f.write(line)
    f.close()


    

def add_printout(i, pop, filename):
    # pop - CMATRIX(nstates, 1)

    f = open(filename,"a")
    line = "step= %4i " % i    

    tot_pop = 0.0
    for st in xrange(pop.num_of_cols):
        pop_o = pop.get(st,st).real
        tot_pop = tot_pop + pop_o
        line = line + " P(%4i)= %8.5f " % (st, pop_o)
    line = line + " Total= %8.5f \n" % (tot_pop)
    f.write(line)
    f.close()


