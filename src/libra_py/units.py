#*********************************************************************************
#* Copyright (C) 2018-2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: units
   :platform: Unix, Windows
   :synopsis: This module defines the units conversion factors and implements the 
       functions for unit conversion

.. moduleauthor:: Alexey V. Akimov

"""


# Unit conversions Parameters
# Length units
Angst = 1.889725989                                # in Bohr

# Energy units
inv_cm2ev = 1.239841870000000847635400992705E-4    # (1.0/8065.54468111324)
ev2Ha = 0.0367498438131637940538752710301          # (1.0/27.211)    # 27.2 ev is 1 Ha 
ev2au = 0.0367498438131637940538752710301          # (1.0/27.211)    # 27.2 ev is 1 Ha 
au2ev = 27.211                                     # 1 Ha = 27.211 eV
inv_cm2Ha = 4.5563995075520960186520193771085E-6   # inv_cm2ev * ev2Ha
hartree = 627.5094709                              # 1 Ha = 627.5.. kcal/mol

# Time units
fs2au = 41.339396444811905746176105828855          # (1.0/0.02419)   # 41.34 a.u. is 1 fs 
au2fs = 0.02419                                    # 41.34 a.u. is 1 fs
ps2au = 41.339396444811905746176105828855E+3       # 
au2ps = 0.02419E-3                                 # 

au2wavn = 219471.53631777237364                    # 27.211 * 8065.54468111324
wavn2au = 4.5563995075520960186520193771085E-6     # 1/au2wavn

# Mass units
amu = 1836.0                                       # chemical mass unit in masses of electron

# Fundamental constants
boltzmann = 1.9872065E-3                           # in kcal/mol*K
kB = 3.166815151251607985889922605788E-6           # Ha/mol*K



def length_converter(inp_units, out_units):
    """Length units conversion factor

    Args:
        inp_units ( int ): defines the input units
        out_units ( int ): defines the output units

    Note:
        Both inp_units and out_units are encoded as follows:

            * 0 - atomic units (Bohrs)
            * 1 - Angstroms

    Returns:
        double: the length units conversion factor
         
    """

    scl = 1.0
    if inp_units==0:  # input in Bohrs
        if out_units==0:  # output in Bohrs
            scl = 1.0
        elif out_units==1:  # output in Angstroms
            scl = 1.0/Angst

    elif inp_units==1:  # input in Angstroms
        if out_units==0:  # output in Bohrs
            scl = Angst
        elif out_units==1:  # output in Angstroms
            scl = 1.0

    return scl
