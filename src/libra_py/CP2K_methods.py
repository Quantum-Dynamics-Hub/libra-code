#*********************************************************************************
#*
#* Copyright (C) 2020 Mohammad Shakiba, Brendan Smith, Alexey V. Akimov
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: CP2K_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing the CP2K inputs and outputs.
.. moduleauthors:: 
       Mohammad Shakiba, Brendan Smith, Alexey V. Akimov 
  
"""


import os
import sys
import math
import re
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *




def ndigits( integer_number: int ):
    """
    This function calculates the number of digits for an integer number.

    Args:

        integer_number ( integer ): An integer number.

    Returns:
 
        digit_num ( integer ): The number of digits for 'num'.
    """

    # Setting up the digit number counter
    digit_num = 0

    while( integer_number > 0 ):

        digit_num = digit_num+1
        integer_number = integer_number//10
        
    return digit_num



def state_num_cp2k(state_num: int):
    """
    This function turns an integer number to a string. This string 
    has the format in five digits and if the number of digits is less than
    five it will append zeros before the number. This is necessary to read
    the .cube files produced by CP2K.

    Args:

        state_num (integer): Integer number which is in fact the energy level.

    Returns:

        state_num_str (string): The number in five digits (or more) in string format.

    """
    # Setting up an initial string
    state_num_str = ''

    # Now count the number of digits for the state sumber
    # and then add 0 before the number in the string

    if ndigits( state_num ) == 1:

        state_num_str = '0000%d' % state_num

    elif ndigits( state_num ) == 2:

        state_num_str = '000%d' % state_num

    elif ndigits( state_num ) == 3:

        state_num_str = '00%d' % state_num

    elif ndigits( state_num ) == 4:

        state_num_str = '0%d' % state_num

    else:
        state_num_str = '%d' % state_num
    
    return state_num_str




def cube_file_names_cp2k( static_or_dynamics, project_name, min_band, max_band, time_step, spin ):
    """
    This function will produce the cube file names produced by CP2K both for energy calculations
    and molecular dynamics. 

    Args:

        static_or_dynamics ( integer ): This variable has two options:
                                      0: static calculations
                                      1: molecular dynamics calculations

        project_name ( string ): The project name used in CP2K input file.

        min_band ( integer ): The minimum state number to be considered.

        max_band ( integer ): The maximum state number to be considered.

        time_step ( integer ): The time step in CP2K calculations.

        spin ( integer ): This variable can get only two options:
                          1 for alpha spin
                          2 for beta spin

    Returns:

        names ( list ): The list of cube file names.
	
    """
    # Define the names of the cube files for time t.
    cube_file_names = []

    # Dynamics calculations
    if static_or_dynamics == 1:

        for i in range( min_band, max_band + 1 ):

            # Using the state_num_cp2k function to obtain the names of the .cube files for state 'i'.
            state_num = state_num_cp2k( i )

            # For time step '( t )'
            cube_file_name_of_band_i = str( '%s-WFN_%s_%d-1_%d.cube' % ( project_name, state_num, spin, time_step ) )
            cube_file_names.append( cube_file_name_of_band_i )

    # Static calculations
    elif static_or_dynamics == 0:

        for i in range( min_band, max_band+1 ):

            # Using the int_to_five_digit_str to obtain the names of the .cube files for state 'i'.
            state_num = state_num_cp2k( i )

            # For time step '( t )'
            cube_file_name_of_band_i = str('%s-%d-WFN_%s_%d-1_0.cube' % ( project_name, time_step, state_num, spin ) )
            cube_file_names.append( cube_file_name_of_band_i )
			
    return cube_file_names



def read_cp2k_tddfpt_log_file( tddfpt_log_file_name, num_ci_states, tolerance ):
    """
    This function reads log files generated from TD-DFPT calculations using CP2K and returns the TD-DFPT
    excitation energies, Slater determinent states in terms of the Kohn-Sham orbital indicies, 
    and the Slater determinent coefficients that comprise the multi-configurational electronic states
    
    Args:

        tddfpt_log_file_name ( string ): name of the log file for this particular timestep

        num_ci_states ( int ): how many ci states to consider

        tolerance ( float ): the tolerance for accepting SDs are members of the CI wavefunctions

    Returns:

       excitation_energies ( list ): The excitation energies in the log file.

       ci_basis ( list ): The CI-basis configuration.

       ci_coefficients ( list ): The coefficients of the CI-states.

    """
    
    f = open(tddfpt_log_file_name,'r')
    lines = f.readlines()
    f.close()

    for i in range( 0, len(lines) ):

        if 'excitation' in lines[i].lower().split():
            if 'analysis' in lines[i].lower().split():

                # When found the line in which contains 'Excitation analysis' in 
                # the log file, append it in the variable exc_anal_line
                exc_anal_line = i
        
        if 'states' in lines[i].lower().split():
            if 'multiplicity' in lines[i].lower().split():

                # Here we search for the line that contains
                # 'R-TDDFPT states of multiplicity 1' in the log file
                r_tddfpt_line = i

    excitation_energies = []
    
    # Start from 5 lines after finding the line contaning 'R-TDDFPT states of multiplicity 1'
    # This is because they contain the energies from that line.
    for i in range( r_tddfpt_line+5, len( lines ) ):
        if len( lines[i].split() ) == 0:
            break
        excitation_energies.append( float( lines[i].split()[2] ) )
        
    # Start from 5 lines after finding the line contaning 'Excitation analysis'
    # From that point we have the state numbers with their configurations.
    # So, we append the lines which contain only 'State number' and stop
    # whenever we reach to a blank line.
    state_num_lines = []
    for i in range( exc_anal_line+5, len( lines ) ):

        if len( lines[i].split() ) == 0:
            state_num_lines.append( i )
            break

        if len( lines[i].split() ) == 1:
            state_num_lines.append( i )

    # Setting up the CI-basis list and their coefficients list
    ci_basis        = []
    ci_coefficients = []
    for i in range( num_ci_states ):
        
        # CI states and their coefficients for each excited state
        tmp_ci_state              = []
        tmp_ci_state_coefficients = []
        for j in range( state_num_lines[i]+1, state_num_lines[i+1] ):

            if abs( float( lines[j].split()[2] ) ) > tolerance:
                tmp_ci_state.append( [ int( lines[j].split()[0] ), int( lines[j].split()[1] ) ]  )
                tmp_ci_state_coefficients.append( float( lines[j].split()[2] ) )

        # Append the CI-basis and and their coefficients for
        # this state into the ci_basis and ci_coefficients lists
        ci_basis.append( tmp_ci_state )
        ci_coefficients.append( tmp_ci_state_coefficients )
        
    return excitation_energies, ci_basis, ci_coefficients

