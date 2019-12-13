#*********************************************************************************                     
#* Copyright (C) 2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: dynamics_hdf5
   :platform: Unix, Windows
   :synopsis: This module implements the read/write functions specifically designed to work with the dynamics module
   Uniquely to this module, all the read/write operations involve HDF5 files
.. moduleauthor:: Alexey V. Akimov

  List of functions:

  
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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
from . import units
from . import data_read



class hdf5_saver:

    def __init__(self, _filename):
        """
        The constructor of the class objects

        Args:
            _filename ( string ): the name of the HDF5 file to be created

        Example:
            saver = hdf5_saver("data.hdf")

        """

        self.filename = _filename

        with h5py.File(self.filename, "w") as f:
            g = f.create_group("default")
            g.create_dataset("data", data=[])


    def add_dataset(self, data_set_name, dim, data_type):
        """

        Args:

            data_set_name ( string ): the name of the data set
            dim  ( tuple ) : dimensions of the data set
            data_type ["C", "R", "I"] : for complex, real, or integer
        
        Example:

            _nsteps, _ntraj, _ndof, _nstates = 10, 1, 2, 2

            saver.add_dataset("time", (_nsteps,) , "R")  # for saving the time axis
            saver.add_dataset("total_energies", (_nsteps, _ntraj), "R") # for saving trajectory-resolved energies
            saver.add_dataset("q", (_nsteps, _ntraj, _ndof), "R")  # for saving trajectory-resolved coordinates
            saver.add_dataset("Cadi", (_nsteps, _nstates, _nstates), "C") # for saving trajectory-resolve TD-SE amplitudes
            saver.add_dataset("Hadi", (_nsteps, _ntraj, _nstates, _nstates), "C")  # for saving trajectory-resolved H matrices


        """

        with h5py.File(self.filename, "a") as f:
            g = f.create_group(data_set_name)
            g.attrs["dim"] = dim
            g.attrs["data_type"] = data_type

            if data_type == "C":            
                g.create_dataset("data", dim, dtype=complex, maxshape=dim, compression="gzip", compression_opts = 4)
            elif data_type == "R":            
                g.create_dataset("data", dim, dtype=float, maxshape=dim, compression="gzip", compression_opts = 4)
            elif data_type == "I":            
                g.create_dataset("data", dim, dtype=int, maxshape=dim, compression="gzip", compression_opts = 9)



    def save_scalar(self, istep, data_name, data):

        with h5py.File(self.filename, "a") as f:
            f[F"{data_name}/data"][istep] = data


    def save_multi_scalar(self, istep, iscal, data_name, data):

        with h5py.File(self.filename, "a") as f:
            f[F"{data_name}/data"][istep, iscal] = data



    
    def save_matrix(self, istep, data_name, data):
        """          
          Add a matrix

          istep ( int ) :  index of the timestep for the data
          data_name ( string ) : how to call this data set internally
          data ( (C)MATRIX(nx, ny) ) : the actual data to save
          mode ( string ): how to operate on the data file
 
              "w" : use when calling this function for the first time - the data set will be initialized
              "a" : append - just update the existing dataset, use when calling the function for the following times

          data_type ( int ):
 
               0 : real-valued
               1 : complex-valued
           
        """

        nx, ny = data.num_of_rows, data.num_of_cols

        with h5py.File(self.filename, "a") as f:

            for i in range(nx):
                for j in range(ny):
                    f[F"{data_name}/data"][istep, i, j] = data.get(i, j)



    def save_multi_matrix(self, istep, imatrix, data_name, data):
        """          
          Add a matrix

          istep ( int ) :  index of the timestep for the data
          data_name ( string ) : how to call this data set internally
          data ( (C)MATRIX(nx, ny) ) : the actual data to save
          mode ( string ): how to operate on the data file
 
              "w" : use when calling this function for the first time - the data set will be initialized
              "a" : append - just update the existing dataset, use when calling the function for the following times

          data_type ( int ):
 
               0 : real-valued
               1 : complex-valued
           
        """

        nx, ny = data.num_of_rows, data.num_of_cols

        with h5py.File(self.filename, "a") as f:

            for i in range(nx):
                for j in range(ny):
                    f[F"{data_name}/data"][istep, imatrix, i, j] = data.get(i, j)




def init_hdf5(saver, hdf5_output_level, _nsteps, _ntraj, _ndof, _nadi, _ndia):


    if hdf5_output_level>=1:

        # Time axis (integer steps)
        saver.add_dataset("timestep", (_nsteps,) , "I")  

        # Time axis
        saver.add_dataset("time", (_nsteps,) , "R")  
        
        # Average kinetic energy
        saver.add_dataset("Ekin_ave", (_nsteps,) , "R")  
        
        # Average potential energy
        saver.add_dataset("Epot_ave", (_nsteps,) , "R")  
        
        # Average total energy
        saver.add_dataset("Etot_ave", (_nsteps,) , "R")  
        
        # Fluctuation of average kinetic energy
        saver.add_dataset("dEkin_ave", (_nsteps,) , "R")  
        
        # Fluctuation of average potential energy
        saver.add_dataset("dEpot_ave", (_nsteps,) , "R")  
        
        # Fluctuation of average total energy
        saver.add_dataset("dEtot_ave", (_nsteps,) , "R")  


    if hdf5_output_level>=2:

        # Trajectory-resolved instantaneous adiabatic states
        saver.add_dataset("states", (_nsteps, _ntraj), "I") 


    if hdf5_output_level>=3:

        # Average adiabatic SH populations
        saver.add_dataset("SH_pop", (_nsteps, _nadi, 1), "R") 

        # Average adiabatic density matrices
        saver.add_dataset("D_adi", (_nsteps, _nadi, _nadi), "C") 

        # Average diabatic density matrices
        saver.add_dataset("D_dia", (_nsteps, _ndia, _ndia), "C") 



        # Trajectory-resolved coordinates
        saver.add_dataset("q", (_nsteps, _ntraj, _ndof), "R") 

        # Trajectory-resolved momenta
        saver.add_dataset("p", (_nsteps, _ntraj, _ndof), "R") 

        # Trajectory-resolved adiabatic TD-SE amplitudes
        saver.add_dataset("Cadi", (_nsteps, _ntraj, _nadi), "C") 

        # Trajectory-resolved diabatic TD-SE amplitudes
        saver.add_dataset("Cdia", (_nsteps, _ntraj, _ndia), "C") 


    if hdf5_output_level>=4:

        # Trajectory-resolved vibronic Hamiltoninans in the adiabatic representation
        saver.add_dataset("hvib_adi", (_nsteps, _ntraj, _nadi, _nadi), "C") 

        # Trajectory-resolved vibronic Hamiltoninans in the diabatic representation
        saver.add_dataset("hvib_dia", (_nsteps, _ntraj, _ndia, _ndia), "C") 

        # Trajectory-resolved time-overlaps of the adiabatic states
        saver.add_dataset("St", (_nsteps, _ntraj, _nadi, _nadi), "C") 

        # Trajectory-resolved diabatic-to-adiabatic transformation matrices 
        saver.add_dataset("basis_transform", (_nsteps, _ntraj, _ndia, _nadi), "C") 

        # Trajectory-resolved projector matrices (from the raw adiabatic to consistent adiabatic)
        saver.add_dataset("projector", (_nsteps, _ntraj, _nadi, _nadi), "C") 



def save_hdf5_1D(saver, i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot):
        # Timestep 
        saver.save_scalar(i, "timestep", i) 

        # Actual time
        saver.save_scalar(i, "time", dt*i)  

        # Average kinetic energy
        saver.save_scalar(i, "Ekin_ave", Ekin)  

        # Average potential energy
        saver.save_scalar(i, "Epot_ave", Epot)  

        # Average total energy
        saver.save_scalar(i, "Etot_ave", Etot)  

        # Fluctuation of average kinetic energy
        saver.save_scalar(i, "dEkin_ave", dEkin)  

        # Fluctuation of average potential energy
        saver.save_scalar(i, "dEpot_ave", dEpot)  

        # Fluctuation average total energy
        saver.save_scalar(i, "dEtot_ave", dEtot)  


def save_hdf5_2D(saver, i, states):

        # Trajectory-resolved instantaneous adiabatic states
        # Format: saver.add_dataset("states", (_nsteps, _ntraj), "I")        
        ntraj = len(states)
        for itraj in range(ntraj):
            saver.save_multi_scalar(i, itraj, "states", states[itraj])




def save_hdf5_3D(saver, i, pops, dm_adi, dm_dia, q, p, Cadi, Cdia):
        # Average adiabatic SH populations
        # Format: saver.add_dataset("SH_pop", (_nsteps, _nadi, 1), "R") 
        saver.save_matrix(i, "SH_pop", pops) 


        # Average adiabatic density matrices
        # Format: saver.add_dataset("D_adi", (_nsteps, _nadi, _nadi), "C") 
        saver.save_matrix(i, "D_adi", dm_adi) 

        # Average diabatic density matrices
        # Format: saver.add_dataset("D_dia", (_nsteps, _ndia, _ndia), "C") 
        saver.save_matrix(i, "D_dia", dm_dia) 


        # Trajectory-resolved coordinates
        # Format: saver.add_dataset("q", (_nsteps, _ntraj, _dof), "R") 
        saver.save_matrix(i, "q", q.T()) 

        # Trajectory-resolved momenta
        # Format: saver.add_dataset("p", (_nsteps, _ntraj, _dof), "R") 
        saver.save_matrix(i, "p", p.T()) 

        # Trajectory-resolved adiabatic TD-SE amplitudes
        # Format: saver.add_dataset("C_adi", (_nsteps, _ntraj, _nadi), "C") 
        saver.save_matrix(i, "Cadi", Cadi.T()) 

        # Trajectory-resolved diabatic TD-SE amplitudes
        # Format: saver.add_dataset("C_dia", (_nsteps, _ntraj, _ndia), "C") 
        saver.save_matrix(i, "Cdia", Cdia.T()) 


def save_hdf5_4D(saver, i, tr, hvib_adi, hvib_dia, St, U, projector):

        # Trajectory-resolved vibronic Hamiltoninans in the adiabatic representation
        # Format: saver.add_dataset("hvib_adi", (_nsteps, _ntraj, _nadi, _nadi), "C") 
        saver.save_multi_matrix(i, tr, "hvib_adi", hvib_adi) 


        # Trajectory-resolved vibronic Hamiltoninans in the diabatic representation
        # Format: saver.add_dataset("hvib_dia", (_nsteps, _ntraj, _ndia, _ndia), "C") 
        saver.save_multi_matrix(i, tr, "hvib_dia", hvib_dia) 

        # Trajectory-resolved time-overlaps of the adiabatic states
        # Format: saver.add_dataset("St", (_nsteps, _ntraj, _nadi, _nadi), "C") 
        saver.save_multi_matrix(i, tr, "St", St) 

        # Trajectory-resolved diabatic-to-adiabatic transformation matrices 
        # Format: saver.add_dataset("basis_transform", (_nsteps, _ntraj, _ndia, _nadi), "C") 
        saver.save_multi_matrix(i, tr, "basis_transform", U) 

        # Trajectory-resolved projector matrices (from the raw adiabatic to consistent adiabatic)
        # Format: saver.add_dataset("projector", (_nsteps, _ntraj, _nadi, _nadi), "C") 
        saver.save_multi_matrix(i, tr, "projector", projector) 


