#*********************************************************************************                     
#* Copyright (C) 2019-2020 Alexey V. Akimov                                                   
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

       List of classes:

           * class mem_saver
           * class hdf5_saver

       List of functions:

           # TSH
           * init_hdf5(saver, hdf5_output_level, _nsteps, _ntraj, _ndof, _nadi, _ndia)
           * save_hdf5_1D(saver, i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot)
           * save_hdf5_2D(saver, i, states)
           * save_hdf5_3D(saver, i, pops, dm_adi, dm_dia, q, p, Cadi, Cdia)
           * save_hdf5_4D(saver, i, tr, hvib_adi, hvib_dia, St, U, projector)

           # HEOM
           * heom_init_hdf5(saver, hdf5_output_level, _nsteps, _nquant)
           * heom_save_hdf5_1D(saver, i, dt)
           * heom_save_hdf5_3D(saver, i, denmat)

           # Exact
           * exact_init_hdf5(saver, hdf5_output_level, _nsteps, _ndof, _nstates, _ngrid)
           * exact_init_custom_hdf5(saver, _nsteps, _ncustom_pops, _nstates)


.. moduleauthor:: Alexey V. Akimov
  
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


class mem_saver:
    """
    This class is needed for saving variables into a dictionary

    Example of usage:
    
        x = mem_saver(["q", "p"])

        x.add_data("q", 1.0)
        x.add_data("q", -1.0)

        print(x.data["q"])

    """


    def __init__(self, _keywords=[]):
        """
        Initializes an object that would store the data
        in different formats - the temporary and the numpy-consistent
        
        """
        
        # Names of the data sets
        self.keywords = list(_keywords)  
        
        # "Raw" data - elements could be of whatever data type
        self.data = {}  
        
        # "Numpy" data - elements are numpy arrays
        self.np_data = {}
        
        # Only initialize the "raw" data, don't touch the numpy
        for keyword in self.keywords:
            self.data[keyword] = []            

            
    def add_data(self, data_name, _data):
        """
        To collect arbitrary data
        
        Args:
            data_name (string): the name of the data set
            _data (anything): the data to be added
        
        
        This function simply appends the data elements
        to the prepared lists. There is no restriction on the 
        data type for the elements added
        """
        if keyword in self.keywords:
            self.data[data_name].append(_data)

                        
    def add_dataset(self, data_name, shape, dtype):
        """
        To initialize the numpy arrays to hold the data
        
        Args:
            data_name (string): the name of the data set
            shape (tuple of ints): the dimensions of the numpy array
            dtype (one of ["I", "R", or "C"]: the tye of data to be stored in the array
            target (int): 0 - np_data, 1 - np_data_current
            
        """
        if data_name in self.keywords:
            pass
        else:
            self.keywords.append(data_name)


        if dtype=="I":
            self.np_data[data_name] = np.empty(shape, int)
            
        elif dtype=="R":
            self.np_data[data_name] = np.empty(shape, float)            
            
        elif dtype=="C":
            self.np_data[data_name] = np.empty(shape, complex)            

            
        else:
            print(F"ERROR: the dtype = {dtype} is not allowed in add_dataset")
                
    
              
    def save_scalar(self, istep, data_name, _data):
        """
        Saves a scalar to 1D array
        """

        if data_name in self.keywords and data_name in self.np_data.keys():
            self.np_data[data_name][istep] = _data

                        
    def save_multi_scalar(self, istep, iscal, data_name, _data):
        """
        Saves a sacalar to 2D array
        """

        if data_name in self.keywords and data_name in self.np_data.keys():            
            self.np_data[data_name][istep, iscal] = _data

        
    def save_matrix(self, istep, data_name, _data):
        """          
        Saves a matrix to a 3D array
        
        Args:

          istep ( int ) :  index of the timestep for the data
          data_name ( string ) : how to call this data set internally
          _data ( (C)MATRIX(nx, ny) ) : the actual data to save          
           
        """

        if data_name in self.keywords and data_name in self.np_data.keys():        
            nx, ny = _data.num_of_rows, _data.num_of_cols        
            for i in range(nx):
                for j in range(ny):     
                    self.np_data[data_name][istep, i, j] = _data.get(i, j)


    def save_multi_matrix(self, istep, imatrix, data_name, _data):
        """          
        Saves one of a series of matrices

        Args:
        
          istep ( int ) :  index of the timestep for the data
          data_name ( string ) : how to call this data set internally
          data ( (C)MATRIX(nx, ny) ) : the actual data to save
           
        """

        if data_name in self.keywords and data_name in self.np_data.keys():        
            nx, ny = _data.num_of_rows, _data.num_of_cols        
            for i in range(nx):
                for j in range(ny):       
                    self.np_data[data_name][istep, imatrix, i, j] = _data.get(i, j)

    
    def save_data(self, filename, data_names, mode):
        """
        To save the numpy data into HDF5 files
        
        Args:
            filename (string): the name of the HDF5 file where to save the res
            data_names (list of strings): the list of the names of the data sets to save
            mode ("w" or "a"): whether to overwrite the file or to append to it
        
        """

        print("In mem_saver.save_data()")
        print("data_name = ", data_names)        
        print("keywords = ", self.keywords)
        print("keys = ", self.np_data.keys() )

        #for key in self.np_data.keys():
        #    print(key, self.np_data[key] )
  
        
        with h5py.File(filename, mode) as f:
            
            for data_name in data_names:                
                if data_name in self.np_data.keys():

                    print(F"Saving the dataset named {data_name}/data")
                                                        
                    g = f.create_group(data_name)                                                            
                    g.create_dataset("data", data = self.np_data[data_name])
                    
                else:
                    print(F"{data_name} is not in the list {self.np_data.keys()}" )


    
    def save_data_txt(self, prefix, data_names, mode, istep):
        """
        To save the numpy data into TXT files in a given directory
        
        Args:
            prefix (string): the name of the directory, where to save the results
            data_names (list of strings): the list of the names of the data sets to save
            mode ("w" or "a"): whether to overwrite the file or to append to it
        
        """

        if mode=="w":

            print("In mem_saver.save_data_txt()")
            print("data_name = ", data_names)        
            print("keywords = ", self.keywords)
            print("keys = ", self.np_data.keys() )
            
            for data_name in data_names:
                if data_name in self.np_data.keys():
                    f = open(F"{prefix}/{data_name}.txt", "w")
                    f.close()

        elif mode=="a":  # append the istep data to the files

            for data_name in data_names:
                if data_name in self.np_data.keys():

                    #print(F"{data_name}")
                    f = open(F"{prefix}/{data_name}.txt", "a")
                    X = self.np_data[data_name]
                    shp = X.shape  # dimentionality of the data

                    if len(shp)==1:
                        x = X[istep]

                        if X.dtype==complex:
                            f.write(F"{x.real:10.5e} {x.imag:10.5e}\n") 
                        elif X.dtype==float:
                            f.write(F"{x:10.5e}\n") 
                        else:
                            f.write(F"{x}\n") 

                    elif len(shp)==2:
                        for i1 in range(shp[1]):
                            x = X[istep, i1]

                            if X.dtype==complex:
                                f.write(F"{x.real:10.5e} {x.imag:10.5e}") 
                            elif X.dtype==float:
                                f.write(F"{x:10.5e}") 
                            else:
                                f.write(F"{x}") 
                            f.write("  ")
                        f.write("\n")

                    elif len(shp)==3:
                        for i1 in range(shp[1]):
                            for i2 in range(shp[2]):
                                x = X[istep, i1, i2]

                                if X.dtype==complex:
                                    f.write(F"{x.real:10.5e} {x.imag:10.5e}") 
                                elif X.dtype==float:
                                    f.write(F"{x:10.5e}") 
                                else:
                                    f.write(F"{x}") 
                                f.write("  ")
                            f.write("    ")
                        f.write("\n")

                    elif len(shp)==4:
                        for i1 in range(shp[1]):
                            for i2 in range(shp[2]):
                                for i3 in range(shp[3]):
                                    x = X[istep, i1, i2, i3]

                                    if X.dtype==complex:
                                        f.write(F"{x.real:10.5e} {x.imag:10.5e}") 
                                    elif X.dtype==float:
                                        f.write(F"{x:10.5e}") 
                                    else:
                                        f.write(F"{x}") 
                                    f.write("  ")
                                f.write("    ")
                            f.write("        ")
                        f.write("\n")

                    f.close()
                    
                else:
                    #print(F"{data_name} is not in the list {self.np_data.keys()}" )
                    pass





class hdf5_saver:

    def __init__(self, _filename, _keywords=[]):
        """
        The constructor of the class objects

        Args:
            _filename ( string ): the name of the HDF5 file to be created
            _keywords ( list of strings ): the names of the data to be added, if the 
                the `data_name` argument in any of the `save_*` functions doesn't exist 
                in the provided list of keywords, the data will not be actually saved into
                the HDF5 file

        Example:
            saver = hdf5_saver("data.hdf")

        """

        self.keywords = list(_keywords)

        self.filename = _filename
        self.use_compression = 1
        self.c_compression_level = 4
        self.r_compression_level = 4
        self.i_compression_level = 9

        with h5py.File(self.filename, "w") as f:
            g = f.create_group("default")
            g.create_dataset("data", data=[])

        print("HDF5 saver is initialized...")
        print(F"the datasets that can be saved are: {self.keywords}")


    def add_keywords(self, _keywords):

        for keyword in _keywords:
 
            if keyword not in self.keywords:
                self.keywords.append(keyword)



    def set_compression_level(self, _use_compression, _compression_level):
        """
        To control how much of data compression we want to exercise

        """

        self.use_compression = _use_compression;
        self.c_compression_level = _compression_level[0]
        self.r_compression_level = _compression_level[1]
        self.i_compression_level = _compression_level[2]



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

            if self.use_compression==1:

                if data_type == "C":            
                    g.create_dataset("data", dim, dtype=complex, maxshape=dim, compression="gzip", compression_opts = self.c_compression_level)
                elif data_type == "R":            
                    g.create_dataset("data", dim, dtype=float, maxshape=dim, compression="gzip", compression_opts = self.r_compression_level)
                elif data_type == "I":            
                    g.create_dataset("data", dim, dtype=int, maxshape=dim, compression="gzip", compression_opts = self.i_compression_level)

            else:
                if data_type == "C":            
                    g.create_dataset("data", dim, dtype=complex, maxshape=dim)
                elif data_type == "R":            
                    g.create_dataset("data", dim, dtype=float, maxshape=dim)
                elif data_type == "I":            
                    g.create_dataset("data", dim, dtype=int, maxshape=dim)


    def save_scalar(self, istep, data_name, data):

        if data_name in self.keywords:
            with h5py.File(self.filename, "a") as f:
                f[F"{data_name}/data"][istep] = data


    def save_multi_scalar(self, istep, iscal, data_name, data):

        if data_name in self.keywords:
            with h5py.File(self.filename, "a") as f:
                f[F"{data_name}/data"][istep, iscal] = data


    
    def save_matrix(self, istep, data_name, data):
        """          
          Add a matrix

          istep ( int ) :  index of the timestep for the data
          data_name ( string ) : how to call this data set internally
          data ( (C)MATRIX(nx, ny) ) : the actual data to save
           
        """

        if data_name in self.keywords:

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
           
        """

        if data_name in self.keywords:

            nx, ny = data.num_of_rows, data.num_of_cols

            with h5py.File(self.filename, "a") as f:

                for i in range(nx):
                    for j in range(ny):
                        f[F"{data_name}/data"][istep, imatrix, i, j] = data.get(i, j)






