#*********************************************************************************                     
#* Copyright (C) 2017 Wei Li, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file calculate_ir.py 
# This script computes the infrared spectrum of alizarin molecule via FT of the 
# autocorrelation function of the time dependent dipole moment information
# the following input files should be presented at the current directory
#
# the equilibrated  MD trajectory  ---> _mol_md.xyz
# the partial charges              ---> charge

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


from libra_py import *



def get_Num_atoms(fname):
# get the number of the atoms from the trajectory file
# fname - the xyz trajectory file 

    with open(fname,'r') as fo:
        Natom = int(fo.next())
		
    return Natom

def read_3d_data(fname,num_atoms):
# read the atom coordinates from the xyz formated file
# fname - the xyz trajectory file
# num_atoms - the numer of atoms

    # here we preset a infinitely larges timestep
    time_step = 0
    coords = [[ [ 0 for i in xrange(3) ] for j in range(num_atoms) ] \
        for k in xrange(1000000)]
	
	
    with open(fname,'r') as fo:
        for line in fo:
            try:
                fo.next()
            except StopIteration:
                break
            for n in range(0,num_atoms):
                line = fo.next()
                coords[time_step][n]= line.split()[1:]
				
            time_step += 1
			
			
    return (coords,time_step)

def read_charge(fname):
# read the partial charge in
# fname - the partial charge file

    charge = []
	
    with open(fname,'r') as fo:
        for line in fo:
            charge.append(float(line))
			
    total_charges = sum(charge)
    print 'total charges is: ', total_charges
    # the total charge should not be strickly zero due to the nurmerical
    # uncertainty
	
    return charge
	
def calculate_dipole(coords,charges,num_atoms,frame):
# calculate the dipole moment for a frame of the trajectory
# we use the calssical approximation to compute the dipole
# Dipole = sum_i {charge * position}
# we assume the partial charge does not change as time moving on
# coords - the list contains the coordinates for a frame
# charges - the list contain the partial charge
# num_atoms - the number of atoms
# frame - the index of which frame to compute


    # calculate the centroid first
    centroid = [0.0] * 3
	
    for i in xrange(num_atoms):
        centroid[0] += float(coords[frame][i][0])
        centroid[1] += float(coords[frame][i][1])
        centroid[2] += float(coords[frame][i][2])
        for j in range(0,3):
            centroid[j] /= num_atoms
	
    # find the relative position of each charge to the centroid,
    # then calculate the dipole
    dipole = VECTOR()
	
    for i in range(0,num_atoms):
        position = [float(coords[frame][i][0]), \
		           float(coords[frame][i][1]), \
		           float(coords[frame][i][2])]
		           
        rel = VECTOR( (position[0] - centroid[0]), \
		              (position[1] - centroid[1]), \
		              (position[2] - centroid[2]) )
		              
        dipole += (float(charges[i]) * rel)
		
		
    return dipole*k	

def calculate_total_dipole(coords,charges,num_atoms,nframes):
# calculate the diple moment of all frames for the given trajectory
# coords - the 3d list contains the coordinates for the trajectory
# charges - the list contain the partial charge
# num_atoms - the number of atoms
# nframes - the number of frames of the trajectory

    total_dipole = []
    dipole_norm_total = []
	
    for i in range(0,nframes):
        dipole = calculate_dipole(coords,charges,num_atoms,i)
        total_dipole.append(dipole)
        dipole_norm_total.append(dipole.length())


    print 'the average dipole (in D):    ',  \
	       sum(dipole_norm_total)/len(dipole_norm_total) 
	
	  
    return total_dipole




if __name__ == '__main__':
  
    # Some physical constants   
    elementary_q = 1.602176487e-19     # elementary charge in C
    debye_conv = 3.33564e-30           # 1 Debye in C*m
    k = elementary_q*1e-10/debye_conv  # conversion factor

    ev2iv_cm = 8065.54468111324

    # input files, can be changed by usrs
    fname_coord="_mol_md.xyz"
    fname_charge="charge"   
    
    
    data = [] # list to store the time-dependent dipole moment

    # Steps to read the coordinate and calculate the dipole
    num_atoms = get_Num_atoms(fname_coord)
    coords,nframes = read_3d_data(fname_coord,num_atoms)
    charges = read_charge(fname_charge)
    data=calculate_total_dipole(coords,charges,num_atoms,nframes)

    # now do the ACF and FT
    acf.recipe1(data, 1.0, 2000.0, 1.0)
    
    print "job was done"
    
    


