import os
import sys
import logging
import numpy as np
import multiprocessing as mp
import glob
import matplotlib.pyplot as plt

import util.libutil as comn
from libra_py import pdos

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler('pdos.log')
file_handler.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def compute_pdos(params):
    """ Computes the pDOS.

    Please note that this function only accommodates pDOS files generated using cp2k.
    If you have run step2 via the nbra workflows, the pDOS files should be stored in
    the all_pdosfiles folder.

    Args:
        params (dict): The parameters dictionary used for containing 
            the variables needed to compute the pDOS 

    Returns:
        None: but produces plots that show the pDOS.
    """

    logger.debug("Entered into the function compute_pdos")

    critical_params = [ ]
    default_params = { "pdos_type":"atom_resolved", "thermal":"False", }
    comn.check_input(params, default_params, critical_params)
    logger.debug("Checked params in the function initialize_step2_jobs")

    pdos_type = params["pdos_type"]
    thermal = params["thermal"]

    if pdos_type == "orbital_resolved":
        logger.debug("Using pdos_type orbital_resolved")
        #===================== Orbital resolved columns
        #                               C, total            C, s               C, p                C, d
        angular_momentum_cols = [ [ list(range(3,12)), list(range(3,4)), list(range(4,7)), list(range(7,12)) ], \
        #                               H, total            H, s               H, p        
                                  [ list(range(3,7)), list(range(3,4)), list(range(4,7)) ] ]
        # This is for orbital resolved
        labels = ['C, s','C, p','C, d','H, s','H, p']
        # Colors, the color orders are based on the labels, Here are chosen for elements
        colors = ['red','blue','green','orange', 'purple']
        outname = 'orbital'

    elif pdos_type == "atom_resolved":
        logger.debug("Using pdos_type atom_resolved")
        #===================== Atom resolved columns
        #                           C, total          C, total
        angular_momentum_cols = [ [ list(range(3,12)), list(range(3,12))],\
        #                           H, total          H, total
                                  [ list(range(3,7)), list(range(3,7)) ]  ]
        # This is for atom resolved
        labels = ['C','H']
        # This row is for atom resolved
        colors = ['blue','red']
        outname = 'atom'

    #===================== Other inputs
    # The time step is MD, For static calculations we only use 0
    time_step = 0
    sigma = 0.1
    coef = 1
    npoints = 2000 # number of points for the grid mesh vector
    energy_conversion = 27.211386 # Hartree to eV
    nprocs = 12 # number of processors

    # create the pool of processors
    pool = mp.Pool(nprocs)

    # This variable will contain the convolved PDOS for each element and each angular momentum list
    Total_results = []
    homos = []
    #========================== Convolution
    logger.debug("Beginning convolution")
    # Here we change the order,
    # In fact the k1 is C, k2 is H
    for k in [1,2]:

        # Create an empty list for summation of the convolved PDOS
        # This is used for Total DOS
        dos_summation_angular = []
        energy_grid_ave = []

        # Define a zero vector to sum the convolved PDOS for each angular momentum list
        # It is based on the number of points
        zero_vec = np.zeros((npoints))
        # Now for each angular momentum append the zero vector
        # The same for enegy grid average since we want the energy grid average
        for i in range(len(angular_momentum_cols[k-1])):
            dos_summation_angular.append(zero_vec)
            energy_grid_ave.append(zero_vec)

        # Create the variables for pool.starmap
        vars_for_pool = []
        # Find all the pdos files of an element in all_pdosfiles for the first trajectory

        if thermal == True:
            DOS_files1 = glob.glob('all_pdosfiles/*k%d-1.pdos'%k)
        else:
            DOS_files1 = glob.glob('all_pdosfiles/*k%d-1.pdos'%k)

        for DOS_file in DOS_files1:
            params2 = {}
            params2["cp2k_pdos_file"] = DOS_file
            params2["time_step"] = time_step
            params2["sigma"] = sigma
            params2["coef"] = coef
            params2["npoints"] = npoints
            params2["energy_conversion"] = energy_conversion
            params2["angular_momentum_cols"] = list(angular_momentum_cols[k-1])
            vars_for_pool.append(params2)

        results_for_k = pool.map(pdos.convolve_cp2k_pdos, vars_for_pool)

        # We initialize all the homos average: homos_ave
        # This variable is the same for each element so it will be repeated but doesn't 
        # change the results
        homos_ave = 0
        # Now we need to take the average of them (The same is for average HOMO eergy level
        # and the energy grid as well)
        # Energy grid is the 0th element, convolved DOS is the 1st element, HOMO energy levels is the 2nd one
        # 1. We add the convolved ones to dos_summation_angular
        # Here we start to take the averages by summing the results.
        for i in range(len(DOS_files1)):
            energy_grid_ave       += results_for_k[i][0] # First element is the energy grid output by the function
            dos_summation_angular += results_for_k[i][1] # Second element is the convolved DOS
            homos_ave             += results_for_k[i][2] # Third element is the HOMO energy level

        # 2. We take the average for that by dividing by len(DOS_files)
        # 3. We append it to Total_results
        Total_results.append(dos_summation_angular/len(DOS_files1))
        # Uncomment only if you use two trajectories
        #Total_results.append(dos_summation_angular/(len(DOS_files1)+len(DOS_files2)))
        energy_grid_ave /= len(DOS_files1)
        # Uncomment only if you use two trajectories
        #energy_grid_ave /= (len(DOS_files1)+len(DOS_files2))
        homos_ave /= len(DOS_files1)
        # Uncomment only if you use two trajectories
        #homos_ave /= (len(DOS_files1)+len(DOS_files2))
    logger.debug("Ending convolution")

    # Close the pool
    pool.close()
    pool.join()
    #================================== End of convolution

    #================================== Total density
    # We first compute the total density of states through the first computed angular momentum 
    # for [3,12] as is shown above - So we only choose the 0 index
    # Make a zero vector of npoints
    total_density = np.zeros((npoints))
    # i here is each element
    for i in range(len(Total_results)):
        # Sum the total by Total_results of zero index angular momentum column for each element
        total_density += Total_results[i][0]

    #============================== Plotting 
    logger.debug("Begin plotting")
    figure = plt.figure(num=None, figsize=(3.21, 2.41), dpi=1200, edgecolor='black', frameon=True)

    # Plot the total density by black color
    plt.plot(energy_grid_ave[0]-homos_ave,total_density,label='Total',color='black', linewidth=2.0)

    # set up a counter for labels
    c = 0
    for i in range(len(Total_results)):
        for j in range(1,len(Total_results[i])):
            plt.plot(energy_grid_ave[0]-homos_ave,Total_results[i][j],label=labels[c],color=colors[c], linewidth=2)
            c += 1
    #plt.ylim(0,60)
    plt.legend(fontsize=6.75, ncol=1, loc='upper center')
    plt.xlabel('Energy, eV',fontsize=12)
    plt.ylabel('DOS, 1/eV',fontsize=12)

    if thermal == True:
        plt.title("Adamantane, 300 K",fontsize=12)
    else:
        plt.title("Adamantane, 0 K",fontsize=12)

    plt.tight_layout()
    plt.show()

    if thermal == True:
        if pdos_type == "orbital_resolved":
            plt.savefig('orbital_DOS_300K_'+outname+'_average.png', dpi=300)
        elif pdos_type == "atom_resolved":
            plt.savefig('atom_DOS_300K_'+outname+'_average.png', dpi=300)

    elif thermal == False:
        if pdos_type == "orbital_resolved":
            plt.savefig('orbital_DOS_0K_'+outname+'.png', dpi=300)
        elif pdos_type == "atom_resolved":
            plt.savefig('atom_DOS_0K_'+outname+'.png', dpi=300)


