#*********************************************************************************                     
#* Copyright (C) 2020 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: dynamics_exact
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing exact on-the-grid dynamics and for their visualization
.. moduleauthor:: Alexey V. Akimov

  List of functions:

  def init_wfc(params, _potential, model_params ):
  def save_data_hdf5(step, wfc, saver, params):  
  def save_data_mem(step, wfc, saver, params):
  def run_dynamics(wfc, params, model_params, savers):
  def run_relaxation(_params, _potential, model_params):
  def plot_mem(res, _params, model_params, plot_params):
  def plot_hdf5(plot_params):

"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"


import os
import sys
import math
import copy
import time

import h5py
import numpy as np
import matplotlib.pyplot as plt   # plots


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
from . import units
from . import data_outs
from . import data_conv
from . import dynamics_io
from . import dynamics_hdf5





def init_wfc(params, _potential, model_params ):
    """
    This is a generic and a bit excessive initialization procedure - some of the parameters
    may not be needed for some of the methods, but it allows us to have a simple structure of the program
    
    istate = [rep, index_of_the_state]
    
    """
            
    critical_params = []
    default_params = { "nsteps":200, "dt":10.0, "progress_frequency":0.1,
                       "rmin":[-15.0], "rmax":[15.0], "dx":[0.1], "nstates":2,
                       "x0":[0.0], "p0":[0.0], "istate":[1,0], "masses":[2000.0], "k":[0.001]
                      }
    comn.check_input(params, default_params, critical_params)
    
    # Grid properties 
    dx = Py2Cpp_double(params["dx"])
    rmin = Py2Cpp_double(params["rmin"])
    rmax = Py2Cpp_double(params["rmax"])
    nstates = params["nstates"]  

    # Dynamical properties 
    nsteps = params["nsteps"]
    dt = params["dt"]    
    
    # Properties of the initial wavefunction
    istate = params["istate"]          
    k = params["k"]
    masses = Py2Cpp_double(params["masses"])
    

    x0 = Py2Cpp_double(params["x0"])
    p0 = Py2Cpp_double(params["p0"])
    ndof = len(x0)
    
    dx0_tmp = []
    for idof in range(ndof):
        sigmax2 = 0.5*math.sqrt(1.0/(k[idof]*masses[idof]))        
        dx0_tmp.append(math.sqrt(sigmax2))
    dx0 = Py2Cpp_double(dx0_tmp)    
        
        
    
    # Here we initialize the grid and wavefunction
    wfc = Wfcgrid2( rmin, rmax,  dx, nstates)
    
    wfc.direct_allocate_tmp_vars(0)    # last time-step wavefunction, diabatic rep
    wfc.direct_allocate_tmp_vars(1)    # last time-step wavefunction, adiabatic rep

    
    wfc.update_Hamiltonian(_potential, model_params, 0)  # update Hamiltonian: diabatic -
                                                        # need to compute diabatic-to-adiabatic transform
    wfc.update_Hamiltonian(_potential, model_params, 1)  # update Hamiltonian: adiabatic, NACs
        
    wfc.update_propagator_H(0.5*dt)                     # copute the dia-to-adi transform + exp(-i* V *dt/2)
    
    wfc.add_wfc_Gau(x0, p0, dx0, istate[1], 1.0+0.0j, istate[0])   # Add to the diabatic or adiabatic state
        
    if istate[0]==0:
        # If we added a diabatic state:
        wfc.update_adiabatic()        # update adiabatic wavefunction    
            
    elif istate[0]==1:
        # If we added an adiabatic state:
        wfc.update_diabatic()        # update diabatic wavefunction            
    
    
    wfc.update_propagator_K(dt, masses) # update reci space propagator in diabatic rep, exp(-iTdt) 
    wfc.update_reciprocal(1)     # update reci of adiabatic function        
    wfc.update_reciprocal(0)     # update reci of diabatic function
        
        
    print( "Norm (dia) = ", wfc.norm(0) )
    print( "Norm (adi) = ", wfc.norm(1) )
    print( "Ekin (dia) = ", wfc.e_kin(masses, 0) )
    print( "Ekin (adi) = ", wfc.e_kin(masses, 1) )
    print( "Epot (dia) = ", wfc.e_pot(0) )
    print( "Epot (adi) = ", wfc.e_pot(1) )

    
    return wfc



def save_data_hdf5(step, wfc, saver, params):
    
    dt = params["dt"]    
    masses = Py2Cpp_double(params["masses"])
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level>=1:
        saver.save_scalar(step, "timestep", step) 
        saver.save_scalar(step, "time", step * dt)        
        saver.save_scalar(step, "Ekin_dia", wfc.e_kin(masses, 0)) 
        saver.save_scalar(step, "Ekin_adi", wfc.e_kin(masses, 1))         
        saver.save_scalar(step, "Epot_dia", wfc.e_pot(0)) 
        saver.save_scalar(step, "Epot_adi", wfc.e_pot(1))         
        saver.save_scalar(step, "Etot_dia", wfc.e_tot(masses, 0)) 
        saver.save_scalar(step, "Etot_adi", wfc.e_tot(masses, 1))                 
        saver.save_scalar(step, "norm_dia", wfc.norm(0) )
        saver.save_scalar(step, "norm_adi", wfc.norm(1) )
        
        
    if hdf5_output_level>=2:    

        saver.save_matrix(step, "pop_dia", wfc.get_pops(0) ) 
        saver.save_matrix(step, "pop_adi", wfc.get_pops(1) )     
        saver.save_matrix(step, "q_dia", wfc.get_pow_q(0, 1) ) 
        saver.save_matrix(step, "q_adi", wfc.get_pow_q(1, 1) )         
        saver.save_matrix(step, "q2_dia", wfc.get_pow_q(0, 2) ) 
        saver.save_matrix(step, "q2_adi", wfc.get_pow_q(1, 2) )         
        saver.save_matrix(step, "p_dia", wfc.get_pow_p(0, 1) ) 
        saver.save_matrix(step, "p_adi", wfc.get_pow_p(1, 1) )         
        saver.save_matrix(step, "p2_dia", wfc.get_pow_p(0, 2) ) 
        saver.save_matrix(step, "p2_adi", wfc.get_pow_p(1, 2) ) 

        ncustom_pops = len(params["custom_pops"])
        for pop_type in range(ncustom_pops):

            pop_prms = params["custom_pops"][pop_type]            
            pop_val = wfc.get_pops(pop_prms[0], Py2Cpp_double(pop_prms[1]), Py2Cpp_double(pop_prms[2])) 

            saver.save_multi_matrix(step, pop_type, "custom_pops",  pop_val )

            #print(F"pop_val.num_rows= {pop_val.num_of_rows}  pop_val.num_cols= {pop_val.num_of_cols}")
            #print(F"rep, bmin, bmax = {pop_prms}")
            #print(F"step= {step} pop_type= {pop_type} P0= {pop_val.get(0, 0)}  P1= {pop_val.get(1, 0)} ")


        
        
    if hdf5_output_level>=3:                
        saver.save_matrix(step, "denmat_dia", wfc.get_den_mat(0) ) 
        saver.save_matrix(step, "denmat_adi", wfc.get_den_mat(1) ) 
        
        
    if hdf5_output_level>=4:                
        """
        This is REALLY-REALLY bad option to go since it is soo time-consuming
        """
        
        for ipt in range(wfc.Npts):
            saver.save_multi_matrix(step, ipt, "PSI_dia", wfc.PSI_dia[ipt]) 
            
        for ipt in range(wfc.Npts):            
            saver.save_multi_matrix(step, ipt, "PSI_adi", wfc.PSI_adi[ipt]) 
            
        for ipt in range(wfc.Npts):
            saver.save_multi_matrix(step, ipt, "reciPSI_dia", wfc.reciPSI_dia[ipt]) 
            
        for ipt in range(wfc.Npts):            
            saver.save_multi_matrix(step, ipt, "reciPSI_adi", wfc.reciPSI_adi[ipt])             


            

def save_data_mem(step, wfc, saver, params):
    
    dt = params["dt"]    
    masses = Py2Cpp_double(params["masses"])
    mem_output_level = params["mem_output_level"]
    
    
    if mem_output_level>=1:         
        saver.add_data("timestep", step)
        saver.add_data("time", step * dt)
        saver.add_data("Ekin_dia", wfc.e_kin(masses, 0) )
        saver.add_data("Ekin_adi", wfc.e_kin(masses, 1) )
        saver.add_data("Epot_dia", wfc.e_pot(0) )
        saver.add_data("Epot_adi", wfc.e_pot(1) )
        saver.add_data("Etot_dia", wfc.e_tot(masses, 0) )
        saver.add_data("Etot_adi", wfc.e_tot(masses, 1) )
        saver.add_data("norm_dia", wfc.norm(0) )
        saver.add_data("norm_adi", wfc.norm(1) )

    if mem_output_level>=2:        
        saver.add_data("pop_dia", wfc.get_pops(0) ) 
        saver.add_data("pop_adi", wfc.get_pops(1) )     
        saver.add_data("q_dia", wfc.get_pow_q(0, 1) ) 
        saver.add_data("q_adi", wfc.get_pow_q(1, 1) )         
        saver.add_data("q2_dia", wfc.get_pow_q(0, 2) ) 
        saver.add_data("q2_adi", wfc.get_pow_q(1, 2) )         
        saver.add_data("p_dia", wfc.get_pow_p(0, 1) ) 
        saver.add_data("p_adi", wfc.get_pow_p(1, 1) )         
        saver.add_data("p2_dia", wfc.get_pow_p(0, 2) ) 
        saver.add_data("p2_adi", wfc.get_pow_p(1, 2) )         

    if mem_output_level>=3:        
        saver.add_data("denmat_dia", wfc.get_den_mat(0) )  
        saver.add_data("denmat_adi", wfc.get_den_mat(1) )         
        
    if mem_output_level>=4:                
        saver.add_data("PSI_dia", wfc.PSI_dia) 
        saver.add_data("PSI_adi", wfc.PSI_adi) 
        saver.add_data("reciPSI_dia", wfc.reciPSI_dia) 
        saver.add_data("reciPSI_adi", wfc.reciPSI_adi) 
        


def run_dynamics(wfc, params, model_params, savers):
    
    integrators_map = {"SOFT": 0,
                       "direct_dia": 1,
                       "direct_adi": 2,
                       "Colbert_Miller_dia":3,
                       "Colbert_Miller_adi":4,
                       "Colbert_Miller_SOFT":5
                      }
    integrator_id = integrators_map[ params["integrator"] ]
        
    nsteps = params["nsteps"]      
    print_freq = int(params["progress_frequency"]*nsteps)    
    dt = params["dt"]    
    masses = Py2Cpp_double(params["masses"])    
    
            
    
    #================ Special setups ===================
    expT, expV = None, None
    if integrator_id == 5:
       
        T = wfc.operator_T(Py2Cpp_int([1]), masses, 1.0+0.0j) 
    
        expT = CMATRIX(wfc.Npts, wfc.Npts)
        exp_matrix(expT, T, -dt*1.0j, 0)
    
        expV = CMATRIXList()
        for ipt in range(wfc.Npts):
            expv = CMATRIX(wfc.nstates, wfc.nstates)
            exp_matrix(expv, wfc.Hdia[ipt], -dt*0.5j, 0);        
            expV.append(CMATRIX(expv))
        
        
    #=============== Propagation ==========================
        
    for step in range(nsteps):
        
        #================ Saving the data ==================
        if step%print_freq==0:
            print(F" step= {step}")
            
        # Save properties
        if savers["hdf5_saver"] != None:            
            save_data_hdf5(step, wfc, savers["hdf5_saver"], params)
            
        if savers["txt_saver"] != None:            
            save_data_txt(step, wfc, savers["txt_saver"], params)
            
        if savers["mem_saver"] != None:            
            save_data_mem(step, wfc, savers["mem_saver"], params)
    
        #================ Integration ==================
    
        if integrator_id==0:  # SOFT            
            wfc.SOFT_propagate()      # evolve the diabatic wfc 
            
        elif integrator_id == 1: # direct_dia
            if step==0:            
                wfc.direct_propagate_dia1(dt, masses)                           
            else:
                wfc.direct_propagate_dia2(dt, masses)                          
                
        elif integrator_id == 2: # direct_adi
            if step==0:            
                wfc.direct_propagate_adi1(dt, masses)                           
            else:
                wfc.direct_propagate_adi2(dt, masses)                                          
                
        elif integrator_id == 3: # Colbert-Miller_dia        
            if step==0:            
                wfc.Colbert_Miller_propagate_dia1(dt, masses)                           
            else:
                wfc.Colbert_Miller_propagate_dia2(dt, masses)                          

        elif integrator_id == 4: # Colbert-Miller_adi        
            if step==0:            
                wfc.Colbert_Miller_propagate_adi1(dt, masses)                           
            else:
                wfc.Colbert_Miller_propagate_adi2(dt, masses)                                          

        elif integrator_id == 5: # Colbert-Miller_SOFT        
            wfc.Colbert_Miller_SOFT(expT, expV, 0)    
            
            
        #============= Update other variables ================
        if integrator_id in [0, 1, 3, 5]:
            wfc.update_adiabatic()  
            
        if integrator_id in [2, 4]:
            wfc.update_diabatic()              
        
        wfc.update_reciprocal(0)  # update reci of diabatic function
        wfc.update_reciprocal(1)  # update reci of adiabatic function
        


def run_relaxation(_params, _potential, model_params):
    """

    # For diabatic populations in boxes
    _params["customs_pops"] = [ [0, [-100], [-10] ], 
                                [0, [-10], [10] ],
                                [0, [10], [100] ]
                              ]

    # For adiabatic populations in boxes
    _params["customs_pops"] = [ [1, [-100], [-10] ], 
                                [1, [-10], [10] ],
                                [1, [10], [100] ]
                              ]

    """
    
        
    params = dict(_params)                
    nstates = len(model_params["E_n"])
    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"]*nsteps)    
    params.update({"nstates":nstates})


    ncustom_pops = 0
    if params["custom_pops"]!=None:
        ncustom_pops = len(params["custom_pops"])
        print(F"ncustom_pops = {ncustom_pops}")
            
    
    print("Run calculations with the dynamical parameters ", params)
    print("Run calculations with the model parameters ", model_params)    
    
        
    wfc = init_wfc(params, _potential, model_params )
    ndof = wfc.ndof
    ngrid = wfc.Npts
    
    
    #================ Create savers ==================    
    # Create an output directory, if not present  
    prefix = params["prefix"]
    
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
            
    # Simulation parameters                    
    f = open(F"{prefix}/_dyn_params.txt","w")
    f.write( str(params) );  f.close()
    
    f = open(F"{prefix}/_model_params.txt","w")
    f.write( str(model_params) );  f.close()    

    
    
    properties_to_save = params["properties_to_save"]
    
    #====== HDF5 ========
    hdf5_saver = None
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level > 0:                
        hdf5_saver = dynamics_hdf5.hdf5_saver(F"{prefix}/data.hdf", properties_to_save) 
        hdf5_saver.set_compression_level(params["use_compression"], params["compression_level"])
        dynamics_hdf5.exact_init_hdf5(hdf5_saver, hdf5_output_level, nsteps, ndof, nstates, ngrid)
        dynamics_hdf5.exact_init_custom_hdf5(hdf5_saver, nsteps, ncustom_pops, nstates)  # boxed populations on adiabatic/diabatic states




    #====== TXT ========
    txt_saver = None
    if params["txt_output_level"] > 0:
        pass
    
    #====== MEM =========
    mem_saver = None
    if params["mem_output_level"] > 0:
        mem_saver =  dynamics_hdf5.mem_saver(properties_to_save)
                         
    savers = {"hdf5_saver":hdf5_saver, "txt_saver":txt_saver, "mem_saver":mem_saver }
    

    #==================== Dynamics ======================    
    start = time.time()                               
    
    res = run_dynamics(wfc, params, model_params, savers)                     
    
    end = time.time()    
    print(F"Calculation time = {end - start} seconds")
    
    
    if mem_saver != None:        
        return mem_saver

        


def plot_mem(res, _params, model_params, plot_params):
    """
    Args:
        res ( mem_saver ): the object with all the suitable parameters
        _params ( dict ): simulation control parameters
        model_params ( dict ): the parameters of the model we compute
        
    """

    colors = {}

    colors.update({"11": "#8b1a0e"})  # red       
    colors.update({"12": "#FF4500"})  # orangered 
    colors.update({"13": "#B22222"})  # firebrick 
    colors.update({"14": "#DC143C"})  # crimson   

    colors.update({"21": "#5e9c36"})  # green
    colors.update({"22": "#006400"})  # darkgreen  
    colors.update({"23": "#228B22"})  # forestgreen
    colors.update({"24": "#808000"})  # olive      

    colors.update({"31": "#8A2BE2"})  # blueviolet
    colors.update({"32": "#00008B"})  # darkblue  

    colors.update({"41": "#2F4F4F"})  # darkslategray

    clrs_index = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]



    # Parameters and dimensions
    critical_params = [  ] 
    default_params = {  "colors":colors, "clrs_index":clrs_index }
    comn.check_input(plot_params, default_params, critical_params)
        
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]



    
    
    params = dict(_params)
    
    nsteps = params["nsteps"]    
    nstates = len(model_params["E_n"])
    prefix = params["prefix"]
    properties_to_save = params["properties_to_save"]
 
    
    t = None
    if "time" in properties_to_save:
        t = res.data["time"]
    
    #=============== Populations ======================
    
    plt.figure(1, figsize=(36, 12)) # dpi=300, frameon=False)    
    plt.subplot(1, 2, 1)
    plt.title('Adiabatic population dynamics' )
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    
    if "pop_adi" in properties_to_save and t != None:
        for i in range(nstates):        
            Pi = data_conv.unpack1(res.data["pop_adi"], i, 0, 0)                
            plt.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
            plt.legend()
    
    
    plt.subplot(1, 2, 2)
    plt.title('Diabatic population dynamics' )
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    
    if "pop_dia" in properties_to_save and t != None:
        for i in range(nstates):
            Pi = data_conv.unpack1(res.data["pop_dia"], i, 0, 0)                        
            plt.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
            plt.legend()
            
    plt.savefig("%s/Fig1.png" % (prefix), dpi=300)
    plt.savefig("%s/Fig1.pdf" % (prefix), dpi=300)
            
                
    #============= Energies =====================
    plt.figure(2, figsize=(36, 12)) # dpi=300, frameon=False)           
    
    plt.subplot(1, 2, 1)
    plt.title('Energies' )
    plt.xlabel('t, a.u.')
    plt.ylabel('Energy, a.u.')
    if "Ekin_dia" in properties_to_save \
       and "Epot_dia" in properties_to_save \
       and "Etot_dia" in properties_to_save \
       and t != None:
        
        Ekin_dia = res.data["Ekin_dia"]  
        Epot_dia = res.data["Epot_dia"]  
        Etot_dia = res.data["Etot_dia"]  
        plt.plot(t, Etot_dia, label='$Etot_{dia}$', linewidth=10, color = colors["11"])   
        plt.plot(t, Ekin_dia, label='$Ekin_{dia}$', linewidth=10, color = colors["21"])   
        plt.plot(t, Epot_dia, label='$Epot_{dia}$', linewidth=10, color = colors["31"])   
        plt.legend()
    
    plt.subplot(1, 2, 2)
    plt.title('Energies' )
    plt.xlabel('t, a.u.')
    plt.ylabel('Energy, a.u.')
    
    if "Ekin_adi" in properties_to_save \
       and "Epot_adi" in properties_to_save \
       and "Etot_adi" in properties_to_save \
       and t != None:
        
        Ekin_adi = res.data["Ekin_adi"]  
        Epot_adi = res.data["Epot_adi"]  
        Etot_adi = res.data["Etot_adi"]  
        plt.plot(t, Etot_adi, label='$Etot_{adi}$', linewidth=10, color = colors["11"])   
        plt.plot(t, Ekin_adi, label='$Ekin_{adi}$', linewidth=10, color = colors["21"])   
        plt.plot(t, Epot_adi, label='$Epot_{adi}$', linewidth=10, color = colors["31"])   
        plt.legend()

    plt.savefig("%s/Fig2.png" % (prefix), dpi=300)
    plt.savefig("%s/Fig2.pdf" % (prefix), dpi=300)
    
  
    #============= Phase spaces & Norms  =====================
    plt.figure(3, figsize=(36, 12)) # dpi=300, frameon=False)           
                
    plt.subplot(1, 2, 1)
    plt.title('Phase space' )
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Momentum, a.u.')   
    
    if "q_dia" in properties_to_save and "p_dia" in properties_to_save:
        ndof = res.data["q_dia"][0].num_of_rows
    
        for idof in range(ndof):
            qi = data_conv.unpack1(res.data["q_dia"], idof, 0, 0)
            pi = data_conv.unpack1(res.data["p_dia"], idof, 0, 0)
    
            plt.plot(qi, pi, label='', linewidth=10, color = colors[clrs_index[i]])   
            plt.legend()
                    

    plt.subplot(1, 2, 2)
    plt.title('Norms' )
    plt.xlabel('Time, a.u.')
    plt.ylabel('Norm')
    
    if "norm_dia" in properties_to_save and "norm_adi" in properties_to_save and t != None:
        
        nrm_dia = res.data["norm_adi"]
        nrm_adi = res.data["norm_dia"]
                
        plt.plot(t, nrm_dia, label='Diabatic', linewidth=10, color = colors["11"])   
        plt.plot(t, nrm_adi, label='Adiabatic', linewidth=10, color = colors["21"])   
        plt.legend()

    plt.savefig("%s/Fig3.png" % (prefix), dpi=300)
    plt.savefig("%s/Fig3.pdf" % (prefix), dpi=300)
    
        
    plt.show()
    plt.close()



def plot_hdf5(plot_params):
    """
    This function is meant to plot the results stored in the hdf files generated by the exact dynamics runs

    Args:

        prefix ( string ): the name of the directory containing the input HDF5 file
            This directory will also be used to output the generated picture files [ default : "out"]
        filename ( string ): name of the HDF5 file to read [ default: "data.hdf"]
        output_level ( int ): the level of info contained in the HDF5 file [ default : 3]        
        which_adi_states ( list of ints ) : indices of the adiabatic states to print [ default: [0] ]
        which_dia_states ( list of ints ) : indices of the diabatic states to print [ default: [0] ]
        colors ( dictionary ): the definition of the colors to use
        clrs_index ( list of strings ) : defines the mapping of the colors on integers and vice versa 
    

    """
    
    plt.rc('axes', titlesize=24)      # fontsize of the axes title
    plt.rc('axes', labelsize=20)      # fontsize of the x and y labels
    plt.rc('legend', fontsize=20)     # legend fontsize
    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels

    plt.rc('figure.subplot', left=0.2)
    plt.rc('figure.subplot', right=0.95)
    plt.rc('figure.subplot', bottom=0.13)
    plt.rc('figure.subplot', top=0.88)


    colors = {}

    colors.update({"11": "#8b1a0e"})  # red       
    colors.update({"12": "#FF4500"})  # orangered 
    colors.update({"13": "#B22222"})  # firebrick 
    colors.update({"14": "#DC143C"})  # crimson   

    colors.update({"21": "#5e9c36"})  # green
    colors.update({"22": "#006400"})  # darkgreen  
    colors.update({"23": "#228B22"})  # forestgreen
    colors.update({"24": "#808000"})  # olive      

    colors.update({"31": "#8A2BE2"})  # blueviolet
    colors.update({"32": "#00008B"})  # darkblue  

    colors.update({"41": "#2F4F4F"})  # darkslategray

    clrs_index = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]



    # Parameters and dimensions
    critical_params = [  ] 
    default_params = {  "prefix":"out", "filename":"data.hdf", "hdf5_output_level":3,
                        "which_dofs":[0], "which_adi_states":[0], "which_dia_states":[0],
                        "colors":colors, "clrs_index":clrs_index,
                        "properties_to_save": 
                          [ "timestep", "time", "Ekin_dia", "Ekin_adi", "Epot_dia", 
                            "Epot_adi", "Etot_dia", "Etot_adi", "norm_dia", "norm_adi",
                            "pop_dia", "pop_adi", "q_dia", "q_adi", "q2_dia", "q2_adi", 
                            "p_dia", "p_adi", "p2_dia", "p2_adi",
                            "denmat_dia", "denmat_adi", 
                            "PSI_dia", "PSI_adi", "reciPSI_dia", "reciPSI_adi" ] 
                     }
    comn.check_input(plot_params, default_params, critical_params)
        

            
    filename = plot_params["filename"]
    prefix = plot_params["prefix"]
    hdf5_output_level = plot_params["hdf5_output_level"]    
    which_dofs = plot_params["which_dofs"]
    which_adi_states = plot_params["which_adi_states"]
    which_dia_states = plot_params["which_dia_states"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    properties_to_save = plot_params["properties_to_save"]

    out_prefix = prefix
    

    with h5py.File(F"{prefix}/{filename}", 'r') as f:
    
        t = None
        if "time" in properties_to_save:                        
            t = list(f["time/data"][:])
    
        #=============== Populations ======================
    
        plt.figure(1, figsize=(36, 12)) # dpi=300, frameon=False)    
        plt.subplot(1, 2, 1)
        plt.title('Adiabatic population dynamics' )
        plt.xlabel('Time, a.u.')
        plt.ylabel('Population')
    
        if "pop_adi" in properties_to_save and t != None:
            nstates = f["pop_adi"].attrs['dim'][1]                        
            for i in range(nstates):        
                if i in which_adi_states:
                    Pi = list(f["pop_adi/data"][:, i, 0])                
                    plt.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
                    plt.legend()
    
    
        plt.subplot(1, 2, 2)
        plt.title('Diabatic population dynamics' )
        plt.xlabel('Time, a.u.')
        plt.ylabel('Population')
    
        if "pop_dia" in properties_to_save and t != None:
            nstates = f["pop_dia"].attrs['dim'][1]                        
            for i in range(nstates):
                if i in which_dia_states:
                    Pi = list(f["pop_dia/data"][:, i, 0])                            
                    plt.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
                    plt.legend()
            
        plt.savefig("%s/Fig1.png" % (prefix), dpi=300)
        plt.savefig("%s/Fig1.pdf" % (prefix), dpi=300)
            
                
        #============= Energies =====================
        plt.figure(2, figsize=(36, 12)) # dpi=300, frameon=False)           
    
        plt.subplot(1, 2, 1)
        plt.title('Energies' )
        plt.xlabel('t, a.u.')
        plt.ylabel('Energy, a.u.')
        if "Ekin_dia" in properties_to_save \
           and "Epot_dia" in properties_to_save \
           and "Etot_dia" in properties_to_save \
           and t != None:
        
            Ekin_dia =  list(f["Ekin_dia/data"][:])                             
            Epot_dia =  list(f["Epot_dia/data"][:])                             
            Etot_dia =  list(f["Etot_dia/data"][:])                             
        
            plt.plot(t, Etot_dia, label='$Etot_{dia}$', linewidth=10, color = colors["11"])   
            plt.plot(t, Ekin_dia, label='$Ekin_{dia}$', linewidth=10, color = colors["21"])   
            plt.plot(t, Epot_dia, label='$Epot_{dia}$', linewidth=10, color = colors["31"])   
            plt.legend()
    
    
        plt.subplot(1, 2, 2)
        plt.title('Energies' )
        plt.xlabel('t, a.u.')
        plt.ylabel('Energy, a.u.')
    
        if "Ekin_adi" in properties_to_save \
           and "Epot_adi" in properties_to_save \
           and "Etot_adi" in properties_to_save \
           and t != None:
        
            Ekin_adi =  list(f["Ekin_adi/data"][:])                             
            Epot_adi =  list(f["Epot_adi/data"][:])                             
            Etot_adi =  list(f["Etot_adi/data"][:])                             
        
            plt.plot(t, Etot_adi, label='$Etot_{adi}$', linewidth=10, color = colors["11"])   
            plt.plot(t, Ekin_adi, label='$Ekin_{adi}$', linewidth=10, color = colors["21"])   
            plt.plot(t, Epot_adi, label='$Epot_{adi}$', linewidth=10, color = colors["31"])   
            plt.legend()

        plt.savefig("%s/Fig2.png" % (prefix), dpi=300)
        plt.savefig("%s/Fig2.pdf" % (prefix), dpi=300)
    
  
        #============= Phase spaces & Norms  =====================
        plt.figure(3, figsize=(36, 12)) # dpi=300, frameon=False)           
                
        plt.subplot(1, 2, 1)
        plt.title('Phase space' )
        plt.xlabel('Coordinate, a.u.')
        plt.ylabel('Momentum, a.u.')   
    
        if "q_dia" in properties_to_save and "p_dia" in properties_to_save:                                
            ndof = f["q_dia"].attrs['dim'][1]                        
    
            for idof in range(ndof):
                if idof in which_dofs:
                    qi =  list(f["q_dia/data"][:, idof, 0])                             
                    pi =  list(f["p_dia/data"][:, idof, 0])                                         
    
                    plt.plot(qi, pi, label='', linewidth=10, color = colors[clrs_index[i]])   
                    plt.legend()
                    
                
        plt.subplot(1, 2, 2)
        plt.title('Norms' )
        plt.xlabel('Time, a.u.')
        plt.ylabel('Norm')
    
        if "norm_dia" in properties_to_save and "norm_adi" in properties_to_save and t != None:
        
            nrm_dia = list(f["norm_dia/data"][:])                             
            nrm_adi = list(f["norm_adi/data"][:])                                     
                
            plt.plot(t, nrm_dia, label='Diabatic', linewidth=10, color = colors["11"])   
            plt.plot(t, nrm_adi, label='Adiabatic', linewidth=10, color = colors["21"])   
            plt.legend()

        plt.savefig("%s/Fig3.png" % (prefix), dpi=300)
        plt.savefig("%s/Fig3.pdf" % (prefix), dpi=300)
    
        
        plt.show()
        plt.close()
        
