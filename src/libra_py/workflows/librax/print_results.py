#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file print_results.py
# This module implements the functions which print out the results of NA-MD calculation.

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def print_one_traj(isnap, iconf, i_ex, itraj, mol, syst, mu, epot, ekin, etot, eext, params):
    # This function prints out the results for only one trajectory.
    #
    # \param[in] isnap   snap index 
    # \param[in] iconf   initial geometry index
    # \param[in] i_ex    initial excitation index
    # \param[in] itraj   TSH trajectory index
    # \param[in] mol     list of objects containing Nuclear DOF's whose size is nconfig*nstates*num_SH_traj. Lists below have the same size as this.
    # \param[in] syst    list of objects containing atomic system information
    # \param[in] mu      is a transition dipole moment matrix (MATRIX object) of a single trajectory
    # \param[in] epot    current electronic state potential energy for one trajectory, so this is a floating point value
    # \param[in] ekin    current kinetic energy of the system for one trajectory
    # \param[in] etot    = epot + ekin for one trajectory
    # \param[in] eext    = etot + (Thermostat energy) for one trajectory
    # \param[in] params  list of input parameters from {gms,qe,g09}_run.py
    #                                                                                                                                                        
    # Used in:  print_results.py/print_ens_traj

    kB = 3.166811429e-6 # Boltzmann constant in hartree unit                                                                                                 
    dt_nucl = params["dt_nucl"]
    #flag_ao = params["flag_ao"]
    MD_type = params["MD_type"]
    print_coherences = params["print_coherences"]
    nstates = len(params["excitations"])
    nstates_init = len(params["excitations_init"])
    ij = isnap*params["Nsteps"]
    
    traj_file_prefix = params["res"]+"md"
    ene_file_prefix = params["res"]+"ene"
    mu_file_prefix = params["res"]+"mu"

    # Re-compute energies, to print
    # WE DON'T WANT TO CALL compute_potential_energy here!
    #    if params["interface"] == "GAMESS":
    #        epot = tot_ene + compute_potential_energy(mol, el, ham, f_pot)
    #    elif params["interface"] == "QE":
    #        epot = tot_ene.get(i_ex,i_ex) + compute_potential_energy(mol, el, ham, f_pot)
    #    print "epot = ", epot    
    #    ekin = compute_kinetic_energy(mol)

    curr_T = 2.0*ekin/(3.0*syst.Number_of_atoms*kB)

    # set file name
    num_tmp = "_"+str(iconf)+"_"+str(i_ex)+"_"+str(itraj)
    ene_file = ene_file_prefix+num_tmp+".txt"
    traj_file = traj_file_prefix+num_tmp+".xyz"
    mu_file = mu_file_prefix+num_tmp+".txt"

    ##print 
    # Geometry
    syst.set_atomic_q(mol.q)
    syst.print_xyz(traj_file,ij)

    # Energy
    fe = open(ene_file,"a")
    fe.write("t= %8.5f ekin= %10.7f  epot= %10.7f  etot= %10.7f  eext= %10.7f curr_T= %8.5f\n" % (ij*dt_nucl, ekin, epot, etot, eext,curr_T))
    fe.close()
    
    if params["interface"] == "GAMESS":       
        # Dipole moment components
        fm = open(mu_file,"a")
        line = "t= %8.5f " % (ij*dt_nucl)
        # *************************************************************************************
        Nao = mu[0].num_of_rows
        for k in xrange(Nao):
            line = line + " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f " % (mu[0].real().get(k,k),mu[0].imag().get(k,k),\
                                                                         mu[1].real().get(k,k),mu[1].imag().get(k,k),\
                                                                         mu[2].real().get(k,k),mu[2].imag().get(k,k))
        # commented out temporarily because mu[k] is a complex number; will be modified later
        # **************************************************************************************

        line = line + "\n"
        fm.write(line)
        fm.close()


def print_ens_traj(isnap,mol,syst,mu,epot,ekin,etot,eext,params):
    # This function prints out auxiliary results for an ensemble of trajectories:
    # one trajctory contains the index of initial geometry, initial excitation state
    # (TSH trajectory, optionally)
    #
    # \param[in] isnap   snap index
    # \param[in] mol     list of objects containing Nuclear DOF's whose size is nconfig*nstates*num_SH_traj. Lists below have the same size as this.
    # \param[in] syst    list of objects containing atomic system information
    # \param[in] mu      list of transition dipole moments (MATRIX objects)
    # \param[in] epot    list of the current electronic state potential energies for each trajectory in the ensemble, so this is a list of N floating point values
    # \param[in] ekin    list of the current kinetic energies of the system
    # \param[in] etot    list of  (epot + ekin)
    # \param[in] eext    list of (etot + Thermostat energy)
    # \param[in] params  list of input parameters from {gms,qe,g09}_run.py
    #
    # Used in:  md.py/run_MD 

    nconfig = params["nconfig"]
    nstates = len(params["excitations"])
    nstates_init = len(params["excitations_init"])
    num_SH_traj = params["num_SH_traj"]

    for iconf in xrange(nconfig):
        for i_ex in xrange(nstates_init):
            cnt = iconf*nstates_init + i_ex
            for itraj in xrange(num_SH_traj):

                cnt_inc_el = iconf*nstates_init*num_SH_traj + i_ex*num_SH_traj + itraj
                if params["do_rescaling"] ==1:
                    cnt = cnt_inc_el
                print_one_traj(isnap,iconf,params["excitations_init"][i_ex],itraj,mol[cnt],syst[cnt],mu[cnt],epot[cnt], ekin[cnt], etot[cnt], eext[cnt],params)



def pops_ave_TSH_traj(i,el,params):
    # This function averages SE and SH populations over TSH trajectories.
    #
    # \param[in] i       snap index
    # \param[in] el      object containing electronic DOF's
    # \param[in] params  list of input parameters from run.py 
    # l_se_pop - list of SE populations whose size is nconfig x nstates. 
    # l_sh_pop - list of SH populations (the same size above)
    # Used in: md.py/run_MD

    nconfig = params["nconfig"]
    nstates = el[0].nstates
    nstates_init = len(params["excitations_init"])
    num_SH_traj = params["num_SH_traj"]
    SH_type = params["tsh_method"]
    ij = i*params["Nsteps"]
    dt_nucl = params["dt_nucl"]
    print_coherences = params["print_coherences"]

    # define prefixes
    se_pop_file_prefix = params["res"]+"se_pop"
    sh_pop_file_prefix = params["res"]+"sh_pop"

    l_se_pop = []
    l_sh_pop = []

    for iconf in xrange(nconfig):
        for i_ex in xrange(nstates_init):

            i_ex_init = params["excitations_init"][i_ex]
            
            #***********  SE population **************

            se_pop = []
            se_coh_re = []
            se_coh_im = []

            # average SE populations
            for st in xrange(nstates):
                p_tmp = 0.0
                for itraj in xrange(num_SH_traj):
                    cnt = iconf*nstates_init*num_SH_traj + i_ex*num_SH_traj + itraj
                    p_tmp += el[cnt].rho(st,st).real
                se_pop.append(p_tmp/float(num_SH_traj))
                l_se_pop.append(p_tmp/float(num_SH_traj)) # store SE pops

            # average coherences
            if print_coherences == 1:
                for st in xrange(nstates):
                    for st1 in xrange(st):
                        c1_tmp = 0.0; c2_tmp = 0.0
                        for itraj in xrange(num_SH_traj):
                            cnt = iconf*nstates_init*num_SH_traj + i_ex*num_SH_traj + itraj
                            c1_tmp += el[cnt].rho(st,st1).real
                            c2_tmp += el[cnt].rho(st,st1).imag                            
                        se_coh_re.append(c1_tmp/float(num_SH_traj))
                        se_coh_im.append(c2_tmp/float(num_SH_traj))

            # set file name
            num_tmp = "_"+str(iconf)+"_"+str(i_ex_init)
            se_pop_file = se_pop_file_prefix+num_tmp+".txt"

            # Populations
            fel = open(se_pop_file,"a")

            # Print time
            line_se = "t= %8.5f " % (ij*dt_nucl)

            # Print populations
            for st in xrange(nstates):
                line_se = line_se + " %8.5f " % se_pop[st]

            if print_coherences == 1:
                # Print coherences
                cnt = 0
                for st in xrange(nstates):
                    for st1 in xrange(st):
                        line_se = line_se + " %8.5f %8.5f " % (se_coh_re[cnt], se_coh_im[cnt])
                        cnt += 1

            line_se = line_se + "\n"

            fel.write(line_se)
            fel.close()

            #***********  SH population **************
            if SH_type>=1:

                # evaluate TSH probabilities
                tsh_probs = [0.0]*nstates
                for itraj in xrange(num_SH_traj): # count number of trajectories
                    cnt = iconf*nstates_init*num_SH_traj + i_ex*num_SH_traj + itraj
                    tsh_probs[el[cnt].istate] += 1.0
            
                for st in xrange(nstates):
                    tsh_probs[st] = tsh_probs[st]/float(num_SH_traj)
                    l_sh_pop.append(tsh_probs[st]) # store SH pops

                # set file name
                num_tmp = "_"+str(iconf)+"_"+str(i_ex_init)
                sh_pop_file = sh_pop_file_prefix+num_tmp+".txt"

                # Populations     
                fel1 = open(sh_pop_file,"a")

                # Print time      
                line_sh = "t= %8.5f " % (ij*dt_nucl)

                #Print populations
                for st in xrange(nstates):
                    line_sh = line_sh + " %8.5f " % (tsh_probs[st])

                line_sh = line_sh + "\n"

                fel1.write(line_sh)
                fel1.close()

    return l_se_pop, l_sh_pop


def pops_ave_geometry(i,nstates,se_pop,sh_pop,params):
    # This function prints out SE and SH populations averaged over initial geometry
    # 
    # \param[in] i       snap index
    # \param[in] nstates number of excitation states
    # \param[in] se_pop  list of SE populations averaged over TSH trajectories; the size is nconfig*nstates
    # \param[in] sh_pop  list of SH populations averaged over TSH trajectories; the size is nconfig*nstates 
    # \param[in] params  list of input parameters from {gms,qe,g09}_run.py 
    #
    # Used in: md.py/run_MD

    nstates_init = len(params["excitations_init"])
    nconfig = params["nconfig"]
    SH_type = params["tsh_method"]
    ij = i*params["Nsteps"]
    dt_nucl = params["dt_nucl"]

    # define prefixes
    se_pop_ex_file_prefix = params["res"]+"se_pop_ex"
    sh_pop_ex_file_prefix = params["res"]+"sh_pop_ex"


    #print "length of se_pop is %d" %(len(se_pop))
    #print "length of sh_pop is %d" %(len(sh_pop))

    for i_ex in xrange(nstates_init):

        i_ex_init = params["excitations_init"][i_ex]

        # set file name                                                                                                                               
        se_pop_file = se_pop_ex_file_prefix+str(i_ex_init)+".txt"
        sh_pop_file = sh_pop_ex_file_prefix+str(i_ex_init)+".txt"

        # Populations                                                                                                            
        fse = open(se_pop_file,"a")
        fsh = open(sh_pop_file,"a")

        # Print time                                                                                                                                  
        line_se = "t= %8.5f " % (ij*dt_nucl)
        line_sh = "t= %8.5f " % (ij*dt_nucl)

        for st in xrange(nstates):

            se_sum = 0.0; sh_sum = 0.0;
            for iconf in xrange(nconfig):
                cnt = iconf*nstates_init*nstates + i_ex*nstates + st

                se_sum += se_pop[cnt]
                sh_sum += sh_pop[cnt]

            #Print populations
            line_se = line_se + " %8.5f " % (se_sum/float(nconfig))
            line_sh = line_sh + " %8.5f " % (sh_sum/float(nconfig))

        line_se = line_se + "\n"
        line_sh = line_sh + "\n"

        fse.write(line_se); fse.close();
        fsh.write(line_sh); fsh.close();
