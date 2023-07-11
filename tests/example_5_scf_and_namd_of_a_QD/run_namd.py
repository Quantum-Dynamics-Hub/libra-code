#*********************************************************************************
#* Copyright (C) 2017-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  This code demonstarates how to run atomistic SCF or NA-MD calculations
#
#
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *



def add_printout(i, pop, tr_adi, tr_dia, filename):
# pop - CMATRIX(nstates, 1)

    f = open(filename,"a")
    line = "step= %4i " % i    

    tot_pop = 0.0
    for st in xrange(pop.num_of_cols):
        pop_o = pop.get(st,st).real
        tot_pop = tot_pop + pop_o
        line = line + " P(%4i)= %8.5f " % (st, pop_o)
    line = line + " Total= %8.5f  Tr(|adi><adi|)= %8.5f Tr(|dia><dia|)= %8.5f\n" % (tot_pop , tr_adi, tr_dia)
    f.write(line)
    f.close()



def add_printout2(i, denmat_se, denmat_sh, H, filename):
# denmat_se, denmat_sh - Density matrices for SE and SH  CMATRIX(nstates, nstates)
#  H - vibronic Hamiltonian CMATRIX(nstates, nstates)
# filename - the output file name

    f = open(filename,"a")

    ntraj = len(H)
    ene_se, ene_sh = 0.0, 0.0
    ene_se2, ene_sh2 = 0.0, 0.0

    for tr in xrange(ntraj):
        ene_se = ene_se + (denmat_se[tr] * CMATRIX(H[tr].real()) ).tr().real
        ene_sh = ene_sh + (denmat_sh[tr] * CMATRIX(H[tr].real()) ).tr().real

        ene_se2 = ene_se2 + (denmat_se[tr] * CMATRIX(H[tr].real()) ).tr().real / denmat_se[tr].tr().real
        ene_sh2 = ene_sh2 + (denmat_sh[tr] * CMATRIX(H[tr].real()) ).tr().real / denmat_sh[tr].tr().real

    ene_se = ene_se/float(ntraj)
    ene_sh = ene_sh/float(ntraj)

    ene_se2 = ene_se2/float(ntraj)
    ene_sh2 = ene_sh2/float(ntraj)


    line = "step= %4i Ene_SE= %8.5f  Ene_SH= %8.5f Ene_SE_norm= %8.5f  Ene_SH_norm= %8.5f \n" % (i, ene_se, ene_sh, ene_se2, ene_sh2)
    f.write(line)
    f.close()




        
def run_scf(syst, min_, max_):

    T = 278.0
    sigma = 0.1
    rnd = Random()

    ham = listHamiltonian_QM("control_parameters_eht.dat", syst)
    ham.compute_scf(syst)

    es = ham.get_electronic_structure() 
    homo = es.Nocc_alp - 1  # index of the HOMO orbital in the adiabatic sub-system

    print "Adiabatic sub-system: number of electrons (alpha) = ", es.Nocc_alp, " homo = ", homo       
    print "Electronic structure of the adiabatic sub-system"
    print "Index  Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
    for j in xrange(es.Norb):
        print "%5i  %12.8f   %12.8f  %12.8f   %12.8f" %(j, es.get_bands_alp(j), es.get_occ_alp(j), es.get_bands_bet(j), es.get_occ_bet(j) )

    orb_adi = range(homo+min_, homo+max_+1)
    
    scf.spectrum(ham, "T_mo.dat", "spectrum.txt")
    scf.print_orbitals(ham, syst, orb_adi, "init", [40,40,40])
    scf.print_pdos(ham, syst, [ ["s",range(0,syst.Number_of_atoms)], 
                                ["p",range(0,syst.Number_of_atoms)], 
                                ["d",range(0,syst.Number_of_atoms)], 
                              ], "pdos" ) 

   


def run_namd(act_space, SD_basis, init_st):
    """
    act_space - list of integers, relative to HOMO (determined automatically)
    SD_basis - in terms of the active space
    init_st - in terms of SD basis states
    """

    #--------------------- Parameters ----------------------

    T = 278.0  # K
    sigma = 0.1  # 

    do_collapse = 0  # 0 - no decoherence, 1 - decoherence
    sh_method = 1    # 0 - MSSH,  1 - FSSH
    num_sh_traj = 1

    nst_adi = len(SD_basis)
        
    #--------------------- Initialization ----------------------

    rnd = Random()

    # Create molecular system and initialize the properties
    syst = init_ensembles.init_systems2(num_sh_traj, "CdSe-qd.ent", "true_pdb", rnd, T, sigma, "elements.dat")

    # Create a force field and then MM Hamiltonians and connect them to the systems
    uff = ForceField({"bond_functional":"Harmonic", "angle_functional":"Fourier",
                      "dihedral_functional":"General0", "oop_functional":"Fourier",
                      "mb_functional":"LJ_Coulomb","R_vdw_on":40.0,"R_vdw_off":55.0 })
    LoadUFF.Load_UFF(uff,"uff.dat")

    ham_mm = init_ensembles.init_mm_Hamiltonians(syst, uff)

    # Electronic. Nuclear DOFs and initialize them
    mol = init_ensembles.init_mols(syst, num_sh_traj, 3*syst[0].Number_of_atoms)

    el = []
    for tr in xrange(num_sh_traj):
        el.append( Electronic(1,0) )
        nve_md.nve_md_init(syst[tr], mol[tr], el[tr], ham_mm[tr])

    # Thermostat
    tparams = {"Temperature":T,"Q":100.0,"thermostat_type":"Nose-Hoover","nu_therm":0.001,"NHC_size":5}
    therm = init_ensembles.init_therms(num_sh_traj, 3*syst[0].Number_of_atoms, tparams)



    # Initialize the QM Hamiltonian at two timesteps
    prms = Control_Parameters(); 
    control_filename = "control_parameters_eht.dat"
    #control_filename = "control_parameters_indo.dat"


    dt = 20.0

    Hvib2_ave = CMATRIX(nst_adi, nst_adi)
    ham_cur, ham_old, Hvib = [], [], []
    active_space = []


    for tr in xrange(num_sh_traj):
        ham_cur.append( listHamiltonian_QM(control_filename, syst[tr] ) ) # this includes the overlap calculation
        ham_cur[tr].compute_scf(syst[tr])                  # this includes the core Hamiltonian calculations

        ham_old.append( listHamiltonian_QM(ham_cur[tr]) )

        es = ham_cur[tr].get_electronic_structure() 

        homo = es.Nocc_alp - 1  # index of the HOMO orbital in the adiabatic sub-system
        print "Adiabatic sub-system: number of electrons (alpha) = ", es.Nocc_alp, " homo = ", homo       

        if tr==0:
            for it in act_space:
                active_space.append( homo + it )


        Hvib.append( namd.compute_Hvib_sd(ham_old[tr], ham_cur[tr], active_space, SD_basis, dt) )



        if tr==num_sh_traj-1:        
            print "Hvib = ";  Hvib[tr].show_matrix();        
            print "Electronic structure of the adiabatic sub-system"
            print "Index  Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
            for j in xrange(es.Norb):
                print "%5i  %12.8f   %12.8f  %12.8f   %12.8f" %(j, es.get_bands_alp(j), es.get_occ_alp(j), es.get_bands_bet(j), es.get_occ_bet(j) )
    


    #=================== Propagation ====================
    ########################## Cooling #################################

    md = MD({"max_step":100,"ensemble":"NVE","integrator":"DLML","terec_exp_size":10,"dt":10.0,"n_medium":1,"n_fast":1,"n_outer":1}); 
    md.show_info()


    # State: system + thermostat + (optional barostat)
#    anneal_schedule = [ {"dt":20.0, "nsteps":5, "ncycles":20}] # for test (very fast)
#    anneal_schedule = [ {"dt":20.0, "nsteps":50, "ncycles":20}] # for test
    anneal_schedule = [ {"dt":10.0, "nsteps":5, "ncycles":20},  {"dt":10.0, "nsteps":50, "ncycles":20},  {"dt":10.0, "nsteps":50, "ncycles":20} ]


    ST = []
    for tr in xrange(num_sh_traj):
        ST.append( State() )
        ST[tr].set_system(syst[tr]); ST[tr].set_thermostat(therm[tr]); ST[tr].set_md(md)
        ST[tr].init_md(mol[tr], el[tr], ham_mm[tr], rnd)

        f = open("_en_cooling-%i.txt" % tr,"w"); f.close()

        for item in anneal_schedule:
            md.dt = item["dt"]; md.max_step = item["nsteps"] 
        
            for i in xrange(item["ncycles"]):
                syst[tr].print_xyz("_mol_cooling-%i.xyz" % tr,i)
                ST[tr].run_md(el[tr], ham_mm[tr])
                ekin = ST[tr].E_kin; epot = ST[tr].E_pot
                ST[tr].cool()
        
                f = open("_en_cooling-%i.txt" % tr,"a")
                f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ekin, epot, ST[tr].E_tot, ST[tr].H_NP, ST[tr].curr_T ))
                f.close()



    ########################## Production MD - thermalization #################################

    for tr in xrange(num_sh_traj):
        syst[tr].init_atom_velocities(T, rnd)  # must be this!!!

        md.max_step = 10;  md.ensemble = "NVT"; md.dt = 40.0;
        f = open("_en_md-therm-%i.txt" % tr,"w"); f.close()

        for i in xrange(50):  # for test
            syst[tr].print_xyz("_mol_md-therm-%i.xyz" % tr,i)
            ST[tr].run_md(el[tr], ham_mm[tr])

            f = open("_en_md-therm-%i.txt" % tr,"a")
            f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ST[tr].E_kin, ST[tr].E_pot, ST[tr].E_tot, ST[tr].H_NP, ST[tr].curr_T ))
            f.close()


    ########################## Production MD - NA-MD #################################

    #====================== Initialize ================================
    ist_adi = []
    Coeff_adi = []
    denmat_adi_se = []
    denmat_adi_sh = []

       
    for tr in xrange(num_sh_traj):
        f = open("_en_md-%i.txt" % tr,"w"); f.close()

        md.max_step = 1;  md.ensemble = "NVE";  #md.dt = 20.0; 
        md.dt = 40.0
  

        ist_adi.append(init_st)
        Coeff_adi.append(CMATRIX(nst_adi, 1)); Coeff_adi[tr].set(init_st, 1.0, 0.0)
        denmat_adi_se  = tsh.amplitudes2denmat(Coeff_adi)
        denmat_adi_sh.append(CMATRIX(nst_adi, nst_adi)); denmat_adi_sh[tr].set(init_st, init_st, 1.0, 0.0)


    f = open("_pop_adi_sh.txt","w");    f.close()
    f = open("_pop_adi_se.txt","w");    f.close()
    f = open("_Hvib.txt","w");          f.close()
    f = open("_energy-adi.txt","w");    f.close()



    #=======================================================================    
    print "The SCF calculations at the thermalized structure!"
    min_ = min(act_space)
    max_ = max(act_space)
    run_scf(syst[0], min_, max_)  

    for i in xrange(1000):        

        t = Timer();  t.start()
        for tr in xrange(num_sh_traj):
            syst[tr].print_xyz("_mol_md-%i.xyz" % tr, i) 
            ST[tr].run_md(el[tr], ham_mm[tr])  # evolve q

            #=============== Now handle the electronic structure ==========================
            ham_old[tr] = listHamiltonian_QM(ham_cur[tr])        
        
            # Update positions of the basis AOs!!!
            for o in xrange(ham_cur[tr].Norb):
                at_indx = ham_cur[tr].ao_to_atom_map[o]  # indices of the atoms in the system bound to this Hamiltonin object, that
                                                         # is syst_adi !!!  so we will need another mapping of the indices of those atoms to
                                                         # the indices to the overall system atoms
                ham_cur[tr].basis_ao[o].set_position(syst[tr].Atoms[at_indx].Atom_RB.rb_cm)


            # Update overlaps and core Hamiltonian
            ham_cur[tr].compute_overlap(syst[tr])
            ham_cur[tr].compute_scf(syst[tr])       # this includes the core Hamiltonian calculations

            Hvib[tr] = namd.compute_Hvib_sd(ham_old[tr], ham_cur[tr], active_space, SD_basis, dt * md.max_step)
            
            datautils.show_matrix_splot(Hvib[tr].imag(), "_Hvib_im.dat")
            datautils.show_matrix_splot(Hvib[tr].real(), "_Hvib_re.dat")


            # Compute the running average
            for ni in xrange(nst_adi):
                for nj in xrange(nst_adi):
                    new_re = math.sqrt( (i*(Hvib2_ave.get(ni,nj).real**2) + (Hvib[tr].get(ni,nj).real**2)/float(num_sh_traj) )/float(i+1) )
                    new_im = math.sqrt( (i*(Hvib2_ave.get(ni,nj).imag**2) + (Hvib[tr].get(ni,nj).imag**2)/float(num_sh_traj) )/float(i+1) )

                    Hvib2_ave.set(ni, nj, new_re, new_im)

             

            #============== TD-SE and surface hopping =================

            # Coherent dynamics
            propagate_electronic(md.dt * md.max_step, Coeff_adi[tr], Hvib[tr])  # propagate in the adiabatic basis

            # Hopping
            ksi = rnd.uniform(0.0, 1.0);  ksi2 = rnd.uniform(0.0, 1.0)
            ist_adi[tr] = tsh.hopping(Coeff_adi[tr], Hvib[tr], ist_adi[tr], sh_method, do_collapse, ksi, ksi2, md.dt *  md.max_step, T)

            denmat_adi_sh[tr] *= 0.0
            denmat_adi_sh[tr].set(ist_adi[tr],ist_adi[tr], 1.0, 0.0)


            f = open("_en_md-%i.txt" % tr,"a")
            f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ST[tr].E_kin, ST[tr].E_pot, ST[tr].E_tot, ST[tr].H_NP, ST[tr].curr_T ))
            f.close()




        #===================== Statistics ==============================

        print "MD step %i took %8.5f seconds " % (i, t.stop())

        tr_adi_adi, tr_adi_dia = 1.0, 1.0
        denmat_adi_se  = tsh.amplitudes2denmat(Coeff_adi)
        ave_pop_adi_sh, ave_pop_adi_se = tsh.ave_pop(denmat_adi_sh, denmat_adi_se)


        #===================== Print out ==============================

        datautils.show_matrix_splot(Hvib2_ave.imag(), "_Hvib_ave_re.dat")
        datautils.show_matrix_splot(Hvib2_ave.real(), "_Hvib_ave_im.dat")


        add_printout(i, ave_pop_adi_sh, tr_adi_adi, tr_adi_dia, "_pop_adi_sh.txt")
        add_printout(i, ave_pop_adi_se, tr_adi_adi, tr_adi_dia, "_pop_adi_se.txt")
        add_printout(i, Hvib[0], -1, -1, "_Hvib.txt")
        add_printout2(i, denmat_adi_se, denmat_adi_sh, Hvib, "_energy-adi.txt")




##==============================================================
## Main simulaiton entry
##==============================================================
##
## SCF
##
"""
rnd = Random()
syst = init_ensembles.init_systems2(1, "CdSe-qd.ent", "true_pdb", rnd, 300.0, 0.0, "elements.dat")

min_ = -1  # HOMO + min_
max_ =  1  # HOMO + max_

run_scf(syst[0], min_, max_)  
"""

##
## SCF
##

act_sp = [0,1]  #
SD_basis = [[0,2], [1,2]]

run_namd(act_sp, SD_basis, 1)


