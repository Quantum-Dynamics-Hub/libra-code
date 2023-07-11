#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  This code demonstrates simple NA-MD simulation in a molecular system
#  Approximations: NBRA
#  Setups: nuclear dynamics: UFF, electronic dynamics: FSSH/MSSH
#          electronic structure: EHT
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





def compute_Hvib(ham_old, ham_cur, orb_adi, dt):
    ##         
    # Computes the vibronic Hamiltonian
    # \param[in] ham_old (Hamiltonian object) - unperturbed at t-dt
    # \param[in] ham_cur (Hamiltonian object) - unperturbed at t
    # \param[in] orb_dia (list of ints) - indices of the orbitals included in the diabatic active space
    # \param[in] orb_adi (list of ints) - indices of the orbitals included in the adiabatic active space
    # \param[in] dt (float) Time step (in a.u.)
    # \param[in] orb_prime (list of ints) - indiced of the orbitals coupled by the perturbation
    # \param[in] H_prime (float) - the magnitude of the perturbation

    N_adi = len(orb_adi)

    es_old = ham_old.get_electronic_structure()
    Nocc = es_old.Nocc_alp
    Fao_old = CMATRIX(es_old.get_Fao_alp() )
    Sao_old = CMATRIX(es_old.get_Sao() )
    res_old = Fock_to_P(Fao_old, Sao_old, Nocc, 1.0, 0.001, 1e-8, 0) 
    E_old = res_old[0] # MO energies
    C_old = res_old[1] # MO-LCAO
    e_old = CMATRIX(N_adi, N_adi)
    pop_submatrix(E_old, e_old, orb_adi, orb_adi)


    es_cur = ham_cur.get_electronic_structure()
    Nocc = es_cur.Nocc_alp
    Fao_cur = CMATRIX(es_cur.get_Fao_alp() )
    Sao_cur = CMATRIX(es_cur.get_Sao() )
    res_cur = Fock_to_P(Fao_cur, Sao_cur, Nocc, 1.0, 0.001, 1e-8, 0) 
    E_cur = res_cur[0] # MO energies
    C_cur = res_cur[1] # MO-LCAO
    e_cur = CMATRIX(N_adi, N_adi)
    pop_submatrix(E_cur, e_cur, orb_adi, orb_adi)



    # <MO(t)|MO(t-dt)>
    dist2 = 10000.0

    D = CMATRIX(N_adi,N_adi)
    MO_overlap(D, ham_old.basis_ao, ham_cur.basis_ao, C_old, C_cur, Py2Cpp_int(orb_adi), Py2Cpp_int(orb_adi), dist2 )  # 
    D = (0.5/dt)*(D - D.H())



    Eadi = 0.5*(e_old + e_cur)        
    Hvib = Eadi - 1.0j*D       # direct adiabatic

    datautils.show_matrix_splot(Hvib.imag(), "_Hvib_im.dat")
    datautils.show_matrix_splot(Hvib.real(), "_Hvib_re.dat")
    
    return Hvib




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




def print_mo_frag(active_orb, sub_ham_cur, syst, prefix):

    nfrags = len(active_orb) 
    # For fragments
    for fr in xrange(nfrags):
        prms = Control_Parameters()
        prms.nx_grid, prms.ny_grid, prms.nz_grid = 40, 40, 40
        prms.charge_density_prefix = prefix+"_fmo_"
        prms.orbs = Py2Cpp_int(active_orb[fr]) 
        el_str_fr = sub_ham_cur[fr].get_electronic_structure()
        charge_density( el_str_fr, syst, sub_ham_cur[fr].basis_ao, prms)



def print_mo_adi(active_orb_adi, ham_cur, syst, prefix):

    prms = Control_Parameters()
    prms.nx_grid, prms.ny_grid, prms.nz_grid = 40, 40, 40
    prms.charge_density_prefix = prefix+"_mo_"
    prms.orbs = Py2Cpp_int(active_orb_adi) 
    el_str_adi = ham_cur.get_electronic_structure()
    charge_density( el_str_adi, syst, ham_cur.basis_ao, prms)





        

def main():

    #--------------------- Parameters ----------------------

    T = 278.0  # K
    sigma = 0.1  # 

    do_collapse = 0  # 0 - no decoherence, 1 - decoherence
    sh_method = 1    # 0 - MSSH,  1 - FSSH
    num_sh_traj = 5
    H_prime = 0.01      # strength of the perturbation

    homo = 29 # need to know this in advance

    # Define diabatic and adiabatic active spaces   
    orb_adi = range(homo-10,homo+10)
    init_st_adi = 11                      

    nst_adi = len(orb_adi)


        
    #--------------------- Initialization ----------------------

    rnd = Random()

    # Create molecular system and initialize the properties
    syst = init_ensembles.init_systems2(num_sh_traj, "2benz_aa.ent", "pdb", rnd, T, sigma, "elements.dat")

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
#    control_filename = "control_parameters_indo.dat"


    dt = 20.0

    ham_cur, ham_old, Hvib = [], [], []

    for tr in xrange(num_sh_traj):
        ham_cur.append( listHamiltonian_QM(control_filename, syst[tr] ) ) # this includes the overlap calculation
        ham_cur[tr].compute_scf(syst[tr])                  # this includes the core Hamiltonian calculations

        ham_old.append( listHamiltonian_QM(ham_cur[tr]) )

        es = ham_cur[tr].get_electronic_structure() 
        if (es.Nocc_alp - 1 != homo):
            print "Error: Re-set HOMO to the proper index = ", es.Nocc_alp - 1
            sys.exit(0)

        homo = es.Nocc_alp - 1  # index of the HOMO orbital in the adiabatic sub-system
        print "Adiabatic sub-system: number of electrons (alpha) = ", es.Nocc_alp, " homo = ", homo       


        Hvib.append( compute_Hvib(ham_old[tr], ham_cur[tr], orb_adi, dt) )


        if tr==num_sh_traj-1:        
            print "Hvib = ";  Hvib[tr].show_matrix();        
            print "Electronic structure of the adiabatic sub-system"
            print "Index  Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
            for j in xrange(es.Norb):
                print "%5i  %12.8f   %12.8f  %12.8f   %12.8f" %(j, es.get_bands_alp(j), es.get_occ_alp(j), es.get_bands_bet(j), es.get_occ_bet(j) )
    


    #=================== Propagation ====================
    ########################## Cooling #################################

    md = MD({"max_step":100,"ensemble":"NVE","integrator":"DLML","terec_exp_size":10,"dt":40.0,"n_medium":1,"n_fast":1,"n_outer":1}); 
    md.show_info()


    # State: system + thermostat + (optional barostat)
#    anneal_schedule = [ {"dt":20.0, "nsteps":5, "ncycles":20}] # for test (very fast)
#    anneal_schedule = [ {"dt":20.0, "nsteps":50, "ncycles":20}] # for test
    anneal_schedule = [ {"dt":20.0, "nsteps":500, "ncycles":20}]


    ST = []
    for tr in xrange(num_sh_traj):
        ST.append( State() )
        ST[tr].set_system(syst[tr]); ST[tr].set_thermostat(therm[tr]); ST[tr].set_md(md)
        ST[tr].init_md(mol[tr], el[tr], ham_mm[tr], rnd)

        f = open("_en_cooling-%i.txt" % tr,"w"); f.close()

        for item in anneal_schedule:
            md.dt = item["dt"]; md.max_step = item["nsteps"] 
        
            for i in xrange(item["ncycles"]):
                syst[tr].set_atomic_q(mol[tr].q)
                syst[tr].print_xyz("_mol_cooling-%i.xyz" % tr,i)
                ST[tr].run_md(mol[tr], el[tr], ham_mm[tr])
                ekin = ST[tr].E_kin; epot = ST[tr].E_pot
                ST[tr].cool()
        
                f = open("_en_cooling-%i.txt" % tr,"a")
                f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ekin, epot, ST[tr].E_tot, ST[tr].H_NP, ST[tr].curr_T ))
                f.close()



    ########################## Production MD - thermalization #################################

    for tr in xrange(num_sh_traj):
        syst[tr].init_atom_velocities(T, rnd)  # must be this!!!

        md.max_step = 10;  md.ensemble = "NVT"; md.dt = 20.0;
        f = open("_en_md-therm-%i.txt" % tr,"w"); f.close()

        for i in xrange(50):  # for test
            syst[tr].set_atomic_q(mol[tr].q)
            syst[tr].print_xyz("_mol_md-therm-%i.xyz" % tr,i)
            ST[tr].run_md(mol[tr], el[tr], ham_mm[tr])

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
        md.dt = 20.0
  

        ist_adi.append(init_st_adi)
        Coeff_adi.append(CMATRIX(nst_adi, 1)); Coeff_adi[tr].set(init_st_adi, 1.0, 0.0)
        denmat_adi_se  = tsh.amplitudes2denmat(Coeff_adi)
        denmat_adi_sh.append(CMATRIX(nst_adi, nst_adi)); denmat_adi_sh[tr].set(init_st_adi, init_st_adi, 1.0, 0.0)


    f = open("_pop_adi_sh.txt","w");    f.close()
    f = open("_pop_adi_se.txt","w");    f.close()
    f = open("_Hvib.txt","w");          f.close()
    f = open("_energy-adi.txt","w");    f.close()



    #=======================================================================    

    for i in xrange(1000):        

        t = Timer();  t.start()
        for tr in xrange(num_sh_traj):

            syst[tr].set_atomic_q(mol[tr].q)
            syst[tr].print_xyz("_mol_md-%i.xyz" % tr, i)
            ST[tr].run_md(mol[tr], el[tr], ham_mm[tr])

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

            Hvib[tr] = compute_Hvib(ham_old[tr], ham_cur[tr], orb_adi, dt)


            #============== TD-SE and surface hopping =================

            # Coherent dynamics
            propagate_electronic(md.dt, Coeff_adi[tr], Hvib[tr])  # propagate in the adiabatic basis

            # Hopping
            ksi = rnd.uniform(0.0, 1.0);  ksi2 = rnd.uniform(0.0, 1.0)
            ist_adi[tr] = tsh.hopping(Coeff_adi[tr], Hvib[tr], ist_adi[tr], sh_method, do_collapse, ksi, ksi2, md.dt, T)

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


        add_printout(i, ave_pop_adi_sh, tr_adi_adi, tr_adi_dia, "_pop_adi_sh.txt")
        add_printout(i, ave_pop_adi_se, tr_adi_adi, tr_adi_dia, "_pop_adi_se.txt")
        add_printout(i, Hvib[0], -1, -1, "_Hvib.txt")
        add_printout2(i, denmat_adi_se, denmat_adi_sh, Hvib, "_energy-adi.txt")



main()
