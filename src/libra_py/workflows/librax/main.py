# *********************************************************************************
# * Copyright (C) 2016-2019 Kosuke Sato, Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 2 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: main
   :platform: Unix, Windows
   :synopsis: This module initialized the calculations and starts the calculations
.. moduleauthor:: Kosuke Sato, Alexey V. Akimov

"""

import md
import extract_qe
import create_input_qe
import x_to_libra_qe
import os
import sys
import math
import copy
import unittest

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.packages.qe.methods as QE_methods

"""
from libra_py import *
from create_input_gms import *
from create_input_qe import *
from create_input_g09 import *
from x_to_libra_gms import *
from x_to_libra_qe import *
from x_to_libra_g09 import *
from md import *
from spin_indx import *
from extract_qe import *
from extract_g09 import *
import include_mm
"""


def construct_active_space(nel, excitations):
    """

    This function constructs the active space for NA-MD calculations

    Args:
        nel ( int ): the number of electrons in our system
        excitations ( list of Libra.excitation objects ): the objects representing the Slater
            Determinats (excitations)

    Returns:
        list of integers: the indices of the orbitals that are used to construct the excitations
            included in the dynamical basis and therefore constitute the active space.
            Returning [-1] means the active space can not be constructed for given number of electrons
            and the list of excitations
    """

    active_space = []
    homo = nel / 2 + nel % 2  # the index of HOMO starting from 1 (not 0!)

    # Find the lowest orbital from which excitation occurs and find the highest orbital to where
    # electron is sent
    min_from = 0
    max_to = 0
    for ex in excitations:

        if min_from > ex.from_orbit[0]:
            min_from = ex.from_orbit[0]
            if min_from + homo <= 0:
                print "Error in construct_active_space: Reconsider the basis of your excitations."
                print "Specifically %5i %5i -> %5i %5i" % (ex.from_orbit[0], ex.from_spin[0], ex.to_orbit[0], ex.from_spin[0])
                print "You are trying to excite from the orbital that is below the lowest possible orbital."
                return [-1]

        if max_to < ex.to_orbit[0]:
            max_to = ex.to_orbit[0]

    # This definition may be a bit excessive, but it is general and correct (at least for single
    # excitations)
    active_space = range(min_from + homo, max_to + homo + 1)

    return active_space


def sanity_check(params):
    """

    This functions runs some sanity check -
    some parameters require specific values of other parameters

    """

    # What happens if you run NVT simulations
    if params["MD_type"] == 1:

        if "therm" in params.keys():
            if params["therm"] is None:
                print "NVT simulations require some valid thermostat. None is given. Exiting..."
                sys.exit(0)
        else:
            print "NVT simulations require thermostat! Use \"therm\" keyword. Exiting..."
            sys.exit(0)


def main(params):
    """

    Finds the keywords and their patterns and extracts the parameters

    Args:
        params ( dictionary ): control parameters

        * **params["dt_nucl"]** ( double ): time step for nuclear dynamics [ units: a.u.; default: 41.0 ]
        * **params["excitations"]** ( list of "excitation" objects ): the dynamical basis
            in the NA-MD calculations
        * **params["excitations_init"]** ( list of ints ): the indices of the initial excitations to
            consider. Index 0 usually corresponds to the ground state

        * **params["tsh_method"]** ( int ): the selector if the NA-MD method



    # \\param[in] params  A list of  input parameters from {gms,qe}_run.py
    #### Returned data:
    #### test_data - the output data for debugging, in the form of dictionary
    #### data - the data extracted from gamess output file, in the form of dictionary
    #
    # This function prepares lists of initial parameters from GAMESS output file
    # and executes run_MD function where NA-MD calculation is done.
    #
    # Used in:  run.py
    """

    # Now try to get parameters from the input
    critical_params = ["excitations"]
    default_params = {"dt_nucl": 41.0, }
    comn.check_input(params, default_params, critical_params)

    dt_nucl = params["dt_nucl"]
    nstates = len(params["excitations"])
    nstates_init = len(params["excitations_init"])
    ninit = params["nconfig"]
    SH_type = params["tsh_method"]
    # nspin = params["nspin"]  This parameter is used only in Libra-QE interface.
    uff = params["ff"]

    kT = params["electronic_smearing"]

    num_SH_traj = 1
    if SH_type >= 1:  # calculate no SH probs.
        num_SH_traj = params["num_SH_traj"]

    ntraj = ninit * nstates_init
    if params["do_rescaling"] == 1:  # nuclear trajectory depends eletronic hops.
        ntraj = ninit * nstates_init * num_SH_traj

    active_space = []
    mo_pool_alp, mo_pool_bet = None, None

    if params["interface"] == "QE":
        # ntraj = nstates*ninit*num_SH_traj
        pass
#        active_space = [5,6,7]  # For C2H4
#    #********** active space is defined here *****************
    elif params["interface"] == "GAMESS":  # or params["interface"]=="G09":
        # ntraj = nstates_init*ninit*num_SH_traj
        for i in range(params["min_shift"], params["max_shift"] + 1):
            active_space.append(i + params["HOMO"] + 1)  # Here MO order start from 1, not 0.
    # *********************************************************

    ################# Step 0: Use the initial file to create a working input file ###############

    #### Step 1: Read initial input, run first calculation, and initialize the "global" variables ####
    # t = Timer()
    # t.start()
    # Initialize variables for a single trajectory first!
    label, Q, R, ao, tot_ene = None, [], None, [], None
    sd_basis = SDList()
    all_grads = []
    e = MATRIX(nstates, nstates)  # contains total energy of excited states!

    if params["interface"] == "GAMESS":

        os.system("cp %s %s" % (params["gms_inp0"], params["gms_inp"]))
        params["gms_inp_templ"] = read_gms_inp_templ(params["gms_inp"])
        exe_gamess(params)
        label, Q, R, grads, E, c, ao, params["nel"] = gms_extract(
            params["gms_out"], params["excitations"], params["min_shift"], active_space, params["debug_gms_unpack"])
        e = MATRIX(E)
        homo = params["nel"] / 2 + params["nel"] % 2

        for ex_st in xrange(nstates):
            mo_pool_alp = CMATRIX(c)
            mo_pool_bet = CMATRIX(c)
            alp, bet = index_spin(params["excitations"][ex_st], active_space, homo)

            # use excitation object to create proper SD object for different excited state
            sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet))
            sd_basis.append(sd)
            all_grads.append(copy.deepcopy(grads))  # newly defined

    elif params["interface"] == "G09":
        os.system("cp %s %s" % (params["g09_inp0"], params["g09_inp"]))
        params["g09_inp_templ"] = read_g09_inp_templ(params["g09_inp"])
        exe_g09(params)
        # while not os.path.exists(params["g09_out"]):
        #    time.sleep(1)
        params["nel"] = g09_extract_first(params["g09_out"], params["debug_g09_unpack"])
        active_space = construct_active_space(params)
        label, Q, R, grads, E, c, ao, params["nel"] = g09_extract(
            params["g09_out"], params["excitations"], params["min_shift"], active_space, params["debug_g09_unpack"])
        e = MATRIX(E)
        homo = params["nel"] / 2 + params["nel"] % 2

        for ex_st in xrange(nstates):
            mo_pool_alp = CMATRIX(c)
            mo_pool_bet = CMATRIX(c)
            alp, bet = index_spin(params["excitations"][ex_st], active_space, homo)
            # alp,bet = index_spin(params["excitations"][0],active_space, homo)
            print "MO_pool", mo_pool_alp.show_matrix()
            # use excitation object to create proper SD object for different excited state
            sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet))
            sd_basis.append(sd)
            all_grads.append(copy.deepcopy(grads))  # newly defined

    elif params["interface"] == "QE":

        info = None

        # AVA: params["qe_inp_templ"] = []
        for ex_st in xrange(nstates):  # modify here #

            print "Starting excitation %i" % (ex_st)

            # Copy the orignal input files into working ones, which will be modified
            # so that we don't screw up the original files if something goes wrong
            os.system("cp x%i.scf.in x%i.scf_wrk.in" % (ex_st, ex_st))
            # AVA: params["qe_inp_templ"].append( read_qe_inp_templ("x%i.scf_wrk.in" % ex_st) )

            print "Starting QE calculations"
            x_to_libra_qe.exe_espresso(ex_st, params)
            status = x_to_libra_qe.check_convergence("x%i.scf.out" % (ex_st))

            print "Status = ", status

            # AVA: flag = 0

            # Basically, here we automatically determine the position of HOMO and will construct
            # the active space before actually using it
            if ex_st == 0:
                info = QE_methods.read_qe_schema("x0.save/data-file-schema.xml")
                # AVA: tot_ene, params["norb"], params["nel"], params["nat"],
                # params["alat"], icoord, iforce = qe_extract_info("x%i.scf.out" % ex_st,
                # ex_st, flag)
                active_space = construct_active_space(info["nelec"], params["excitations"])

            excitation = params["excitations"][ex_st]
            nspin = params["nspin"]
            nel = info["nelec"]
            norbs = info["nbnd"] / nspin  # the number of orbitals per spin channel

            # AVA: occ, occ_alp, occ_bet = excitation_to_qe_occ(params, excitation)
            occ, occ_alp, occ_bet = create_input_qe.excitation_to_qe_occ(norbs, nel, excitation)

            print "Occupation numbers :", occ, occ_alp, occ_bet

            # In the case we didn't converge, lets do the iterations with the occupation numbers
            cnt = 0
            restart_flag = 1  # 1 - change scf_iter, 2 - change scf_iter + restart from the previous pot and wfc

            while status == False:
                print "In the while loop count = %i" % (cnt)

                create_input_qe.write_qe_input_first("x%i.scf_wrk.in" % (ex_st), "x%i.scf_wrk.in" % (ex_st),
                                                     occ, occ_alp, occ_bet, nspin, params["scf_iter"], restart_flag)

                # AVA:  write_qe_input_first("x%i.scf_wrk.in"%ex_st,occ,occ_alp,occ_bet,nspin,params,restart_flag)

                print "starting QE"
                x_to_libra_qe.exe_espresso(ex_st, params)
                status = x_to_libra_qe.check_convergence("x%i.scf.out" % (ex_st))

                print "Status = ", status
                info = QE_methods.read_qe_schema("x%i.save/data-file-schema.xml" % (ex_st), 0)
                print "Energy = ", info["etot"]

                if cnt >= 0:
                    restart_flag = 2

                nbnd = info["nbnd"]

                if nspin == 1:
                    info1, all_e = QE_methods.read_qe_index("x%i.export/index.xml" % (ex_st), range(1, nbnd + 1), 0)

                    bnds = []
                    for i in xrange(nbnd):
                        bnds.append(all_e[0].get(i, i).real)
                    occ = extract_qe.excited_populations(bnds, info["nelec"], nspin, kT, 1)

                    print "Updated occ = ", occ

                elif nspin == 2:
                    nbnd = info["nbnd"]
                    info1, all_e = QE_methods.read_qe_index(
                        "x%i.export/index.xml" %
                        (ex_st), range(
                            1, nbnd / 2 + 1), 0)

                    bnds_alp, bnds_bet = [], []
                    for i in xrange(nbnd / 2):
                        bnds_alp.append(all_e[0].get(i, i).real)
                        bnds_bet.append(all_e[1].get(i, i).real)
                    occ_alp = extract_qe.excited_populations(bnds_alp, info["nelec"], nspin, kT, 1)
                    occ_bet = extract_qe.excited_populations(bnds_bet, info["nelec"], nspin, kT, 1)

                    print "Energies = ", bnds_alp, bnds_bet
                    print "Updated occ = ", occ_alp, occ_bet

                cnt = cnt + 1

    ###########################################
        print "All is done"
        sys.exit(0)

        """   AVA: temporarily comment for debug
            tot_ene, label, R, grads, mo_pool_alp, mo_pool_bet, params["norb"], params["nel"], params["nat"], params["alat"] = qe_extract("x%i.scf.out" % ex_st, active_space, ex_st, nspin, flag)



            homo = params["nel"]/2 +  params["nel"] % 2
            #t.start()
            alp,bet = index_spin(params["excitations"][ex_st], active_space, homo)
            # use excitation object to create proper SD object for different excited state
            sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet) )

            #sd_basis2 = SDList()
            sd_basis.append(sd)
            #t.stop()
            #print "Time to cleade SD object",t.show(),"sec"
            all_grads.append(grads)
            e.set(ex_st, ex_st, tot_ene)

        """

    # Now, clone the single-trajectory variables, to initialize the bunch of such parameters
    # for all trajectories
    sd_basis_list = []
    ham_el_list = []
    ao_list = []
    e_list = []
    grad_list = []
    label_list = []
    Q_list = []
    R_list = []

    for i in xrange(ntraj):
        # AO
        ao_tmp = []
        for x in ao:
            ao_tmp.append(AO(x))
        ao_list.append(ao_tmp)

        # E and C
        e_list.append(MATRIX(e))

        # Slater determinants
        # eventually, the ordering is this: sd_basis_list[traj][ex_st] - type CMATRIX
        sd_basis_list.append(sd_basis)

        # Gradients
        # eventually, the ordering is this: grad_list[traj][ex_st][n_atom] - type VECTOR
        grd2 = []
        for grad in all_grads:  # over all excited states
            grd1 = []
            for g in grad:     # over all atoms
                grd1.append(VECTOR(g))
            grd2.append(grd1)
        grad_list.append(grd2)

        # Coords
        rr = []
        for r in R:
            rr.append(VECTOR(r))
        R_list.append(rr)

        # Labels
        lab = []
        for i in xrange(len(label)):
            lab.append(label[i])
        label_list.append(lab)

        # Q
        qq = []
        for q in Q:
            qq.append(q)
        Q_list.append(qq)

    ################## Step 2: Initialize molecular system and run MD part with TD-SE and SH####

    rnd = Random()  # random number generator object
    syst = []
    syst_mm = []
    el = []

    Ttemp = 0.0
    if params["MD_type"] == 1:  # start NA-MD interacting with thermostat @ t=0
        Ttemp = params["Temperature"]

    # create system object
    Nsh = 1
    if params["do_rescaling"] == 1:
        Nsh = num_SH_traj
    for i in xrange(ninit):
        for i_ex in params["excitations_init"]:
            for itraj in xrange(Nsh):
                df = 0  # debug flag
                # Here we use libra_py module!
                # Utilize the gradients on the ground (0) excited state
                if params["interface"] == "QE":  # For Delta-SCF-NAMD, initialize with excited state (i_ex) gradients.
                    x = init_system.init_system(
                        label_list[i],
                        R_list[i],
                        grad_list[i][i_ex],
                        rnd,
                        Ttemp,
                        params["sigma_pos"],
                        df,
                        "elements.txt")
                else:
                    x = init_system.init_system(
                        label_list[i],
                        R_list[i],
                        grad_list[i][0],
                        rnd,
                        Ttemp,
                        params["sigma_pos"],
                        df,
                        "elements.txt")

                # Add the connectivity - needed if we plan to use MM
                if params["is_MM"]:
                    if os.path.exists(params["ent_file"]):
                        LoadMolecule.Load_Molecule(params["U"], x, params["ent_file"], "pdb")
                    else:
                        print "is_MM is set to 1, which means you need to provide the \
                        a file containing the connectivity information for your MM system. \
                        As of now, such file is called = ", params["ent_file"], " but it cannot \
                        be found on your system. Please check if it exists or set is_MM to 0. \
                        Exiting..."
                        sys.exit(0)

                syst.append(x)

    # all excitations for each nuclear configuration
    for i in xrange(ninit):
        for i_ex in params["excitations_init"]:
            for itraj in xrange(num_SH_traj):
                el.append(Electronic(nstates, i_ex))

    sys.exit(0)  # debug
    # set list of SH state trajectories
    print "run MD"
    md.run_MD(syst, el, ao_list, e_list, sd_basis_list, params, label_list, Q_list, active_space)
    print "MD is done"


class TestConstructActiveSpace(unittest.TestCase):
    def test_active_space(self):
        """Tests the construction of the active space"""

        # For the examples below:  H - HOMO, L - LUMO, a = alpha,b = beta
        #
        #                       GS = H(a)->H(a)
        b1 = [excitation(0, 1, 0, 1)]
        #
        #                        H(a)->H(a)             H(a)->L(a)          H-1(a) -> L(a)
        b2 = [excitation(0, 1, 0, 1), excitation(0, 1, 1, 1), excitation(-1, 1, 1, 1)]
        #
        #                        H(a)->H(a)             H(a)->L+1(a)          H-1(a) -> L(a)
        b3 = [excitation(0, 1, 0, 1), excitation(0, 1, 2, 1), excitation(-1, 1, 1, -1)]

        all_params, exp_res = [], []
        #
        all_params.append({"excitations": b1, "nel": 4})
        exp_res.append([2])
        all_params.append({"excitations": b2, "nel": 4})
        exp_res.append([1, 2, 3])
        all_params.append({"excitations": b3, "nel": 4})
        exp_res.append([1, 2, 3, 4])

        all_params.append({"excitations": b1, "nel": 3})
        exp_res.append([2])
        all_params.append({"excitations": b2, "nel": 3})
        exp_res.append([1, 2, 3])
        all_params.append({"excitations": b3, "nel": 3})
        exp_res.append([1, 2, 3, 4])

        all_params.append({"excitations": b1, "nel": 2})
        exp_res.append([1])
        all_params.append({"excitations": b2, "nel": 2})
        exp_res.append([-1])
        all_params.append({"excitations": b3, "nel": 2})
        exp_res.append([-1])

        for i in range(9):
            print "Running the test with the parameters", all_params[i]
            res = construct_active_space(all_params[i])
            self.assertEqual(res, exp_res[i])


if __name__ == '__main__':
    unittest.main()
