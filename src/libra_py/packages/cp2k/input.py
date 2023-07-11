#***********************************************************
# * Copyright (C) 2023 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import sys
if sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn

def get_kind(elt):
    """
    Args:
        * elt (string): chemical symbol of the atom 

    Returns: 
      dict: the dictionary containing the definition of the key parameters for a given atom type. Namely, the following
         entries are possible:

         * "element" (string): the chemical symbol of the atom of a given type [ e.g. "H"]
         * "basis_set" (string): the type of the basis and the name of the file containing that basis [ e.g. "ORB DZVP-MOLOPT-GTH" ]
         * "potential" (string): the name of the file containing the pseudopotentials [e.g. "GTH-PBE-q4"]
         * "fit_basis_set" (string): the name of the basis set for fitting [e.g. "cFIT3"]
         * "dft_plus_u" (list): if present, DFT+U calculations for this atom kind is requested. The list contains two elements:
           [int, float] - the int is the angular momentum of the orbitals for which the correction is applied, the float is the magnitude of on-site
           correction in atomic units (Ha)

    """
    res = {}

    if elt=="H":
        res = {"element": "H", "basis_set":"ORB DZVP-MOLOPT-GTH", "potential":"GTH-PBE-q4", "fit_basis_set":"cFIT3" }
    elif elt=="Ti":
        res = {"element": "Ti", "basis_set":"ORB DZVP-MOLOPT-SR-GTH", "potential":"GTH-PBE-q12", "fit_basis_set":"cFIT10" }
    return res


def generate(_params):
    """
    Args: 
        _params ( dict ): dictionary of the parameters controlling the calculations

        * **_params["input_filename"]** ( string ) : the name of the input file generated [ default: "md.inp"]
        * **_params["project"]** (string): the project name to be defined in the input file, this variable will
            be used to define the output filenames [ defaults: "Ti17" - this is just for historical reasons]
        * **_params["run_type"]** (string): the run type to be done, see the CP2K manual for the possible options.
            Here, we mainly will use "ENERGY" (for single-point calculations; default) or "MD" (for molecular dynamics)
        * **_params["print_level"]** (string): how much of output to produce by CP2K; see the CP2K manual for the possible options
            [default: "LOW"]. Possible: "MEDIUM", "HIGH"
        * **_params["charge"]** (int): the charge of the system [ default: 0]
        * **_params["multiplicity"]** (int): spin multiplicity of the system: 1 - singlet, 2 - doublet, 3 - triplet, and so on 
            [default: 0]
        * **_params["uks"]** (string): whether to spin-unrestricted calculations (if ".TRUE.") or spin-restrictedo one (if ".FALSE.")
            [default: ".FALSE."]
        * **_params["method"]** (string): the type functional to use. This option affects certain options. Currently possible:
            - "PBE" (for PBE functional)
            - "BEEF" (for BEEF functional), "BEEFWDV"
            - "BP"
            - "LD"
            - "B3LYP"
            - "PBE0"
            - "HSE06"
            - "TPSS"
            - "CAM-B3LYP"
            - "HSE06", "HSE12"
            - "xTB" (for xTB)
        * **_params["max_scf"]** (int): the maximal number of SCF iterations [ default: 100 ]
        * **_params["solver"]** (string): the type of the solver for the SCF procedure. Available options:
            - "OT" (orbital transformation) - can not be use with extra MOs or smearing, only for the ground state calculations
            - "DIAG" (Davidson diagonalization) - the option for excited state calcualations (TD-DFT or sTDA), or for calculations with
              smearing [ default ]
        * **_params["ot.preconditioner"]** (string): the preconditioner for OT calculations, only used if `solver == OT`. 
              Available options: "FULL_ALL", "FULL_SINGLE_INVERSE" [ default ]
        * **_params["ot.minimizer"]** (string): the minimization algorithm for OT calculations, only used if `solver == OT`.
              Available options: "CG" (conjugate gradient), "DIIS" (direct inversion of iterative subspace, default)
        * **_params["ot.linesearch"]** (string): the algorithm for the minimization on an interval in the OT calculations. 
              Only used if `solver == OT`. Available options: "2PNT" [default], "3PNT"
        * **_params["ot.energygap"]** (float): the assumed energy gap for the OT preconditioner in Ha [default: 0.01 Ha] Only used
              if `solver == OT`
        * **_params["diag.preconditioner"]** (string): preconditioner for DIAG SCF solver. Only used inf `solver == DIAG`
              See the CP2K manual for available options. [default: "FULL_ALL"] 
        * **_params["diag.energygap"]** (float): assumed energy gap for the DIAG preconditioner in Ha. Only used if `solver == DIAG`
             [default: 0.01 Ha]
        * **_params["added_mos"]** (int): how many extra MOs (unoccupied) to include in calculations [ default: 20]. Requires 
             `solver == DIAG`, can not be used with "OT" solver.
        * **_params["smearing"]** (Boolean): a flag that enables (if True) fractional occupations of orbitals [default: False]. Requires 
             `solver == DIAG` and `added_mo` to be sufficiently large number to accommodate the populations (need more of the added MOs for
             smaller band gap systems and for larger `smearing_electronic_temperature` parameter). Can not be used with the "OT" solver.
        * **_params["smearing.method"]** (string): the type of smearing to use. See the CP2K manual for options. [defualt: "FERMI_DIRAC"]
        * **_params["smearing.electronic_temperature"]** (float): electronic temperature parameter used in Fermi-Dirac distribution [ default: 300 K]
        * **_params["istate"]** (int): index of the state for which to conduct the calculations [default : 0 - ground state] . For instance, if forces
             are requested in calculations, the forces (and total energy) printed out would refer to this particular state. For any value larger 
             than 0, this keyword turns on the TD-DFT or sTDA calculations. 
        * **_params["nstates"]** (int): how many excited states to include in excited-states calculations with sTDA or TD-DFT [default: 2]
        * **_params["tddft_kernel"]** (string): how to compute the excited states. Only relevant if `istate` > 0. Options include:
            - "FULL" (with DFT methods only, not with xTB) [ default ]
            - "STDA" (with either xTB or DFT)
        * **_params["cell.A"]** (list of 3 floats): the A vector defining the simulation cell
        * **_params["cell.B"]** (list of 3 floats): the B vector defining the simulation cell
        * **_params["cell.C"]** (list of 3 floats): the C vector defining the simulation cell
        * **_params["cell.periodic"]** (string): what kind of periodicity we have. Options: "X", "Y", "Z", "XY", "XZ", "YZ" and "XYZ" [ default : "XYZ"]
        * **_params["xyz_file"]** (string): the name of the xyz file containing the geometry of the system [ default: "inp.xyz"]
        * **_params["kinds"]** (list of dictionaries): definitions of the bases and auxiliary bases for atoms all types
             Each `kind` entry is a dictionary containing the following entries:
             - "element" (string): the chemical symbol of the atom of a given type [ e.g. "H"]
             - "basis_set" (string): the type of the basis and the name of the file containing that basis [ e.g. "ORB DZVP-MOLOPT-GTH" ] 
             - "potential" (string): the name of the file containing the pseudopotentials [e.g. "GTH-PBE-q4"]
             - "fit_basis_set" (string): the name of the basis set for fitting [e.g. "cFIT3"]
             Use the `get_kind` function to facilitate with the construction of such dictionaries 
    """

    h_kind = get_kind("H")
    ti_kind = get_kind("Ti"); ti_kind.update({"dft_plus_u":[2, 0.0]})

    params = dict(_params)
    critical_params = []
    default_params = { "input_filename":"md.inp",
                       "project":"Ti17", "run_type":"ENERGY", "print_level":"LOW",

                       "charge":0, "multiplicity":1, "uks":".FALSE.",

                       "method":"PBE", "max_scf":100,
                       "solver":"DIAG",
                       "ot.preconditioner":"FULL_SINGLE_INVERSE", "ot.minimizer":"DIIS", "ot.linesearch":"2PNT", "ot.energygap":0.01,
                       "diag.preconditioner":"FULL_SINGLE_INVERSE", "diag.energygap":0.01,
                       "added_mos":20, "smearing":False,
                       "smearing.method":"FERMI_DIRAC", "smearing.electronic_temperature":300.0,

                       "istate":0, "nstates":2, "tddft_kernel":"FULL",

                       "cell.A":[30.0, 0.0, 0.0], "cell.B":[0.0, 30.0, 0.0], "cell.C":[0.0, 0.0, 30.0], "cell.periodic":"XYZ",
                       "xyz_file":"inp.xyz", 

                       "kinds": [  h_kind, ti_kind  ]
                     }



#.update({"dft_plus_u": {"L":2, "U_minus_J":0.1 } })

    comn.check_input(params, default_params, critical_params)


    # Unpack
    input_filename = params["input_filename"]
 
    #>============== GLOBAL =================
    project = params["project"]
    run_type = params["run_type"]
    print_level = params["print_level"]

    global_inp = F"""
&GLOBAL
  PROJECT {project}
  RUN_TYPE {run_type}
  PRINT_LEVEL {print_level}
&END GLOBAL
"""

    #>=========== FORCE_EVAL =============

    #>>================== DFT =================
    #>>> ============ XC ===================
    method = params["method"]
    qs_method = "GPW"
    
    xc_input = ""
    xtb_input = ""

    # Shortcut methods
    if method in ["PBE", "BEEF", "BEEFWDV", "BP", "LD", "B3LYP", "PBE0", "TPSS"]:
        xc_input = F"""
    &XC
      &XC_FUNCTIONAL {method}
      &END XC_FUNCTIONAL
    &END XC
        """

    # More advanced onese
    elif method in ["CAM-B3LYP", "HSE06", "HSE12"]:
        func = "HYB_GGA_XC_CAM_B3LYP"
        if method == "CAM-B3LYP":
            func = "HYB_GGA_XC_CAM_B3LYP"
        elif method == "HSE06":
            func = "HYB_GGA_XC_HSE06"
        elif method == "HSE12":
            func = "HYB_GGA_XC_HSE12"
        xc_input = F"""
    &XC
      &XC_FUNCTIONAL
        &{func}
        &END {func}
      &END XC_FUNCTIONAL
      &HF
        &SCREENING
          EPS_SCHWARZ 1.0E-10
        &END
        FRACTION 0.20
      &END HF
    &END XC
        """

    # xTB 
    elif method in ["xTB"]:
        qs_method = "xTB"
        xtb_input = """
      &xTB
        DO_EWALD  F
        CHECK_ATOMIC_CHARGES  F
        COULOMB_INTERACTION T
        &PARAMETER
          DISPERSION_PARAMETER_FILE dftd3.dat
          PARAM_FILE_NAME xTB_parameters
        &END PARAMETER
      &END xTB
        """

    # Explicit names of the functionals
    else:
        xc_input = F"""
    &XC
      &XC_FUNCTIONAL
        &{method}
        &END {method}
      &END XC_FUNCTIONAL
    &END XC
        """     

    #>>> ============ XC END=================

    #>>> ============ AUX ===================
    aux_input = ""
#    if method in ["HSE06", "HSE12" "B3LYP", "CAM-B3LYP", "PBE0", "TPSS"]:
#        aux_input=F"""
#!    &AUXILIARY_DENSITY_MATRIX_METHOD
#!      ! recommended, i.e. use a smaller basis for HFX
#!      ! each kind will need an AUX_FIT_BASIS_SET.
#!      METHOD BASIS_PROJECTION
#!      ! recommended, this method is stable and allows for MD. 
#!      ! can be expensive for large systems
#!      ADMM_PURIFICATION_METHOD MO_DIAG
#!    &END
#        """

    #>>> ============ AUX END ===============

    #>>> ============ QS ====================
    #>>>> ============ SCF ==================

    #>>>>> ============ SOLVER ==============
    solver = params["solver"]

    if solver=="OT":    
        precond = params["ot.preconditioner"]
        minimizer = params["ot.minimizer"]
        lsearch = params["ot.linesearch"]
        egap = params["ot.energygap"]

        solver_input = F"""
      &OT
        PRECONDITIONER {precond}
        MINIMIZER {minimizer}
        ENERGY_GAP {egap}
        LINESEARCH {lsearch}
      &END OT
        """
    elif solver=="DIAG":
        precond = params["diag.preconditioner"]
        egap = params["diag.energygap"]
        solver_input = F"""
      &DIAGONALIZATION
        &DAVIDSON
          PRECONDITIONER {precond}
          ENERGY_GAP {egap}
        &END
      &END DIAGONALIZATION
        """

    #<<<<< ============ SOLVER END ==========

    added_mos = params["added_mos"]
    added_mos_input = ""
    if solver=="DIAG":
        added_mos_input = F"""
      ADDED_MOS {added_mos}
      """    

    smearing = params["smearing"]
    smearing_input = ""
    if solver=="DIAG" and smearing==True:
        sm_meth = params["smearing.method"]        
        sm_temp = params["smearing.electronic_temperature"]
        smearing_input = F"""
      &SMEAR
        METHOD {sm_meth}
        ELECTRONIC_TEMPERATURE {sm_temp}
      &END
      """

    max_scf = params["max_scf"]    

    scf_input = F"""
    &SCF
      MAX_SCF {max_scf}
      SCF_GUESS RESTART
      EPS_SCF 1.e-6        
      {solver_input}
      {added_mos_input}
      {smearing_input}
    &END SCF
    """
    #<<<< ============ SCF END ==================

    qs_input = F"""
    &QS
      METHOD {qs_method}
      {xtb_input}
    &END QS
    """
    #<<<============ QS END ====================

    #>>>============ POISSON ===================
    poisson_input = F"""
    &POISSON
      POISSON_SOLVER PERIODIC
      PERIODIC XYZ
    &END POISSON
    """
    #<<<============ POISSON END ===============

    #>>> ============ EXCITED_STATES ===========
    istate = params["istate"]
    excited_states_input = ""
    if istate > 0:
        excited_states_input = F"""
    &EXCITED_STATES
      STATE {istate}
    &END
        """

    #>>> ============ EXCITED_STATES ===========


    #<<<============ PRINT =====================
    print_input = F"""
    &PRINT
      &MO
        ENERGIES .TRUE.
        OCCUPATION_NUMBERS .TRUE.
        NDIGITS 8
        &EACH
          QS_SCF 0
        &END
      &END
    &END
    """   
    #<<<============ PRINT END =================

    charge = params["charge"]
    multiplicity = params["multiplicity"]
    uks = params["uks"]

    dft_input = F"""
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    BASIS_SET_FILE_NAME BASIS_ADMM
    BASIS_SET_FILE_NAME BASIS_ADMM_MOLOPT
    BASIS_SET_FILE_NAME GTH_BASIS_SETS
    POTENTIAL_FILE_NAME GTH_POTENTIALS

    CHARGE {charge}
    MULTIPLICITY {multiplicity}
    UKS {uks}

    {xc_input}
    {aux_input}
    {scf_input}
    {qs_input}
    {poisson_input}
    {excited_states_input}
    {print_input}
  &END DFT

    """
    #<<================== DFT END ==============


    #>>============ PROPERTIES =================
    istate = params["istate"]
    properties_input = ""

    if istate > 0:  # REQUEST ACTUAL EXCITED STATES CALCULATIONS 
        nstates = params["nstates"]
        kernel = params["tddft_kernel"]
        if method=="xTB":
            kernel = "STDA"

        properties_input = F"""
  &PROPERTIES
    &TDDFPT
      KERNEL {kernel}
      NSTATES     {nstates}     # number of excited states
      MAX_ITER    200           # maximum number of Davidson iterations
      CONVERGENCE [eV] 1.0e-5   # convergence on maximum energy change between iterations

      &MGRID
        NGRIDS 4
        CUTOFF 500 # separate cutoff for TDDFPT calc
      &END
!        Only in case you have a tdwfn file from previous calculations
!       RESTART     .TRUE.
!       WFN_RESTART_FILE_NAME RESTART.tdwfn
    &END TDDFPT
  &END PROPERTIES
    """
    #<<============ PROPERTIES END =============


    #>>============ SUBSYS =====================
    Ax, Ay, Az = params["cell.A"][0], params["cell.A"][1], params["cell.A"][2]
    Bx, By, Bz = params["cell.B"][0], params["cell.B"][1], params["cell.B"][2]
    Cx, Cy, Cz = params["cell.C"][0], params["cell.C"][1], params["cell.C"][2]
    xyz = params["cell.periodic"]

    cell_input = F"""
    &CELL
      A {Ax} {Ay} {Az}
      B {Bx} {By} {Bz}
      C {Cx} {Cy} {Cz}
      PERIODIC {xyz}
    &END CELL
    """  

    xyz_file = params["xyz_file"]
    topo_input = F"""
    &TOPOLOGY   
      COORD_FILE_NAME {xyz_file}
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES T
      &END
    &END    
    """

    kinds = params["kinds"]
    print(kinds)
    kind_input = ""

    for ikind in kinds:
        elt = ikind["element"]
        basis = ikind["basis_set"]
        fit_basis = ikind["fit_basis_set"]
        pot = ikind["potential"]
        dft_plus_u_input = ""
        if "dft_plus_u" in ikind.keys() and method is not "xTB":
            l_val = ikind["dft_plus_u"][0]
            u_val = ikind["dft_plus_u"][1]
            dft_plus_u_input = F"""
      &DFT_PLUS_U
        L {l_val}
        U_MINUS_J {u_val}
      &END DFT_PLUS_U   
            """
        kind_input = F"""{kind_input}
    &KIND {elt}
      ELEMENT {elt}
      BASIS_SET {basis}
      BASIS_SET AUX_FIT {fit_basis} 
      POTENTIAL {pot}
      {dft_plus_u_input}
    &END
        """

    subsys_input = F"""
  &SUBSYS
    {cell_input}
    {topo_input}
    {kind_input}
  &END SUBSYS
    """
    #<<============ SUBSYS END =================

    #<============= FORCE_EVAL END ================

    force_eval_input = F"""
&FORCE_EVAL
  {dft_input}
  {properties_input}
  {subsys_input}
&END FORCE_EVAL
"""


    #============ OVERALL ===============  
    input = F"""{global_inp} {force_eval_input}"""


    f = open(input_filename, "w+")
    f.write(input)
    f.close()



""" Examples of usage:

generate({"input_filename":"default1.inp", "solver":"OT"})
generate({"input_filename":"default2.inp", "solver":"DIAG"})
generate({"input_filename":"default3.inp", "solver":"OT", "smearing":True, "added_mos":50 })
generate({"input_filename":"default4.inp", "solver":"DIAG", "smearing":True, "added_mos":50 })                                                                             
generate({"input_filename":"default5.inp", "solver":"DIAG", "smearing":True, "added_mos":50, "istate":1 })
generate({"input_filename":"default6.inp", "solver":"DIAG", "method":"CAM-B3LYP" })
generate({"input_filename":"default7.inp", "solver":"DIAG", "method":"HSE06" })

generate({"input_filename":"xtb1.inp", "method":"xTB", "solver":"OT"})
generate({"input_filename":"xtb2.inp", "method":"xTB", "solver":"DIAG"})
generate({"input_filename":"xtb3.inp", "method":"xTB", "solver":"DIAG", "istate":1 })

"""
