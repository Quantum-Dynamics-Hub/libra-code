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
  res = {}

  if elt=="H":
      res = {"element": "H", "basis_set":"ORB DZVP-MOLOPT-GTH", "potential":"GTH-PBE-q4" }
  elif elt=="Ti":
      res = {"element": "Ti", "basis_set":"ORB DZVP-MOLOPT-SR-GTH", "potential":"GTH-PBE-q12" }
  return res


def generate(_params):
    """
    method : PBE, BEEF, xTB

    solver:  OT, DIAG

    ot.preconditioner: FULL_ALL, FULL_SINGLE_INVERSE
    ot.minimizer: "CG", "DIIS"
    ot.linesearch: 2PNT, 3PNT
    ot.energygap: 0.001

    diag.preconditioner: FULL_ALL
    diag.energygap: 0.001

    tddft_kernel: FULL (no xTB), STDA (xTB and DFT)
    """

    h_kind = get_kind("H")
    ti_kind = get_kind("Ti"); ti_kind.update({"dft_plus_u":[2, 0.0]})

    params = dict(_params)
    critical_params = []
    default_params = { "input_filename":"md.inp",
                       "project":"Ti17", "run_type":"ENERGY", "print_level":"LOW",

                       "charge":0, "multiplicity":1, "uks":".FALSE.",

                       "method":"PBE", "max_scf":100,
                       "solver":"diag",
                       "ot.preconditioner":"FULL_SINGLE_INVERSE", "ot.minimizer":"DIIS", "ot.linesearch":"2PNT", "ot.energygap":0.001,
                       "diag.preconditioner":"FULL_ALL", "diag.energygap":0.001,
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
    if method in ["PBE", "BEEF"]:
        xc_input = F"""
    &XC
      &XC_FUNCTIONAL {method}
    &END XC_FUNCTIONAL
        """
    elif method in ["HSE06"]:
        xc_input = F"""
    &XC
      &XC_FUNCTIONAL
        &XWPBE
          SCALE_X -0.25
          SCALE_X0 1.0 
          OMEGA 0.11
        &END XWPBE
        &PBE
          SCALE_X 0.0
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        &SCREENING
          EPS_SCHWARZ 1.0E-6
          EPS_SCHWARZ_FORCES 1.0E-5
          SCREEN_ON_INITIAL_P FALSE
        &END SCREENING
        &INTERACTION_POTENTIAL
          CUTOFF_RADIUS 10
          POTENTIAL_TYPE SHORTRANGE
          OMEGA 0.11
        &END INTERACTION_POTENTIAL
        FRACTION 0.25
      &END HF
    &END XC
        """
    elif method in ["B3LYP", "CAM-B3LYP"]:
        func = "XC_HYB_GGA_XC_B3LYP"
        if method == "CAM-B3LYP":
            func = "XC_HYB_GGA_XC_CAM_B3LYP"
        xc_input = F"""
    &XC
      &XC_FUNCTIONAL
        &LIBXC
          FUNCTIONAL {func}
        &END LIBXC
      &END XC_FUNCTIONAL
      &HF
        &SCREENING
          EPS_SCHWARZ 1.0E-10
        &END
        FRACTION 0.20
      &END HF
    &END XC
        """
 
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
    #>>> ============ XC END=================

    #>>> ============ AUX ===================
    aux_input = ""
    if method in ["HSE06", "B3LYP", "CAM-B3LYP"]:
        aux_input=F"""
    &AUXILIARY_DENSITY_MATRIX_METHOD
      ! recommended, i.e. use a smaller basis for HFX
      ! each kind will need an AUX_FIT_BASIS_SET.
      METHOD BASIS_PROJECTION
      ! recommended, this method is stable and allows for MD. 
      ! can be expensive for large systems
      ADMM_PURIFICATION_METHOD MO_DIAG
    &END
        """

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

      &POISSON
        POISSON_SOLVER PERIODIC
        PERIODIC XYZ      
      &END POISSON  
      {scf_input}
    &END QS
    """
    #<<<============ QS END ====================

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


    dft_input = F"""
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    BASIS_SET_FILE_NAME BASIS_ADMM
    BASIS_SET_FILE_NAME BASIS_ADMM_MOLOPT
    BASIS_SET_FILE_NAME GTH_BASIS_SETS
    POTENTIAL_FILE_NAME GTH_POTENTIALS

    {xc_input}
    {aux_input}
    {qs_input}
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
      BASIS_SET AUX_FIT cFIT6  ! this value is taken from the Mohammad's example for Pb
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
    charge = params["charge"]
    multiplicity = params["multiplicity"]
    uks = params["uks"]

    force_eval_input = F"""
&FORCE_EVAL
  CHARGE {charge}
  MULTIPLICITY {multiplicity}
  UKS {uks}
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
