from liblibra_core import *

def load(dyn_general):

    ############################ HOW TO COMPUTE PROPERTIES #########################
    #====== How to update Hamiltonian ===================
    # Options:
    # - 0: don't update any Hamiltonians
    # - 1: recompute only diabatic Hamiltonian, common choice for model Hamiltonians [ default ]
    # - 2: recompute only adiabatic Hamiltonian, use with file-based or on-the-fly workflows
    dyn_general.update({"ham_update_method":1}) 


    #====== How to transform the Hamiltonians between representations ============
    # Options:
    # - 0: don't do any transforms;  usually for NBRA or on-the-fly workflows, so you don't override the read values
    # - 1: diabatic->adiabatic according to internal diagonalization [ default ]
    # - 2: diabatic->adiabatic according to internally stored basis transformation matrix
    # - 3: adiabatic->diabatic according to internally stored basis transformation matrix
    # - 4: adiabatic->diabatic according to local diabatization method
    dyn_general.update( {"ham_transform_method":1 }) 


    #====== How do get the time-overlaps in the dynamics ========
    # Options:
    #  - 0: don't compute it (perhaps because it was already pre-computed or read)
    #  - 1: explicitly compute it from the wavefunctions (the Hamiltonian shall have the basis_transform variables updated)  [ default ]
    dyn_general.update( {"time_overlap_method":1 }) 


    #================== How to update NACs ===============================
    # Options:
    #  - 0: don't update them (e.g. for simplest NAC)
    #  - 1: update according to changed momentum and existing derivative couplings [ default ]
    #  - 2: update according to time-overlaps (only time-derivative NACs)
    dyn_general.update({"nac_update_method":1}) 

    #================ How to compute time-derivative NACs ================
    # Options:
    #  - [-1]: don't, e.g. we use NACs from somewhere else [ default ]
    #  -  0: use HST formula (if nac_update_method==2)
    #  -  1: use NPI of Meek and Levine (if nac_update_method==2)
    dyn_general.update({"nac_algo":-1})    


    #============== How to update vibronic Hamiltonian ==============
    # Options:
    #  - 0: don't update them (e.g. if it is read externally, maybe because we read it from files)
    #  - 1: update according to regular formula: Hvib = Ham - i * hbar * NAC [ default ]
    dyn_general.update( {"hvib_update_method":1 })


    #=========== Phase correction of SSY =================
    # Whether to modify the Hamiltonian in the dynamics according the Shenvi-Subotnik-Yang (SSY)
    # method, see my Chapter Eq. 3.27
    # Note, that this is only applied in the adiabatic representation
    #
    # Options:
    #  - 0: don't [ default ]
    #  - 1: do
    dyn_general.update({"do_ssy":0 }) 


    ############################# FORCES #################################

    #======== In which representation to compute forces ==================
    # To clarify - the forces
    # in both representations shall be equivalent, so this flag actually
    # selects the type of the properties needed to compute such forces.
    # For instance, if it is set to 0, we may be using the derivatives
    # of the diabatic Hamiltonians, diabatic states overlaps, etc.
    #
    # Options:
    #  - 0: diabatic
    #  - 1: adiabatic [ default ]
    dyn_general.update( {"rep_force":1} ) 

    #================ How to compute forces in the dynamics =============
    # Options:
    #  - 0: don't compute forces at all - e.g. we do not really need them
    #  - 1: state-specific  as in the TSH or adiabatic (including adiabatic excited states) [ default ]
    #  - 2: Ehrenfest and MMST (SQC)
    #  - 3: QTSH force
    #  - 4: KC-RPMD force
    dyn_general.update( {"force_method":1})

    #============= Nonclassical force options in the adiabatic QTSH =========
    # Only used with `use_qtsh == 1`
    #
    #  Options:
    #    - 0: Considering only the first-order force, i.e., off-diagonal Ehrenfest force
    #    - 1: The whole force including the second-order term is used [default]
    dyn_general.update({"qtsh_force_option": 0})


    #========= Whether to use the decoherence force in MQCXF ====================
    # The corresponding electronic propagation 
    # is adjusted for the energy conservation. Only used with `decoherence_algo == 6`
    #
    # Options:
    #  - 0: don't use it, so for XF-methods this is only Ehrenfest-like force; EhXF [ default ]
    #  - 1: The whole force including the XF-correction; MQCXF 
    dyn_general.update({"use_xf_force": 0})


    ############################ SURFACE HOPPING #########################

    #============ Surface hop proposal options =================
    #Surface hop proposal methodology.
    #
    #Options:
    #  - [-1]: adiabatic dynamics, no hops [ default ]
    #  - 0: Fewest Switches Surface Hopping (FSSH)
    #  - 1: Global Flux Surface Hopping (GFSH)
    #  - 2: Markov-State Surface Hopping (MSSH)
    #  - 3: Landau-Zener (LZ) options
    #  - 4: Zhu-Nakamura (ZN) options
    #  - 5: DISH
    #  - 6: MASH
    #  - 7: FSSH2
    #  - 8: FSSH3
    #  - 9: GFSH (original)
    dyn_general.update({"tsh_method":9 })


    #Whether to use QTSH - this replaces standard FSSH via the QTSH
    #
    #Options:
    #
    #  - 0: don't apply [ default ]
    #  - 1: use it
    dyn_general.update({"use_qtsh":0 })


    #============ Surface hop acceptance criteria =================
    #Options to control the acceptance of the proposed hops.
    #Options:
    #  - 0: accept all proposed hops  [ default ]
    #
    #  - 10: based on adiabatic energy - accept only those hops that can obey the energy conservation with 
    #        adiabatic potential energies
    #  - 11: based on diabatic energy - same as 10, but we use diabatic potential energies
    #
    #  - 20: based on derivative coupling vectors - accept only those hops that can obey the energy conservation
    #        by rescaling nuclear velocities along the directions of derivative couplings for the quantum nuclear DOFs                   
    #  - 21: based on difference of state-specific forces - same as 20, but the rescaling is done along the vector
    #        parallel to the difference of adiabatic forces on initial and target states
    #
    #  - 31: accept hops with the probability taken from the quantum Boltzmann distribution
    #  - 32: accept hops with the probability taken from the classical Maxwell-Boltzmann distribution
    #  - 33: accept hops with the probability taken from the updated quantum Boltzmann distribution (experimental)
    #
    #  - 40: based on possibility to conserve energy using tcnbra_ekin variables (for TC-NBRA)
    dyn_general.update({"hop_acceptance_algo":10 })


    #=============== Momentum rescaling options =======================
    #Options to control nuclear momenta changes upon successful or frustrated hops.
    #
    #Options:
    #
    #  - 0: don't rescale [ default ]
    #
    #  - 100: based on adiabatic energy, don't reverse on frustrated hops
    #  - 101: based on adiabatic energy, reverse on frustrated hops
    #  - 110: based on diabatic energy, don't reverse on frustrated hops
    #  - 111: based on diabatic energy, reverse on frustrated hops
    #
    #  - 200: along derivative coupling vectors, don't reverse on frustrated hops
    #  - 201: along derivative coupling vectors, reverse on frustrated hops
    #  - 210: along difference of state-specific forces, don't reverse on frustrated hops
    #  - 211: along difference of state-specific forces, reverse on frustrated hops
    #
    #  - 40: does not rescale velocities, but rescales  tcnbra_ekin variables
    dyn_general.update({"momenta_rescaling_algo":0 })  # accept and rescale based on force differences, reverse on frustrated


    #=============== Jasper-Truhlar criterion for momentum reversal ==================
    #A flag to turn on/off the Jasper-Truhlar criterion for reversal of velocities on frustrated hops.
    #According to: Jasper, A. W.; Truhlar, D. G. Chem. Phys. Lett. 2003, 369, 60− 67
    #
    #Options:
    #
    #  - 0: don't use this criterion (naive handling)
    #  - 1: use it [ default ] - the velocities are reversed along the direction d_{a,j} if
    #    a) (F_a * d_{a,j}) * (F_j * d_{a,j}) < 0 and b) (v * d_{a,j}) * (F_j * d_{a,j}) < 0 
    #    where a - is the active state index;  Only in effect, if `momenta_rescaling_algo == 201`
    dyn_general.update({"use_Jasper_Truhlar_criterion":1 })


    ############################ DECOHERENCE #########################

    #=========== Decoherence options =================
    #Selector of the method to incorporate decoherence.
    # 
    #Options:
    #  - [-1]: no decoherence [ default ]
    #  - 0: SDM and alike
    #  - 1: instantaneous decoherence options (ID-S, ID-A, ID-C)
    #  - 2: AFSSH
    #  - 3: BCSH of Linjun Wang
    #  - 4: MF-SD of Bedard-Hearn, Larsen, Schwartz
    #  - 5: SHXF of Min
    #  - 6: MQCXF
    #  - 7: DISH, rev2023
    #  - 8: diabatic IDA, experimental
    #  - 9: simple decoherence, experimental
    dyn_general.update({ "decoherence_algo":7}) 


    #==== Option to control the instantaneous decoherence methodology ========
    #only used with decoherence_algo == 1
    #
    #  - 0: ID-S
    #  - 1: ID-A [default] - if the proposed hop is not successful, we project back to the initial state
    #                        if the proposed hop is accepted - we project onto that state
    #  - 2: ID-C - consistent ID - an experimental algorithm
    #  - 3: ID-A, new: if the proposed hop is not successful, we project out the proposed states
    #                  if the proposed hop is accepted - we project onto that state
    #  - 4: ID-F, new: if the proposed hop is not successful, we project out the proposed states
    #                  but we don't do anything if the hop is successful
    dyn_general.update({ "instantaneous_decoherence_variant":1})


    #=========== Decoherence times options ==================
    #Type of dephasing times/rates calculation:
    #
    #  - [-1]: set all dephasing rates to zero, infinite decoherence times [ default ]
    #  - 0: use the rates read out from the input 
    #  - 1: use the energy-based decoherence method (EDC)    
    #  - 2: Schwartz - mean-field Force-based decoherence (Schwartz 1), using inv_alpha
    #  - 3: Schwartz - pair-wise-based decoherence, (Schwartz 2), using inv_alpha
    #  - 4: Schwartz - mean-field Force-based decoherence (Schwartz 1), but using interaction width
    #  - 5: Gu-Franco 
    dyn_general.update({"decoherence_times_type":3 })

    #================ Decoherence time parameters ===================    
    # For "decoherence_times_type":1
    dyn_general.update( { "decoherence_C_param": 1.0, "decoherence_eps_param":0.1 } )
   
    # For "decoherence_times_type":2 and "decoherence_times_type":3, "decoherence_times_type":4
    #MATRIX(ndof, 1) of 1/alpha - the parameters used in GWP in
    #computing decoherence rates [ default: NULL ]
    A = MATRIX(1,1); A.set(0, 0, 1.0); # dimension and values are arbitrary!!!
    dyn_general.update( { "schwartz_decoherence_inv_alpha":A } )

    #MATRIX(ndof, 1) - the parameters for the spatial extent of NAC in
    #computing decoherence rates [ default: NULL ]
    B = MATRIX(1,1); B.set(0,0, 1.0)
    dyn_general.update( { "schwartz_interaction_width":B } )

    #Reorganization energy of the bath, Ha  - for "decoherence_times_type":5
    #default: 0.0 Ha
    dyn_general.update( { "reorg_eergy":0.0 } )

    #A flag to apply the dephasing-informed approach of Sifain et al 
    #to correct dephasing times: 
    # 
    #  - 0: don't apply [ default ]
    #  - 1: use it 
    dyn_general.update( {"dephasing_informed":0})

    # Other options
    dyn_general.update({"decoherence_rates":MATRIX(2,2), "ave_gaps":MATRIX(2,2) } )
                               


    ############################ INTEGRATORS, TRACKING, PHASES  #########################
    #the method for electronic TD-SE integration:
    #
    #rep_tdse = 0 (diabatic): 1** - with NBRA
    #
    #    -1              - No propagation
    #
    #     0              - Lowdin exp_ with 2-point Hvib_dia 
    #     1              - based on QTAG propagator
    #     2              - based on modified QTAG propagator (Z at two times)
    #     3              - non-Hermitian integrator with 2-point Hvib_dia
    #
    #rep_tdse = 1 (adiabatic):  1** - with NBRA
    #
    #    -1              -  No propagation
    #
    #     0              -  ld, with crude splitting,  with exp_  [ default ]
    #     1              -  ld, with symmetric splitting, with exp_
    #     2              -  ld, original, with exp_
    #     3              -  1-point, Hvib integration, with exp_
    #     4              -  2-points, Hvib integration, with exp_
    #     5              -  3-points, Hvib, integration with the second-point correction of Hvib, with exp_
    #     6              -  same as 4, but without projection matrices (T_new = I)
    #
    #    10              -  same as 0, but with rotations
    #    11              -  same as 1, but with rotations
    #    12              -  same as 2, but with rotations
    #    13              -  same as 3, but with rotations
    #    14              -  same as 4, but with rotations
    #    15              -  same as 5, but with rotations
    #
    #rep_tdse = 2 ( diabatic, density matrix formalism): 1** - with NBRA
    #
    #     0              -  mid-point Hvib with the second-point correction of Hvib
    #
    #rep_tdse = 3 ( adiabatic, density matrix formalism): 1** - with NBRA
    #
    #     0              -  mid-point Hvib with the second-point correction of Hvib
    #     1              -  Zhu Liouvillian
    #
    #    10              -  same as 0, but with rotations
    dyn_general.update({"rep_tdse":1, "electronic_integrator":2 })   # ld, original, with exp_

    #=========== Disable state tracking and phase corrections explicitly for the LD integrators ===============
    #State tracking algorithm:
    #  - -1: use LD approach, it includes phase correction too [ default ]
    #  - 0: no state tracking
    #  - 1: method of Kosuke Sato (may fail by getting trapped into an infinite loop)
    #  - 2: Munkres-Kuhn (Hungarian) algorithm 
    #  - 21: ChatGPT-generated Munkres-Kuhn (Hungarian) algorithm
    #  - 3: experimental stochastic algorithm, the original version with elimination (known problems)
    #  - 32: experimental stochastic algorithms with all permutations (too expensive)
    #  - 33: the improved stochastic algorithm with good scaling and performance, on par with the mincost
    #  - 4: new, experimental force-based tracking
    dyn_general.update({"state_tracking_algo":-1 })

    #The algorithm to correct phases on adiabatic states
    #
    #Options: 
    #  - 0: no phase correction
    #  - 1: according to our phase correction algorithm [ default ]
    # phase correction doesn't matter, if we use the LD integrator
    dyn_general.update({"do_phase_correction":0 })

    #New phase correction, directly applied to NACs. Intended to be used mostly with state_tracking_algo == 4,
    #although can be useful with other state treacking algorithms. Should not be used together with 
    #`do_phase_correction`
    #Options:
    #  - 0: no correction [ default ]
    #  - 1: do this correction
    dyn_general.update({"do_nac_phase_correction":0 })


    #If set to True (1), we will force the reprojection matrix T_new to be the identity matrix. This effectively
    #removes basis-reprojection (local diabatization) approach and turns on the "naive" approach where
    #no trivial crossings exist.
    #Options:
    #  - 0: No - we do want to use the LD approaches by default. [ default]
    #  - 1: Yes - one may need to turn on additional state tracking and phase correction methods
    dyn_general.update({"assume_always_consistent":0 })

