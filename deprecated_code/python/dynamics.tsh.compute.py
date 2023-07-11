
def init_nuclear_dyn_var(Q, P, M, params, rnd):
    """
    Args:
        Q ( list of doubles ): the mean values of coordinates for all DOFs [ units: a.u.]
        P ( list of doubles ): the mean values of momenta for all DOFs [ units: a.u. ]
        M ( list of doubles ): masses of all nuclear DOFs [ units: a.u. ]

        params ( dictionary ): control parameters
 
            * **params["init_type"]** ( int ): the type of sampling of nuclear DOFs
     
                - 0 : initialize ```ntraj``` identical copies of coordinates and momenta
 
                - 1 : keep all coordinates identical, but sample momenta from the normal 
                    distribution with a given width in each dimension

                - 2 : keep all momenta identical, but sample coordinates from the normal 
                    distribution with a given width in each dimension

                - 3 : sample both coordinates and momenta from the normal 
                    distributions with given widths in each dimension

            * **params["force_constant"]** ( list of double ): force constants involved in the Harmonic
                oscillator model: U = (1/2) * k * x^2, and omega = sqrt( k / m )
                These parameters define the Harmonic oscillator ground state wavefunctions from which
                the coordinates and momenta are sampled. [ units: a.u. = Ha/Bohr^2, default: [0.001] ]
                The length should be consistent with the length of Q, P and M arguments
                The ground state of the corresponding HO is given by:
                 psi_0 = (alp/pi)^(1/4)  exp(-1/2 * alp * x^2), with alp = m * omega / hbar
                The corresponding probability density is distributed as:
                Tully uses psi(x) ~ exp(-(x/sigma)^2) so 1/sigma^2 = alp/2 => alp = 2/sigma^2 = m sqrt (k/m) = sqrt( k * m )

                To make sigma = 20 / p0 => k * m = 4 /sigma^4 = 4/ (20/p0)^4 = 4* p^4 / 20^4 =  0.000025 * p^4 
                So: k =  0.000025 * p^4  / m
                

            * **params["ntraj"]** ( int ): the number of trajectories - the parameter defines the
                number of columns the output matrices will have [ default: 1 ]



        rnd ( Random ): random numbers generator object


    Returns:
        q, p, iM:  where:

            * q ( MATRIX(ndof, ntraj) ) : coordinates for all trajectories
            * p ( MATRIX(ndof, ntraj) ) : momenta for all trajectories
            * iM ( MATRIX(ndof, 1) ) : inverse masses of all DOFs (same across the trajectories)

    """


    critical_params = [ ]
    default_params = { "init_type":0, "force_constant":[0.001], "ntraj":1  }
    comn.check_input(params, default_params, critical_params)
        
    init_type = params["init_type"]  
    force_constant = params["force_constant"]
    ntraj = params["ntraj"]    


    if init_type not in [0, 1, 2, 3]:
        print(F"WARNINIG in init_nuclear_dyn_var: \
              the init_type = {init_type} is not known\
              Allowed values are: [0, 1, 2, 3]" )

    if len(Q)!=len(P):
        print(F"ERROR in init_nuclear_dyn_var: \
              the length of input Q is = {len(Q)}, \
              the length of input P is = {len(P)}, but they should be equal to each other" )
        sys.exit(0)

    if len(Q)!=len(M):
        print(F"ERROR in init_nuclear_dyn_var: \
              the length of input Q is = {len(Q)}, \
              the length of input M is = {len(M)}, but they should be equal to each other" )
        sys.exit(0)

    if len(Q)!=len(force_constant):
        print(F"ERROR in init_nuclear_dyn_var: \
              the length of input Q is = {len(Q)}, \
              the length of input P is = {len(force_constant)}, but they should be equal to each other" )
        sys.exit(0)

    if ntraj < 1:
        print(F"ERROR in init_nuclear_dyn_var: \
              the ntraj is= {ntraj}, should be at least 1" )
        sys.exit(0)



    # At this point, it is safe to define ndof:
    ndof = len(Q)
    q, p, iM = MATRIX(ndof,ntraj), MATRIX(ndof,ntraj), MATRIX(ndof,1)
    for dof in range(ndof):
        iM.set(dof, 0, 1.0/M[dof])


    # Mean values
    mean_q, mean_p = MATRIX(ndof,1), MATRIX(ndof, 1)
    for dof in range(ndof):
        mean_q.set(dof, 0, Q[dof])
        mean_p.set(dof, 0, P[dof])

    # Deviations
    sigma_q, sigma_p = MATRIX(ndof,1), MATRIX(ndof, 1)
    if init_type == 0:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, 0.0)
            sigma_p.set(dof, 0, 0.0)

    elif init_type == 1:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, 0.0)
            sigma_p.set(dof, 0, math.sqrt( 0.5*math.sqrt((force_constant[dof]*M[dof])) ))

    elif init_type == 2:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, math.sqrt( 0.5*math.sqrt(1.0/(force_constant[dof]*M[dof])) ))
            sigma_p.set(dof, 0, 0.0 )

    elif init_type == 3:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, math.sqrt( 0.5*math.sqrt(1.0/(force_constant[dof]*M[dof])) ))
            sigma_p.set(dof, 0, math.sqrt( 0.5*math.sqrt((force_constant[dof]*M[dof])) ))

    # Now sample the values
    tsh.sample(q, mean_q, sigma_q, rnd)
    tsh.sample(p, mean_p, sigma_p, rnd)

    return q, p, iM



def init_electronic_dyn_var(params, rnd):
    """
    Args:

        params ( dictionary ): control parameters
 
            * **params["init_type"]** ( int ): the type of sampling of electronic DOFs
     
                - 0 : initialize all states according to "istate" and sets all
                    amplitudes to 1.0  [ default ]
 
                - 1 : initialize all states according to "istate" but sets 
                    amplitudes to exp(-2*pi*i*rnd), where rnd is a random 
                    number uniformly distributed on the [0, 1] interval

                - 2 : initialize all states according to "istates" - the 
                    integer indices are selected randomly according to populations provided in
                    variable "istates", the amplitudes are set to be sqrt(istates[i]), but
                    their phases are identical (like in the option 0)

                - 3 : initialize all states according to "istates" - the 
                    integer indices are selected randomly according to populations provided in
                    variable "istates", the amplitudes are set to be sqrt(istates[i]), but
                    their phases are set to exp(-2*pi*i*rnd), where rnd is a random 
                    number uniformly distributed on the [0, 1] interval (line in option 1)

            * **params["nstates"]** ( int ): the number of electronic states in the basis
                [ default: 1 ]

            * **params["istate"]** ( int ): the index of the initial electronic state, used 
                only when **params["init_type"]** is 0 or 1, in which case it defines on 
                which state the amplitudes will be initialized [ default: 0]

            * **params["istates"]** ( list of ints ): the list of the populations on all electronic
                states included in the basis, used only when **params["init_type"]** is 2 or 3, 
                in which case it defines on which states the amplitudes will be initialized.
                The length of this list should be consistent with ```nstates``` variable. And the sum
                of all entries should be 1.0 [ default: [1.0], meaning that only the lowest state
                is occupied ]

            * **params["rep"]** ( int ): defines for which repersentation we generate the amplitudes
 
                - 0 : diabatic, the corresponding matrix is initialized, the adiabatic amplitudes are zeroes
                - 1 : adiabatic, the corresponding matrix is initialized, the diabatic amplitudes are zeroes [ default ]

            * **params["ntraj"]** ( int ): the number of trajectories - the parameter defines the
                number of columns the output matrices will have [ default: 1]

            * **params["is_nbra"]** (int): A flag for NBRA type calculations. If it is set to 1 then
                                          the Hamiltonian related properties are only computed for one trajectory [ default : 0]

        rnd ( Random ): random numbers generator object


    Returns:
        Cdia, Cadi, states:  where:

            * Cdia ( CMATRIX(nstates, ntraj) ) : amplitudes on all diabatic states for all trajectories
            * Cadi ( CMATRIX(nstates, ntraj) ) : amplitudes on all adiabatic states for all trajectories
            * states ( list of ints ) : state indices for each trajectory

    """

    # Read the parameters
    critical_params = [  ]
    default_params = { "init_type":0, "nstates":1, "istate":0, "istates":[1.0], "rep":1,  "ntraj":1, "is_nbra":0, "verbosity":0  }
    comn.check_input(params, default_params, critical_params)

    init_type = params["init_type"]
    nstates = params["nstates"]  
    istate = params["istate"]          
    istates = params["istates"]
    rep = params["rep"]  
    ntraj = params["ntraj"]
    is_nbra = params["is_nbra"]
    verbosity = params["verbosity"]

    # Sanity check
    if rep not in [0, 1]:
        print(F"WARNINIG in init_electronic_dyn_var: \
              the rep = {rep} is not known\
              Allowed values are: [0, 1]" )

    if init_type not in [0, 1, 2, 3]:
        print(F"WARNINIG in init_electronic_dyn_var: \
              the init_type = {init_type} is not known\
              Allowed values are: [0, 1, 2, 3]" )

    if ntraj < 1:
        print(F"ERROR in init_electronic_dyn_var: \
              the ntraj is= {ntraj}, should be at least 1" )
        sys.exit(0)
    
    if init_type in [0, 1]:
        if istate >= nstates:
            print(F"ERROR in init_electronic_dyn_var: \
                  the istate is= {istate}, but should be less than {nstates}" )
            sys.exit(0)
             
    if init_type in [2, 3]:
        if len(istates)!= nstates:
            print(F"ERROR in init_electronic_dyn_var: \
                  the istates array is of length {len(istates)}, but should \
                  be of length {nstates}")
            sys.exit(0)
        
        summ = sum(istates)
        if math.fabs( summ - 1.0 ) > 1e-5:
            print(F"ERROR in init_electronic_dyn_var: \
                  the sum of the entries in the istates array is {summ}, but should be 1.0")
            sys.exit(0)

               
                
    # Dynamical variables
    Cdia, Cadi = CMATRIX(nstates, ntraj), CMATRIX(nstates, ntraj)
    states = intList() 


    for traj in range(ntraj):

        if init_type==0:
            if verbosity > 0:
                print(F"======= Initialization type is {init_type} ========\n")
                print(F"setting representation {rep} coefficient C_{istate} to 1.0")
 
            if rep==0:
                Cdia.set(istate, traj, 1.0+0.0j);  
            elif rep==1:
                Cadi.set(istate, traj, 1.0+0.0j);  
            states.append(istate) 

        elif init_type==1:
            if verbosity > 0:
                print(F"======= Initialization type is {init_type} ========\n")
                print(F"setting representation {rep} coefficient to complex C_{istate} such that |C_{istate}|^2 = 1.0")

            ksi = rnd.uniform(0.0, 1.0)
            ampl = math.cos(2*math.pi*ksi) + 1.0j*math.sin(2.0*math.pi*ksi)

            if rep==0:
                Cdia.set(istate, traj, ampl);  
            elif rep==1:
                Cadi.set(istate, traj, ampl);  
            states.append(istate) 


        elif init_type==2:
            if verbosity > 0:
                print(F"======= Initialization type is {init_type} ========\n")
                print(F"setting representation {rep} coefficients C_i for all i to sqrt( target populations) ")

            for state, pop in enumerate(istates):
                ampl = math.sqrt( pop )

                if rep==0:
                    Cdia.set(state, traj, ampl);  
                elif rep==1:
                    Cadi.set(state, traj, ampl);  

            ksi = rnd.uniform(0.0, 1.0)
            states.append(tsh.set_random_state(istates, ksi)) 


        elif init_type==3:

            if verbosity > 0:
                print(F"======= Initialization type is {init_type} ========\n")
                print(F"setting representation {rep} coefficients C_i for all i to complex numbers such that |C_i|^2  = target populations ")

            for state, pop in enumerate(istates):

                ksi = rnd.uniform(0.0, 1.0)
                ampl = math.cos(2*math.pi*ksi) + 1.0j*math.sin(2.0*math.pi*ksi)
                ampl = ampl * math.sqrt( pop )

                if rep==0:
                    Cdia.set(state, traj, ampl);  
                elif rep==1:
                    Cadi.set(state, traj, ampl);  

            ksi = rnd.uniform(0.0, 1.0)
            states.append(tsh.set_random_state(istates, ksi)) 

    projections = CMATRIXList()
    if is_nbra==1:
        projections.append( CMATRIX(nstates, nstates) )
        projections[0].identity()
    else:
        for traj in range(ntraj):
            projections.append( CMATRIX(nstates, nstates) )
            projections[traj].identity()

    if verbosity > 0:
        print("========== Cdia ===============\n")
        Cdia.show_matrix()
        print("========== Cadi ===============\n")
        Cadi.show_matrix()
        print("========== states ===============\n")
        print( Cpp2Py(states) )

        if verbosity > 1:
            print("========== projectors ===============\n")
            for traj in range(ntraj):
                print(F"========== projector for trajectory {traj} =============\n")
                projections[traj].show_matrix()


    return Cdia, Cadi, projections, states



def init_amplitudes(q, Cdia, Cadi, dyn_params, compute_model, model_params, transform_direction=0):
    """
    Args:
        q ( MATRIX(ndof, ntraj) ): coordinates 
        Cdia ( MATRIX(ndia, ntraj) ): amplitudes of diabatic states
        Cadi ( MATRIX(nadi, ntraj) ): amplitudes of adiabatic states

        dyn_params ( Python dictionary ): control of the dynamics
        compute_model ( Python function ): the function that does the calculations
        model_params ( Python dictionary ): parameters of the computational model

        transform_direction ( int ): type of transformation

            - 0: diabatic to adiabatic
            - 1: adiabatic to diabatic


    Returns:
        None: but changes Cadi or Cdia, according to the transform_direction

    """

    # Dimensions
    ndia = Cdia.num_of_rows
    nadi = Cadi.num_of_rows
    nnucl= q.num_of_rows
    ntraj= q.num_of_cols


    # Prepare the Hamiltonian's hierarchy
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.add_new_children(ndia, nadi, nnucl, ntraj)
    ham.init_all(2,1)

    # Compute the Hamiltonian transformation
    update_Hamiltonian_q(dyn_params, q, ham, compute_model, model_params)

    # Do the transformations
    if transform_direction==0:  # dia -> adi
        Cadi = transform_amplitudes(dyn_params, Cdia, ham)
    elif rep_tdse==1:  # adi -> dia
        Cdia = transform_amplitudes(dyn_params, Cadi, ham)

