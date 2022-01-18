def BAT(univ,nstate,qpasi,qpasj,params,libra_model):
	"""Returns the (complex) value for the potential *v* on an energetic surface specified by *nsurf* from two basis functions defined by their parameters *qpasi* and *qpasj*, respectively. The computation employs the Bra-ket Averaged Taylor expansion (BAT) about each basis center, which requires a potential function *pot* as well as its first derivative. The potential parameters are stored in the dict 'model_params' as defined in qtag_config.

        Args:
		univ (dictionary): Dictionary containing various system parameters.

		nstate (integer): Integer specifying the potential surface to be calculated (0 = ground, 1 = first excited, ...)

                qpasi (list): The ndof-by-4 parameter list of the i-th basis function. Each entry is an ndof-by-1 column MATRIX.

                qpasj (list): The ndof-by-4 parameter list of the j-th basis function. Each entry is an ndof-by-1 column MATRIX.

		params (dictionary): Dictionary containing the potential parameters.

                libra_model (function object): Function object containing the Libra potential model for the ground and excited states.

        Returns:
                v (complex): The complex potential v computed via the bra-ket averaged Taylor expansion between basis functions i and j.
        """

    ndof=univ['ndof']
    nstates=univ['nstates']

    q1,p1,a1,s1=qpasi[0], qpasi[1], qpasi[2], qpasi[3]
    q2,p2,a2,s2=qpasj[0], qpasj[1], qpasj[2], qpasj[3]
    
    obj1 = libra_model(q1, params, full_id)
    obj2 = libra_model(q2, params, full_id)

    dvx1_sum=0.0;dvx2_sum=0.0
    vx1=obj1.ham_dia.get(nstate,nstate)
    vx2=obj2.ham_dia.get(nstate,nstate)

    for i in range(ndof):
    	dvx1 = obj1.d1ham_dia[i].get(nstate,nstate)
        dvx2 = obj2.d1ham_dia[i].get(nstate,nstate)

    	q1_rr1_q2=complex(a2.get(i)*(q2.get(i)-q1.get(i)),p2.get(i)-p1.get(i))/(a1.get(i)+a2.get(i))
    	q1_rr2_q2=complex(a1.get(i)*(q1.get(i)-q2.get(i)),p2.get(i)-p1.get(i))/(a1.get(i)+a2.get(i))

    	dvx1_sum+=dvx1*q1_rr1_q2
    	dvx2_sum+=dvx2*q1_rr2_q2

    v=0.5*(vx1+vx2+dvx1_sum+dvx2_sum)

    return(v)

def LHA(univ,nstate,qpasi,qpasj,params,libra_model):
	"""Returns the (complex) value for the potential *v* on an energetic surface specified by *nsurf* from two basis functions defined by their parameters *qpasi* and *qpasj*, respectively. The computation employs the Local Harmonic Approximation (LHA), which requires a potential function *pot* as well as its first and second derivatives. The potential parameters are stored in the dict 'model_params' as defined in qtag_config.

        Args:
		univ (dictionary): Dictionary containing various system parameters.

		nstate (integer): Integer specifying the potential surface to be calculated (0 = ground, 1 = first excited, ...)

                qpasi (list): The ndof-by-4 parameter list of the i-th basis function. Each entry is an ndof-by-1 column MATRIX.

                qpasj (list): The ndof-by-4 parameter list of the j-th basis function. Each entry is an ndof-by-1 column MATRIX.

		params (dictionary): Dictionary containing the potential parameters.

                libra_model (function object): Function object containing the Libra potential model for the ground and excited states.

        Returns:
                v (complex): The complex potential v computed via the local harmonic approximation between basis functions i and j.
        """

    ndof=univ['ndof']
    nstates=univ['nstates']

    obj1, obj2 = CMATRIX(nstates,nstates), CMATRIX(nstates,nstates)

    q1,p1,a1,s1=qpasi[0], qpasi[1], qpasi[2], qpasi[3]
    q2,p2,a2,s2=qpasj[0], qpasj[1], qpasj[2], qpasj[3]
    
    obj1=libra_model(q1, params, full_id)
    obj2=libra_model(q2, params, full_id)

    v=complex(0.0,0.0)
    vx1=obj1.ham_dia.get(nstate,nstate)
    vx2=obj2.ham_dia.get(nstate,nstate)
    for i in range(ndof):
        dvx1,d2vx1=obj1.d1ham_dia[i].get(nstate,nstate),obj1.d2ham_dia[i].get(nstate,nstate)
    	        
        z=complex(a1.get(i)*q1.get(i)+a2.get(i)*q2.get(i),p2.get(i)-p1.get(i))/(a1.get(i)+a2.get(i))
        vv0=-dvx1[i]*q1.get(i)+d2vx1[i]/2.0*q1.get(i)**2
        vv1=-d2vx1[i]*q1.get(i)+dvx1[i]
        vv2=d2vx1[i]/2.0
        v+=0.5*(vv0+vv1*z+vv2*(z**2+1.0/(a1.get(i)+a2.get(i))))

    for i in range(ndof):
        dvx2,d2vx2=obj2.d1ham_dia[0].get(nstate,nstate),obj2.d2ham_dia[0].get(nstate,nstate)

        z=complex(a1.get(i)*q1.get(i)+a2.get(i)*q2.get(i),p2.get(i)-p1.get(i))/(a1.get(i)+a2.get(i))
        vv0=-dvx2[i]*q2.get(i)+d2vx2[i]/2.0*q2.get(i)**2
        vv1=-d2vx2[i]*q2.get(i)+dvx2[i]
        vv2=d2vx2[i]/2.0
        v+=0.5*(vv0+vv1*z+vv2*(z**2+1.0/(a1.get(i)+a2.get(i))))

    v+=0.5*(vx1+vx2)
    	    
    return(v)

def super_hamiltonian(univ,qpas,ov,nstates,vcalc,libra_model,model_params):
    """Calculates the single-surface Hamiltonian matrix elements H_ij=<gi|KE+V|gj>, computed using the basis parameters stored in *qpas*. This requires the single-surface overlap matrix *ov* as well as the potential function *pot*, which is specified in qtag_config. The surface is designated by *nsurf*. Returns the single-surface Hamiltonian *H*.

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        qpas (list): List of {q,p,a,s} MATRIX objects and a list of integers for their assigned states. (0=ground, 1=first excited...)

        ov (CMATRIX): The full ntraj-by-ntraj super-overlap matrix.

        nstates (integer): Integer specifying the number of states (surfaces) in the system.

	vcalc (function object): Function object containing the method for approximating potential matrix elements (LHA or BAT).

        libra_model (function object): Function object containing the Libra potential model for various surfaces.

	model_params (dictionary): Dictionary containing the potential parameters to be implemented in libra_model.

    Returns:
        H (CMATRIX): The full ntraj-by-ntraj super-Hamiltonian.

    """

    ndof,ntraj=univ['ndof'],univ['ntraj']
    mass=univ['mass']
    H = CMATRIX(ntraj,ntraj)
    vtot = CMATRIX(nstates,nstates)
    iMasses=MATRIX(ndof,1)

    for i in range(ndof):
    	iMasses.set(i,0,mass[i])

#Extract the components of the qpas object into their constituent parts: q, p, a, s, surface IDs.
    qvals,pvals,avals,svals=MATRIX(qpas[0]),MATRIX(qpas[1]),MATRIX(qpas[2]),MATRIX(qpas[3])
    surf_ids=qpas[4]

#Loop over the number of states, initializing the ntraj_on_surf list so that a lower bound of 0 for the first state can be used...
    ntraj_on_surf=[0]
    for n in range(nstates):
#Extract the relevent trajectories based on the state (n) under consideration...
    	ntraj_on_surf.append(surf_ids.count(nstates))
    	pop_submatrix(qvals,qvals_surf,[dof for dof in range(ndof)],[ntraj_on_surf[nstates]+traj for traj in range(ntraj_on_surf[nstates+1])])
        pop_submatrix(pvals,pvals_surf,[dof for dof in range(ndof)],[ntraj_on_surf[nstates]+traj for traj in range(ntraj_on_surf[nstates+1])])
        pop_submatrix(avals,avals_surf,[dof for dof in range(ndof)],[ntraj_on_surf[nstates]+traj for traj in range(ntraj_on_surf[nstates+1])])
        pop_submatrix(svals,svals_surf,[dof for dof in range(ndof)],[ntraj_on_surf[nstates]+traj for traj in range(ntraj_on_surf[nstates+1])])

#Loop over the trajectories on the specific state, computing the Hamiltonian elements along the (block) diagonal...
        for i in range(ntraj_on_surf[n]):
            q1,p1,a1,s1=qvals_surf.col(i),pvals_surf.col(i),avals_surf.col(i),svals_surf.col(i)
            for j in range(ntraj_on_surf[n]):
                q2,p2,a2,s2=qvals_surf.col(j),pvals_surf.col(j),avals_surf.col(j),svals_surf.col(j)

                ke=gwp_kinetic(q1,p1,s1,a1/2,q2,p2,s2,a2/2,iMasses)
                qpasi=[q1,p1,a1,s1];qpasj=[q2,p2,a2,s2]
                v=vcalc(univ,n,qpasi,qpasj,model_params,libra_model)
                H.set(i+ntraj_on_surf[n],j+ntraj_on_surf[n],ke+v*ov.get(i+ntraj_on_surf[n],j+ntraj_on_surf[n]))

    return(H)

def basis_diag(m,dt,H,ov,b):
    """Returns the updated basis coefficients for both surfaces, stored in a single vector *c_new*, computed as c_new=Z*exp(-i*dt*eps)*Z_dagger*b. The variables eps and Z are the eigenvalues and eigenvectors obtained from diagonalizing the full *m*-by-*m* Hamiltonian matrix *H* using the solve_eigen internal function. Note that the projection vector *b* and the full overlap matrix *ov* are required.

    Args:
        m (integer): Dimension of the eigenvalue and eigenvector matrices. Usually equal to 2*ntraj.

        dt (real): The timestep parameter.

        H (CMATRIX): Total system Hamiltonian, dimension 2*ntraj-by-2*ntraj.

        ov (CMATRIX): Total overlap matrix, dimension 2*ntraj-by-2*ntraj.

        b (CMATRIX): Projection vector, dimensioned 2*ntraj-by-1.

    Returns:

        c_new (CMATRIX): Updated basis coefficient vector for both surfaces. Computed as c_new=Z*exp(-i*dt*eps)*Z_dagger*b.

    """

    evals=CMATRIX(m,m)
    evecs=CMATRIX(m,m)
    solve_eigen(H,ov,evals,evecs,0)
    ct=evecs.T().conj()
    c_new=evecs*(exp_(evals,-dt*1.0j))*ct*b

    return(c_new)