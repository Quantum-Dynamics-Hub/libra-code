from liblibra_core import *

def vapprox(ndof,nstates,n1,n2,qi,pi,ai,si,qj,pj,aj,sj,approx_method,params,compute_model):

    if approx_method == 0:
        v = BAT(ndof,nstates,n1,n2,qi,pi,ai,si,qj,pj,aj,sj,params,compute_model)

    elif approx_method == 1:
        v = LHA(ndof,nstates,n1,n2,qi,pi,ai,si,qj,pj,aj,sj,params,compute_model)

    return(v)

def BAT(ndof,nstates,n1,n2,qi,pi,ai,si,qj,pj,aj,sj,params,compute_model):
    """Returns the (complex) value for the potential *v* on an energetic surface specified by *nsurf* from two basis 
       functions defined by their parameters *qpasi* and *qpasj*, respectively. The computation employs the Bra-ket 
       Averaged Taylor expansion (BAT) about each basis center, which requires a potential function *pot* as well as 
       its first derivative. The potential parameters are stored in the dict 'model_params' as defined in qtag_config.

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        nstate (integer): Integer specifying the potential surface to be calculated (0 = ground, 1 = first excited, ...)

        qpasi (list): The ndof-by-4 parameter list of the i-th basis function. Each entry is an ndof-by-1 column MATRIX.

        qpasj (list): The ndof-by-4 parameter list of the j-th basis function. Each entry is an ndof-by-1 column MATRIX.

        params (dictionary): Dictionary containing the potential parameters.

        libra_model (function object): Function object containing the Libra potential model for the ground and excited states.

    Returns:
        complex : v - The complex potential v computed via the bra-ket averaged Taylor expansion between basis functions i and j.

    """
    obj1 = CMATRIX(nstates,nstates)
    obj2 = CMATRIX(nstates,nstates)

    obj1 = compute_model(qi, params, full_id=0)
    obj2 = compute_model(qj, params, full_id=0)

    vx1=obj1.ham_dia.get(n1,n2)
    vx2=obj2.ham_dia.get(n1,n2)
    v = 0.5 * (vx1 + vx2)    

    dvx1_sum, dvx2_sum = 0.0, 0.0
    for dof in range(ndof):
        dvx1 = obj1.d1ham_dia[dof].get(n1,n2)
        dvx2 = obj2.d1ham_dia[dof].get(n1,n2)

        q1 = qi.get(dof)
        q2 = qj.get(dof)

        p1 = pi.get(dof)
        p2 = pj.get(dof)

        a1 = ai.get(dof)
        a2 = aj.get(dof)

        dq = q2 - q1
        dp = p2 - p1
        denom = a1 + a2 

        q1_rr1_q2 = complex(a2*dq, dp) / denom
        q1_rr2_q2 = complex(-a1*dq, dp) / denom

        dvx1_sum += dvx1*q1_rr1_q2
        dvx2_sum += dvx2*q1_rr2_q2

    v += 0.5*(dvx1_sum + dvx2_sum)

    return(v)


def LHA(ndof,nstates,n1,n2,qi,pi,ai,si,qj,pj,aj,sj,params,compute_model):
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

    obj1 = CMATRIX(nstates,nstates)
    obj2 = CMATRIX(nstates,nstates)

#    full_id = None
    obj1 = compute_model(qi,params,full_id=None)
    obj2 = compute_model(qj,params,full_id=None)

    vxi = obj1.ham_dia.get(n1,n2)
    vxj = obj2.ham_dia.get(n1,n2)

    v = complex(0.0,0.0)
    for dof in range(ndof):
        q1 = qi.get(dof)
        q2 = qj.get(dof)

        p1 = pi.get(dof)
        p2 = pj.get(dof)

        a1 = ai.get(dof)
        a2 = aj.get(dof)

        s1 = si.get(dof)
        s2 = sj.get(dof)

        aqp = a1*q1+a2*q2
        dp = p2-p1
        asum = a1+a2
        z = complex(aqp,dp)/asum

        dvxi = obj1.d1ham_dia[dof].get(n1,n2)
        dvxj = obj2.d1ham_dia[dof].get(n1,n2)

        d2vxi = obj1.d2ham_dia[dof].get(n1,n2)
        d2vxj = obj2.d2ham_dia[dof].get(n1,n2)

        vv0i = -dvxi*q1+d2vxi/2.0*q1**2
        vv0j = -dvxj*q2+d2vxj/2.0*q2**2

        vv1i = -d2vxi*q1+dvxi
        vv1j = -d2vxj*q2+dvxj

        vv2i = d2vxi/2.0
        vv2j = d2vxj/2.0

        v+=0.5*(vv0i+vv0j+(vv1i+vv1j)*z+(vv2i+vv2j)*(z**2+1.0/asum))

    v+=0.5*(vx1+vx2)

    return(v)


def super_hamiltonian(ndof,ntraj,nstates,mass,qpas,approx_method,compute_model,model_params):
    """Calculates the single-surface Hamiltonian matrix elements H_ij=<gi|KE+V|gj>, computed using the basis 
       parameters stored in *qpas*. This requires the single-surface overlap matrix *ov* as well as the potential 
       function *pot*, which is specified in qtag_config. The surface is designated by *nsurf*. 
       Returns the single-surface Hamiltonian *H*.

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        qpas (list): List of {q,p,a,s} MATRIX objects and a list of integers for their assigned states. (0=ground, 1=first excited...)

        ov (CMATRIX): The full ntraj-by-ntraj super-overlap matrix.

        nstates (integer): Integer specifying the number of states (surfaces) in the system.

	vapprox (function object): Function object containing the method for approximating potential matrix elements (LHA or BAT).

        libra_model (function object): Function object containing the Libra potential model for various surfaces.

	model_params (dictionary): Dictionary containing the potential parameters to be implemented in libra_model.

    Returns:
        H (CMATRIX): The full ntraj-by-ntraj super-Hamiltonian.

    """

    super_ham = CMATRIX(ntraj,ntraj)

    #Extract the components of the qpas object into their constituent parts: q, p, a, s, surface IDs.
    qvals,pvals,avals,svals=MATRIX(qpas[0]),MATRIX(qpas[1]),MATRIX(qpas[2]),MATRIX(qpas[3])
    surf_ids=qpas[4]

    dof_dim = list(range(ndof))    

    for n1 in range(nstates):

        #Extract the relevent trajectories for the n1-th state...
        traj_on_surf_n1 = [index for index, traj_id in enumerate(surf_ids) if traj_id == n1]
        ntraj_on_surf_n1 = len(traj_on_surf_n1)

        qvals_surf_n1 = MATRIX(ndof, ntraj_on_surf_n1)
        pvals_surf_n1 = MATRIX(ndof, ntraj_on_surf_n1)
        avals_surf_n1 = MATRIX(ndof, ntraj_on_surf_n1)
        svals_surf_n1 = MATRIX(ndof, ntraj_on_surf_n1)

        pop_submatrix(qvals,qvals_surf_n1, dof_dim, traj_on_surf_n1)
        pop_submatrix(pvals,pvals_surf_n1, dof_dim, traj_on_surf_n1)
        pop_submatrix(avals,avals_surf_n1, dof_dim, traj_on_surf_n1)
        pop_submatrix(svals,svals_surf_n1, dof_dim, traj_on_surf_n1)

        for n2 in range(n1+1):

            #Extract the relevant trajectories for the n2-th state...
            traj_on_surf_n2 = [index for index, traj_id in enumerate(surf_ids) if traj_id == n2]
            ntraj_on_surf_n2 = len(traj_on_surf_n2)

            qvals_surf_n2 = MATRIX(ndof,ntraj_on_surf_n2)
            pvals_surf_n2 = MATRIX(ndof,ntraj_on_surf_n2)
            avals_surf_n2 = MATRIX(ndof,ntraj_on_surf_n2)
            svals_surf_n2 = MATRIX(ndof,ntraj_on_surf_n2)

            pop_submatrix(qvals,qvals_surf_n2, dof_dim, traj_on_surf_n2)
            pop_submatrix(pvals,pvals_surf_n2, dof_dim, traj_on_surf_n2)
            pop_submatrix(avals,avals_surf_n2, dof_dim, traj_on_surf_n2)
            pop_submatrix(svals,svals_surf_n2, dof_dim, traj_on_surf_n2)


            n12_mat = CMATRIX(ntraj_on_surf_n1, ntraj_on_surf_n2)

            n12_mat = state_hamiltonian(qvals_surf_n1, pvals_surf_n1, avals_surf_n1, svals_surf_n1, n1, \
                                        qvals_surf_n2, pvals_surf_n2, avals_surf_n2, svals_surf_n2, n2, \
                                        nstates, mass, approx_method, compute_model, model_params)


            push_submatrix(super_ham, n12_mat, traj_on_surf_n1, traj_on_surf_n2)
            if n1!=n2:
                push_submatrix(super_ham, n12_mat.H(), traj_on_surf_n2, traj_on_surf_n1)


    return super_ham

def state_hamiltonian(qvals_surf_n1, pvals_surf_n1, avals_surf_n1, svals_surf_n1, n1, \
                      qvals_surf_n2, pvals_surf_n2, avals_surf_n2, svals_surf_n2, n2, \
                      nstates, mass, approx_method, compute_model, model_params):
    """
    This function computes a block hamiltonian for the two sets of trajectories - belonging to two surfaces
    (including to the same surface)

    
    """

    ndof = qvals_surf_n1.num_of_rows
    ntraj1 = qvals_surf_n1.num_of_cols  # the number of trajectories on state 1 (itraj)
    ntraj2 = qvals_surf_n2.num_of_cols  # the number of trajectories on state 2 (jtraj)
    state_ham = CMATRIX(ntraj1, ntraj2)

    invM = MATRIX(ndof,1)
    for dof in range(ndof):
        invM.set(dof,0,1.0/mass[dof])


    for i in range(ntraj1):
        qi = qvals_surf_n1.col(i)
        pi = pvals_surf_n1.col(i)
        ai = avals_surf_n1.col(i)
        si = svals_surf_n1.col(i)

        for j in range(ntraj2):
            qj = qvals_surf_n2.col(j)
            pj = pvals_surf_n2.col(j)
            aj = avals_surf_n2.col(j)
            sj = svals_surf_n2.col(j)

            v = vapprox(ndof,nstates,n1,n2,qi,pi,ai,si,qj,pj,aj,sj,approx_method,model_params,compute_model)
            ov = gwp_overlap(qi,pi,si,0.5*ai,qj,pj,sj,0.5*aj)
            state_ham.set(i, j, v*ov)

            # trajectories on the same state - then add kintetic energy term
            if n1 == n2:
                ke = gwp_kinetic(qi,pi,si,0.5*ai,qj,pj,sj,0.5*aj,invM)
                state_ham.add(i, j, ke)

    return state_ham



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

    evals = CMATRIX(m,m)
    evecs = CMATRIX(m,m)
    solve_eigen(H,ov,evals,evecs,0)
    ct = evecs.H()
    c_new = evecs*(exp_(evals,-dt*1.0j))*ct*b

    return(c_new)
