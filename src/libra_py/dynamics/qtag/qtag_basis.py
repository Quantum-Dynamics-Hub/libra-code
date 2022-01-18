"""
..module:: qtag_basis
  :platform: Unix, Windows
  :synopsis: This module contains update functions for the {q,p,a,s} basis components.

..moduleauthors :: Matthew Dutra
"""

"""

TO DO: Modify the behavior of the propagate functions depending on the type of propagation
instead of choosing different functions

"""

from liblibra_core import *

def frozen(univ, i, param_in, *args, **kwargs):
    """The non-update function.
    Args:
        univ (dictionary): Dictionary containing various system parameters.

        i (integer) : The DoF being considered

        param_in (MATRIX): 1-by-ntraj basis parameter matrix for q, p, a, or s.
    Returns:
        param_in (MATRIX): The same input parameter matrix q, p, a, or s.

    """

    return(param_in)


def q_update(univ, ndim, q_in, mom_in, mss_prop):
    """Returns the updated position *q_out* from an input position *q_in* and the momentum *mom*, 
       calculated via q_out = q_in+mom*dt/m (i.e. symplectic).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        ndim (integer): The dimension along which q is being updated.

        q_in (MATRIX): The 1-by-ntraj basis position MATRIX q to be updated.

        mom_in (MATRIX): The 1-by-ntraj momentum MATRIX p corresponding to q.

        mss_prop (function object): Currently unused.

    Returns:
        q_out (MATRIX): Updated 1-by-ntraj basis position MATRIX q.

    """

#	if mss_prop=='two_surf':
#		ntraj=2*univ['ntraj']
#	else:

    ntraj=univ['ntraj']
    dt = univ['dt']
    mass = univ['mass']

    return (q_in+mom_in*dt/mass[ndim])


def p_update(mom):
    """Returns the updated basis parameter *p_out* from the momentum *mom*. Currently these two things are equal, but they don't necessarily need to be.

    Args:
        mom (MATRIX): The 1-by-ntraj basis momentum MATRIX.

    Returns:
        p_out (MATRIX): Updated 1-by-ntraj basis parameter MATRIX p.

    """

    return 1.0*mom


def a_update(univ, ndim, a_in, gmom, mss_prop):
    """Returns the updated basis width *a_out*, calculated from the old basis width *a_in* and the momentum gradient *gmom*.

    Args:

        univ (dictionary): Dictionary containing various system parameters.

        ndim (integer): The dimension along which a is being updated.
 
        a_in (MATRIX): The 1-by-ntraj basis width MATRIX a to be updated.

        gmom (MATRIX): The 1-by-ntraj momentum gradient MATRIX storing the values 
            of (dp/dx) at the basis positions corresponding to MATRIX a.

        mss_prop (function object): Currently unused.

    Returns:
        a_out (MATRIX): Updated 1-by-ntraj basis width MATRIX a.

    """

#	if mss_prop=='two_surf':
#		ntraj=2*univ['ntraj']
#	else:

    ntraj=univ['ntraj']
    dt=univ['dt']
    mass=univ['mass']

    a_out, an1 = MATRIX(1,ntraj), MATRIX(1,ntraj)
    an1.dot_product(a_in, gmom)

    return a_in - 2.0*an1*dt/mass[ndim]


def cls_force_q(univ,q_in,p_in, model, model_params):
    """Returns the updated position *q_out* from an input position *q_in* and the momentum *mom*, 
       calculated via q_out = q_in+mom*dt/m (i.e. symplectic).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        q_in (MATRIX): The 1-by-ntraj basis position MATRIX q to be updated.

        p_in (MATRIX): The 1-by-ntraj momentum MATRIX p corresponding to q.

        model_params (dictionary): Dictionary containing potential parameters.

    Returns:
        q_out (MATRIX): Updated 1-by-ntraj basis position MATRIX q.

		p_out (MATRIX): Updated 1-by-ntraj basis momentum MATRIX p.
	"""
    ntraj=univ['ntraj']
    ndof=univ['ndof']
    dt=univ['dt']
    mass=univ['mass']
    pot=univ['pot_fxn']

    q_out=MATRIX(ndof,ntraj)
    p_out=MATRIX(ndof,ntraj)

#    q_out = q_in + iM * p_in * dt
    for i in range(ntraj):
        x, dvx, d2vx = pot(ndof,q_in.col(i),2,model_params)

        """

        d2vx is not used!  

        Make the dvx and d2vx MATRIX(ndof, ntraj) objects - then simplify
        the double loop
        """
   
        for j in range(ndof):
            p_out.set(j,i, p_in.get(j,i)-dvx[j]*dt )
            q_out.set(j,i, q_in.get(j,i)+(p_in.get(j,i)*dt-dvx[j]*dt**2)/mass[j])

    return (q_out,p_out)


def param_check(basis, cls_chk):
     """Checks the type of propagation requested for each basis parameter {q,p,a,s} from the basis dictionary.

     Args:
         basis (dictionary): The input dictionary containing keywords for the basis parameter propagation.

         cls_chk (character): Keyword determining if classical propagation is used in the mss.

     Returns:
         props (list of function objects): List containing the functions for {q,p,a,s} parameter updates.

    """

    if basis['qtype'] == 'adpt':
        qprop=q_update    
    elif basis['qtype'] == 'frzn':
        qprop=frozen
    else:
        sys.exit("Unrecognized keyword in basis qtype!")


    if basis['ptype'] == 'adpt':
        pprop=p_update
    elif basis['ptype'] == 'frzn':
        pprop=frozen
    else:
        sys.exit("Unrecognized keyword in basis ptype!")

    if basis['atype'] == 'adpt':
        aprop=a_update
    elif basis['atype'] == 'frzn':
        aprop=frozen
    else:
        sys.exit("Unrecognized keyword in basis atype!")

    if basis['stype'] == 'adpt':
        print("Adaptable s not implemented yet! Reverting to frozen...")
        sprop=frozen
    elif basis['stype'] == 'frzn':
        sprop=frozen
    else:
        sys.exit("Unrecognized keyword in basis stype!")

    if cls_chk == 'cls_force':
        cls_prop=cls_force_q
    else:
        cls_prop=frozen

    props = [qprop,pprop,aprop,sprop,cls_prop]

    return (props)
