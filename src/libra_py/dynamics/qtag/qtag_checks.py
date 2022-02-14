import sys

from . import qtag_init
from . import qtag_basis
from . import qtag_ham
from . import qtag_mom
from . import qtag_prop

def user_input(univ,wf0,traj0,mss,model,model_params):
    """Runs checks to ensure the validity of the user input.

    Args:
        univ (dictionary): Dictionary containing system parameters.

        wf0 (dictionary): Dictionary containing initial wavepacket conditions.

        traj0 (dictionary): Dictionary containing initialization parameters and keywords.

        mss (dictionary): Dictionary containing multi-surface scheme keywords.

        model (dictionary): Dictionary containing containing keywords for the potential energy model.

        model_params (dictionary): Dictionary containing the potential energy parameters.

    Returns:
        None.
    """

    ntraj,ndof=univ['ntraj'],univ['ndof']

    if len(univ['mass']) != ndof:
        sys.exit("List of masses does not match given number of DoFs!")

    if len(wf0['q']) != ndof:
        sys.exit("Initial WF q's do not match given number of DoFs!")

    if len(wf0['p']) != ndof:
        sys.exit("Initial WF p's do not match given number of DoFs!")

    if len(wf0['a']) != ndof:
        sys.exit("Initial WF a's do not match given number of DoFs!")

    if len(wf0['s']) != ndof:
        sys.exit("Initial WF s's do not match given number of DoFs!")

    if len(traj0['grid_dims']) != ndof and traj0['placement']=='grid':
        sys.exit("List of grid_dims does not match given number of DoFs!")

    if len(traj0['a0']) != ndof:
        sys.exit("List of a0 values does not match given number of DoFs!")

    if model['rep'] == 'diabatic':
        model['rep'] = 'diab'

    if model['rep'] == 'adiabatic':
        model['rep'] = 'adiab'

    if model['coupling'] == 'exact':
        model['coupling'] = 'exact_gauss_cpl'

    if model['coupling'] == 'none':
        model['coupling'] = 'no_coupling'

#    if model['pot_type'] == 'HO':
#        for i in range(ndof):
#            model_params['d2'][i]=model_params['d2'][i]*np.sqrt(model_params['k1'][i])

#    if ntraj%nstates != 0:
#        print("Number of trajectories does not evenly distribute across states!")
#        dyn_params['ntraj']+=(nstates-ntraj%nstates)
#        print("Adjusting by adding ", nstates-ntraj%nstates, " trajectories!")

    return()

def assign_fobj(qtag_params):

    basis = qtag_params['basis']
    prop_method = str(qtag_params['mss']['prop_method'])
    placement = str(qtag_params['traj0']['placement'])

#Assemble propagation types for {q,p,a,s} basis parameters...
    props = basis_fobj(basis,prop_method)

#Obtain the keywords for the calculation type from the various dictionaries...
    try:
        initialize = getattr(qtag_init,placement)
    except AttributeError:
        sys.exit("Error in traj0 dictionary: 'placement' keyword not recognized!")

    try:
        vapprox = getattr(qtag_ham,str(qtag_params['v_approx']))
    except AttributeError:
        sys.exit("Error in qtag dictionary: 'v_approx' keyword not recognized!")

#    try:
#        pot=getattr(qtag_pots,model_name)
#        if coupling_name == 'LHA' or coupling_name == 'BAT':
#            cplg=getattr(qtag_ham,coupling_name)
#        else:
#            cplg=getattr(qtag_pots,coupling_name)
#        univ['pot_fxn']=pot;univ['cplg_fxn']=cplg
#    except AttributeError:
#        print("Potential Model: ",model_name)
#        print("Coupling: ",coupling_name)
#        sys.exit("Error in model dictionary, check that pot_type and rep are compatible!")

    try:
        mom_calc=getattr(qtag_mom,str(qtag_params['mom_params']['adjust']))
    except AttributeError:
        sys.exit("Error in mom_params dictionary: 'adjust' keyword not recognized!")

    try:
        propagate=getattr(qtag_prop,str(qtag_params['mss']['prop_method']))
    except AttributeError:
        sys.exit("Error in mss dictionary: 'prop_method' keyword not recognized!")

    return(initialize,props,vapprox,mom_calc,propagate)

def basis_fobj(basis, cls_chk):
    """Checks the type of propagation requested for each basis parameter {q,p,a,s} from the basis dictionary.

    Args:
        basis (dictionary): The input dictionary containing keywords for the basis parameter propagation.

        cls_chk (character): Keyword determining if classical propagation is used in the mss.

    Returns:
        props (list of function objects): List containing the functions for {q,p,a,s} parameter updates.
    """

    if basis['qtype'] == 'adpt':
        qprop=qtag_basis.q_update
    elif basis['qtype'] == 'frzn':
        qprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis qtype!")


    if basis['ptype'] == 'adpt':
        pprop=qtag_basis.p_update
    elif basis['ptype'] == 'frzn':
        pprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis ptype!")

    if basis['atype'] == 'adpt':
        aprop=qtag_basis.a_update
    elif basis['atype'] == 'frzn':
        aprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis atype!")

    if basis['stype'] == 'adpt':
        print("Adaptable s not implemented yet! Reverting to frozen...")
        sprop=qtag_basis.frozen
    elif basis['stype'] == 'frzn':
        sprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis stype!")

    if cls_chk == 'cls_force':
        cls_prop=qtag_basis.cls_force_q
    else:
        cls_prop=qtag_basis.frozen

    props = [qprop,pprop,aprop,sprop,cls_prop]

    return (props)
