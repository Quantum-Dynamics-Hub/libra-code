import sys

from . import qtag_init
from . import qtag_basis
from . import qtag_ham
from . import qtag_mom
from . import qtag_prop

def user_input(dyn_params,qtag_params,model_params):
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

    ndof=dyn_params['ndof']

    if len(dyn_params['mass']) != ndof:
        sys.exit("List of masses does not match given number of DoFs!")

    if len(qtag_params['wfq0']) != ndof:
        sys.exit("Initial WF q's do not match given number of DoFs!")

    if len(qtag_params['wfp0']) != ndof:
        sys.exit("Initial WF p's do not match given number of DoFs!")

    if len(qtag_params['wfa0']) != ndof:
        sys.exit("Initial WF a's do not match given number of DoFs!")

    if len(qtag_params['wfs0']) != ndof:
        sys.exit("Initial WF s's do not match given number of DoFs!")

    if len(qtag_params['grid_dims']) != ndof and qtag_params['init_placement']=='grid':
        sys.exit("List of grid_dims does not match given number of DoFs!")

    if len(qtag_params['a0']) != ndof:
        sys.exit("List of a0 values does not match given number of DoFs!")

#    if model['coupling'] == 'exact':
#        model['coupling'] = 'exact_gauss_cpl'

#    if model['coupling'] == 'none':
#        model['coupling'] = 'no_coupling'

#    if model['pot_type'] == 'HO':
#        for i in range(ndof):
#            model_params['d2'][i]=model_params['d2'][i]*np.sqrt(model_params['k1'][i])

    return()


#def assign_fobj(qtag_params):

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

#    return(initialize,props,vapprox,mom_calc,propagate)



def set_basis_updates(qtag_params):
    """Checks the type of propagation requested for each basis parameter {q,p,a,s} from the basis dictionary.

    Args:
        basis (dictionary): The input dictionary containing keywords for the basis parameter propagation.

        cls_chk (character): Keyword determining if classical propagation is used in the mss.

    Returns:
        props (list of function objects): List containing the functions for {q,p,a,s} parameter updates.
    """

    
    if qtag_params['basis_qtype'] == 1:
        qprop=qtag_basis.q_update
    elif qtag_params['basis_qtype'] == 0:
        qprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis qtype!")


    if qtag_params['basis_ptype'] == 1:
        pprop=qtag_basis.p_update
    elif qtag_params['basis_ptype'] == 0:
        pprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis ptype!")

    if qtag_params['basis_atype'] == 1:
        aprop=qtag_basis.a_update
    elif qtag_params['basis_atype'] == 0:
        aprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis atype!")

    if qtag_params['basis_stype'] == 1:
        print("Adaptable s not implemented yet! Reverting to frozen...")
        sprop=qtag_basis.frozen
    elif qtag_params['basis_stype'] == 0:
        sprop=qtag_basis.frozen
    else:
        sys.exit("Unrecognized keyword in basis stype!")

    if qtag_params['prop_method'] == 'cls_force':
        cls_prop=qtag_basis.cls_force_q
    else:
        cls_prop=qtag_basis.frozen

    basis_props = [qprop,pprop,aprop,sprop,cls_prop]

    return (basis_props)
