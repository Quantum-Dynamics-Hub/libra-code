"""
..module:: QTAG_assembler
  :platform: Unix, Windows
  :synopsis: This module contains assembly routines used to select appropriate functions based on the input parameters stored in QTAG_config.

..moduleauthors :: Matthew Dutra
"""

"""
#.. py:function:: assemble_init(traj0)
#   Returns a function object to initialize the basis parameters 
#   {q,p,a,s} based on the 'traj0' dict specified in QTAG_config.
"""

def assemble_init(traj0):
	"""Returns a function object to initialize the basis parameters {q,p,a,s} based on the 'traj0' dict specified in QTAG_config.

	Args:
		traj0 (dictionary): Dictionary containing keywords designating how the basis should be initialized.

	Returns:
		init_basis (function object): Function for initializing basis parameters.
	"""

	if traj0['type'] == 'grid':
		from QTAG_init import grid as init_basis
	elif traj0['type'] == 'gaus':
		from QTAG_init import gaus as init_basis
	else:
		sys.exit("Unrecognized keyword in traj0 type!")
	return(init_basis)

"""
.. py:function:: assemble_basis(basis)
   Returns function objects for adaptable and frozen {q,p,a,s} 
   updates based on the 'basis' dict specified in QTAG_config.
"""

def assemble_basis(basis):
	"""Returns function objects for adaptable and frozen {q,p,a,s} updates based on the 'basis' dict specified in QTAG_config.

	Args:
                basis (dictionary): Dictionary containing keywords designating basis parameters as adaptable or frozen.
        Returns:
                props (list): List of function objects for computing basis updates.
        """

	if basis['qtype'] == 'adpt':
		from QTAG_basis import q_adapt as qprop
	elif basis['qtype'] == 'frzn':
		from QTAG_basis import frozen as qprop

	if basis['ptype'] == 'adpt':
		from QTAG_basis import p_adapt as pprop
	elif basis['ptype'] == 'frzn':
		from QTAG_basis import frozen as pprop

	if basis['atype'] == 'adpt':
		from QTAG_basis import a_adapt as aprop
	elif basis['atype'] == 'frzn':
		from QTAG_basis import frozen as aprop

	if basis['stype'] == 'adpt':
		print("Adaptable s not implemented yet! Reverting to frozen...")
		from QTAG_basis import frozen as sprop
	elif basis['stype'] == 'frzn':
		from QTAG_basis import frozen as sprop

	props=[qprop,pprop,aprop,sprop]
	return(props)

"""
.. py:function:: assemble_mom(mom)
   Returns a function object for the type of momentum 
   calculation to be performed, based on the 'mom' dict 
   specified in QTAG_config.
"""

def assemble_mom(mom):
	"""Returns a function object for the type of momentum calculation to be performed, based on the 'mom' dict specified in QTAG_config.

        Args:
                mom (dictionary): Dictionary containing keywords designating how momentum should be calculated.

        Returns:
                mom_calc (function object): Function for computing momentum.
        """

	if mom['type'] == 'raw':
		from QTAG_mom import mom_raw as mom_calc
	elif mom['type'] == 'lin_fit':
		from QTAG_mom import mom_fit as mom_calc
	elif mom['type'] == 'avg':
		from QTAG_mom import mom_avg as mom_calc
	elif mom['type'] == 'conv':
		from QTAG_mom import mom_conv as mom_calc
	else:
		sys.exit("Unrecognized keyword in mom type!")	

	return(mom_calc)

"""
.. py:function:: assemble_prop(mss)
   Returns a function object for the multi-surface scheme
   to be implemented, based on the 'mss' dict specified 
   in QTAG_config.
"""

def assemble_prop(mss):
	"""Returns a function object for the multi-surface scheme to be implemented, based on the 'mss' dict specified in QTAG_config.

        Args:
                mss (dictionary): Dictionary containing keywords designating the multi-surface scheme to be used.
        Returns:
                propagate (function object): Function for propagating basis functions.
        """

	if mss['type'] == 'sync':
		from QTAG_prop import mss_sync as propagate
	elif mss['type'] == 'fixed':
		print("Space-fixed MSS not implemented yet!")
	elif mss['type'] == 'mf':
		from QTAG_prop import mss_mf as propagate
	else:
		sys.exit("Unrecognized keyword in mss!")

	return(propagate)

"""
.. py:function:: assemble_pot(model)
   Returns function objects for the potential and surface 
   coupling, based on the 'model' dict specified in QTAG_config.
"""

def assemble_pot(model):
	"""Returns function objects for the potential and surface coupling, based on the 'model' dict specified in QTAG_config.

        Args:
                model (dictionary): Dictionary containing keywords specifying the model potential.
        Returns:
                pot (function object): Function for single-surface potential calculations.

		cplg (function object): Function for the potential coupling ground and excited states.
        """

	if model['pottype'] == 'HO':
		from QTAG_pots import model_ho_diab as pot
	elif model['pottype'] == 'T1':
        	from QTAG_pots import model_t1_diab as pot
	elif model['pottype'] == 'T2':
        	from QTAG_pots import model_t2_diab as pot
	elif model['pottype'] == 'T3':
		from QTAG_pots import model_t3_adiab as pot
	else:
		sys.exit("Unrecognized keyword in pottype!")

	if model['cplg'] == 'exact':
		from QTAG_pots import exact_gauss_cpl as cplg
	elif model['cplg'] == 'LHA':
		print("Approximating coupling potential as LHA...")
		from QTAG_ham import LHA as cplg
	elif model['cplg'] == 'none':
		from QTAG_pots import no_coupling as cplg
	else:
		sys.exit("Unrecognized keyword in cplg!")

	return(pot,cplg)
