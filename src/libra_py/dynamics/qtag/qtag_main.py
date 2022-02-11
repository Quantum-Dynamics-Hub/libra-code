"""
..module:: qtag_main
  :platform: Unix, Windows
  :synopsis: This module is the main body of the QTAG method, responsible for running cycles of the BOT algorithm.
..moduleauthors :: Matthew Dutra
"""

import sys
import os
from liblibra_core import *
import util.libutil as comn
import time
import datetime

from libra_py import data_outs
import numpy as np

import qtag_init
import qtag_basis
import qtag_ham
import qtag_mom
import qtag_pots
import qtag_prop
import qtag_calc

"""
univ = {"ndof" : 1, "ntraj" : 35, "dt" : 0.01, "niter" : 800, "mass" : [1.0], "n_snapshots" : 4, "n_data_out" : 1}
wf0 = {"q" : [-2.0], "p" : [0.0], "a" : [1.0], "s" : [0.0]}
traj0 = {"placement" : "grid", "grid_dims" : [35], "rho" : 1e-12, "a0" : [18.0]}
basis = {"qtype" : "adpt", "ptype" : "adpt", "atype" : "frzn", "stype" : "frzn"}
mss = {"prop_method" : "mean_field", "decpl" : 0.3, "mirror" : "true"}
mom_params = {"adjust" : "average", "beta" : 1e-3}
model = {"pot_type" : "HO", "rep" : "diab", "calc_type" : "BAT", "coupling" : "exact"}
model_params = {"rNaF" : 3.779, "k1" : [10.0], "k2" : [10.0], "x0" : [1.0], "d1" : [1.0], "d2" : [0.5], "d3" : [2.0]}
"""

def run_qtag(univ,wf0,traj0,basis,mss,mom_params,model,model_params):
#Create output file objects...
    f_traj1=open('gbc1.txt', 'w')
    f_traj2=open('gbc2.txt', 'w')
    f_mom1=open('p1.txt', 'w')
    f_mom2=open('p2.txt', 'w')
    f_a1=open('a1.txt', 'w')
    f_a2=open('a2.txt', 'w')
    f_obs=open('obs.txt', 'w')
    f_wf1=open('wft1.txt', 'w')
    f_wf2=open('wft2.txt', 'w')

#Libra + QTAG checks for input dictionaries...
    comn.check_input(univ,{"ndof":1,"ntraj":1,"mass":[1836.0],"dt":0.1,"niter":1,"n_data_out":1,"n_snapshots":1},[])
    comn.check_input(wf0,{},["q","p","a","s"])
    comn.check_input(traj0,{"rho":1e-10},["placement","grid_dims","a0"])
    comn.check_input(basis,{"qtype":"adpt","ptype":"adpt","atype":"adpt","stype":"frzn"},[])
    comn.check_input(mss,{"decpl":0.9},["prop_method"])
    comn.check_input(mom_params,{"adjust":"unmodified","beta":1e-1},[])
    comn.check_input(model,{"calc_type":"LHA","coupling":"LHA"},["pot_type","rep"])
    qtag_init.qtag_checks(univ,wf0,traj0,mss,model,model_params)

#Create potential and coupling module names from input...
    model_name="model_"+str(model['pot_type']).lower()+"_"+str(model['rep']).lower()+"_nD"
    coupling_name=str(model['coupling'])

#Assemble propagation types for {q,p,a,s} basis parameters...
    props=qtag_basis.param_check(qtag_params[basis],qtag_params[mss]['prop_method'])

#Obtain the keywords for the calculation type from the various dictionaries...
    try:
        vcalc=getattr(qtag_ham,str(model['calc_type']))
    except AttributeError:
        sys.exit("Error in model dictionary: 'calc_type' keyword not recognized!")

    try:
        pot=getattr(qtag_pots,model_name)
        if coupling_name == 'LHA' or coupling_name == 'BAT':
            cplg=getattr(qtag_ham,coupling_name)
        else:
            cplg=getattr(qtag_pots,coupling_name)
        univ['pot_fxn']=pot;univ['cplg_fxn']=cplg
    except AttributeError:
        print("Potential Model: ",model_name)
        print("Coupling: ",coupling_name)
        sys.exit("Error in model dictionary, check that pot_type and rep are compatible!")

    try:
        mom_calc=getattr(qtag_mom,str(mom_params['adjust']))
    except AttributeError:
        sys.exit("Error in mom_params dictionary: 'adjust' keyword not recognized!")

    try:
        propagate=getattr(qtag_prop,str(mss['prop_method']))
    except AttributeError:
        sys.exit("Error in mss dictionary: 'prop_method' keyword not recognized!")

#Rename variables locally for convenience...
    ndof,ntraj,niter,dt=univ['ndof'],univ['ntraj'],univ['niter'],univ['dt']
    n_data_out,n_snapshots=univ['n_data_out'],univ['n_snapshots']
    n_wf_print=int(niter/n_snapshots)
    if n_wf_print==0:
        n_wf_print=1

#Define MATRIX and CMATRIX objects...
    b1=CMATRIX(ntraj,1);b2=CMATRIX(ntraj,1)
    ct_old=CMATRIX(ntraj,1)
    ct_new=CMATRIX(ntraj,1)

#Create list used in nonad vector assignments...
    pop_list=[]
    for i in range(ntraj):
        pop_list.append(i)

#Initialize the basis parameters {q,p,a,s} and create initial b-matrix...
    try:
        for n in range nstates:
            qpas=getattr(qtag_init,str(traj0['placement']))(ndof,ntraj,traj0,wf0)
            qpas2=getattr(qtag_init,str(traj0['placement']))(ndof,ntraj,traj0,wf0)
    except AttributeError:
        sys.exit("Error in traj0 dictionary: 'placement' keyword not recognized!")

#Create initial projection vectors...
    b1=qtag_init.coeffs(ndof,ntraj,wf0,qpas1,1)
    b2=qtag_init.coeffs(ndof,ntraj,wf0,qpas2,2)
    nsurf1=[1 for i in range(ntraj)]
    nsurf2=[2 for i in range(ntraj)]

    bt=qtag_calc.nonad_assemble("cplx",ntraj,1,b1,b2,dummy)

#Begin the iterative process...
	t=0.0
	t_start=time.perf_counter()
	for iter in range(niter):
		if iter%100*n_data_out==0:
			print("Iteration = ",iter, " of ", niter)

		ov1=qtag_calc.overlap(ntraj,qpas1,qpas1)
		ov2=qtag_calc.overlap(ntraj,qpas2,qpas2)
		ov12=qtag_calc.overlap(ntraj,qpas1,qpas2)

		H1=qtag_ham.hamiltonian(univ,qpas1,ov1,1,vcalc,pot,model_params)
		H2=qtag_ham.hamiltonian(univ,qpas2,ov2,2,vcalc,pot,model_params)
		Hcpl=qtag_ham.coupling(univ,qpas1,qpas2,ov12,3,cplg,pot,model_params)

		if mss['prop_method'] == 'two_surf':
			ovt=qtag_calc.nonad_assemble("cplx",ntraj,ntraj,ov1,ov2,dummy)
		else:
			ovt=qtag_calc.nonad_assemble("cplx",ntraj,ntraj,ov1,ov2,dummy)

#		res=FullPivLU_rank_invertible(ovt)
#		print(F" matrix ovt: rank = {res[0]}  is invertible = {res[1]}")

		Ht=qtag_calc.nonad_assemble("cplx",ntraj,ntraj,H1,H2,Hcpl)

		ct_new=qtag_ham.basis_diag(2*ntraj,dt,Ht,ovt,bt)
		pop_submatrix(ct_new,c1_new,pop_list,[0])
		pop_submatrix(ct_new,c2_new,[ntraj+i for i in pop_list],[0])

		norm1=qtag_calc.norm(c1_new,ov1)
		norm2=qtag_calc.norm(c2_new,ov2)
		energy=qtag_calc.energy(ct_new,Ht)
		if iter%100*n_data_out==0:
			print("SURFACE 1 NORM= %8.6f, SURFACE 2 NORM= %8.6f, TOTAL NORM= %8.6f" %(norm1,norm2,norm1+norm2))
			print("SYSTEM ENERGY= ", energy)

		qpas1n,qpas2n,b1,b2=propagate(univ,mss,mom_calc,props,model_params,qpas1,c1_new,qpas2,c2_new,norm2,mom_params['beta'])

		for vals in range(4):
			qpas1[vals]=MATRIX(qpas1n[vals])
			qpas2[vals]=MATRIX(qpas2n[vals])

		bt=qtag_calc.nonad_assemble("cplx",ntraj,1,b1,b2,dummy)

		if iter%n_data_out==0:
			out1=list(qpas1[0].get(i,j) for i in range(ndof) for j in range(ntraj))
			out2=list(qpas2[0].get(i,j) for i in range(ndof) for j in range(ntraj))
			print(t,*out1,sep=' ',end='\n',file=f_traj1)
			print(t,*out2,sep=' ',end='\n',file=f_traj2)
			print(t,norm1,norm2,energy,sep=' ',end='\n',file=f_obs)

			out1=list(qpas1[1].get(i,j) for i in range(ndof) for j in range(ntraj))
			out2=list(qpas2[1].get(i,j) for i in range(ndof) for j in range(ntraj))
			print(t,*out1,sep=' ',end='\n',file=f_mom1)
			print(t,*out2,sep=' ',end='\n',file=f_mom2)

			out1=list(qpas1[2].get(i,j) for i in range(ndof) for j in range(ntraj))
			out2=list(qpas2[2].get(i,j) for i in range(ndof) for j in range(ntraj))
			print(t,*out1,sep=' ',end='\n',file=f_a1)
			print(t,*out2,sep=' ',end='\n',file=f_a2)

		if iter%n_wf_print==0:
			if ndof == 1:
				qtag_calc.wf_print_1D(ntraj,qpas1,c1_new,f_wf1)
				qtag_calc.wf_print_1D(ntraj,qpas2,c2_new,f_wf2)
			elif ndof == 2:
				qtag_calc.wf_print_2D(ntraj,qpas1,c1_new,f_wf1)
				qtag_calc.wf_print_2D(ntraj,qpas2,c2_new,f_wf2)
			else:
				print("WF printing for ndof>2 not yet implemented!")

		t+=dt

	t_end=time.perf_counter()
	print("Avg iter time =", (t_end-t_start)/niter, "seconds")
	print("Total time =", (t_end-t_start), "seconds")
	print("Total time =", datetime.timedelta(seconds=(t_end-t_start)))

#Print final WF output on each surface...
	if niter%n_wf_print != 0:
		if ndof == 1:
			qtag_calc.wf_print_1D(ntraj,qpas1,c1_new,f_wf1)
			qtag_calc.wf_print_1D(ntraj,qpas2,c2_new,f_wf2)
		elif ndof == 2:
			qtag_calc.wf_print_2D(ntraj,qpas1,c1_new,f_wf1)
			qtag_calc.wf_print_2D(ntraj,qpas2,c2_new,f_wf2)
		else:
			print("WF printing for ndof>2 not yet implemented!")

#Close output files...
	f_traj1.close();f_traj2.close();f_obs.close();
	f_mom1.close(),f_mom2.close()
	f_wf1.close();f_wf2.close()
	f_a1.close();f_a2.close()

#run_qtag(univ,wf0,traj0,basis,mss,mom_params,model,model_params)
