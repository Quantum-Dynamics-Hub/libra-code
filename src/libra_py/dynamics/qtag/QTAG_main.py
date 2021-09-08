"""
..module:: QTAG_main
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
import matplotlib.pyplot as plt
import numpy as np

from QTAG_config import *

import QTAG_init
import QTAG_basis
import QTAG_ham
import QTAG_mom
import QTAG_pots
import nafh_pot
import QTAG_prop
import QTAG_calc

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
comn.check_input(univ,{"dt":0.1,"niter":1,"n_data_out":1,"n_snapshots":1},["ndof","ntraj","mass"])
comn.check_input(wf0,{},["q","p","a","s"])
comn.check_input(traj0,{"rho":1e-10},["placement","grid_dims","a0"])
comn.check_input(basis,{"qtype":"adpt","ptype":"adpt","atype":"adpt","stype":"frzn"},[])
comn.check_input(mss,{"decpl":0.9},["prop_method"])
comn.check_input(mom_params,{"adjust":"raw","beta":1e-1},[])
comn.check_input(model,{"calc_type":"LHA","coupling":"LHA"},["pot_type","rep"])
QTAG_init.qtag_checks(univ,wf0,traj0,model,model_params)

#Create potential and coupling module names from input...
model_name="model_"+str(model['pot_type']).lower()+"_"+str(model['rep']).lower()+"_nD"
coupling_name=str(model['coupling'])

#Assemble propagation types for {q,p,a,s} basis parameters...
props=QTAG_basis.param_check(basis)

#Obtain the keywords for the calculation type from the various dictionaries...
try:
	vcalc=getattr(QTAG_ham,str(model['calc_type']))
except AttributeError:
	sys.exit("Error in model dictionary: 'calc_type' keyword not recognized!")

try:
	pot=getattr(QTAG_pots,model_name)
	if coupling_name == 'LHA' or coupling_name == 'BAT':
		cplg=getattr(QTAG_ham,coupling_name)
	else:
		cplg=getattr(QTAG_pots,coupling_name)
except AttributeError:
	print("Potential Model: ",model_name)
	print("Coupling: ",coupling_name)
	sys.exit("Error in model dictionary, check that pot_type and rep are compatible!")

try:
	mom_calc=getattr(QTAG_mom,str(mom_params['adjust']))
except AttributeError:
	sys.exit("Error in mom_params dictionary: 'adjust' keyword not recognized!")

try:
	propagate=getattr(QTAG_prop,str(mss['prop_method']))
except AttributeError:
	sys.exit("Error in mss dictionary: 'prop_method' keyword not recognized!")

#Rename variables locally for convenience...
ndof,ntraj,niter,dt=univ['ndof'],univ['ntraj'],univ['niter'],univ['dt']
n_data_out,n_snapshots=univ['n_data_out'],univ['n_snapshots']
n_wf_print=int(niter/n_snapshots)

#Define MATRIX and CMATRIX objects...
b1=CMATRIX(ntraj,1);b2=CMATRIX(ntraj,1)
c1_new=CMATRIX(ntraj,1);c2_new=CMATRIX(ntraj,1)
ct_old=CMATRIX(2*ntraj,1);ct_new=CMATRIX(2*ntraj,1)
dummy=CMATRIX(ntraj,ntraj)

#Create list used in nonad vector assignments...
pop_list=[]
for i in range(ntraj):
	pop_list.append(i)

#Initialize the basis parameters {q,p,a,s} and create initial b-matrix...
try:
	qpas1=getattr(QTAG_init,str(traj0['placement']))(ndof,ntraj,traj0,wf0)
	qpas2=getattr(QTAG_init,str(traj0['placement']))(ndof,ntraj,traj0,wf0)
except AttributeError:
	sys.exit("Error in traj0 dictionary: 'placement' keyword not recognized!")

#Create initial projection vectors...
b1=QTAG_init.coeffs(ndof,ntraj,wf0,qpas1,1)
b2=QTAG_init.coeffs(ndof,ntraj,wf0,qpas2,2)
bt=QTAG_ham.nonad_assemble("cplx",ntraj,1,b1,b2,dummy)

#Begin the iterative process...
t=0.0
t_start=time.perf_counter()
for iter in range(niter):
	if iter%(niter/20)==0:
		print("Iteration = ",iter, " of ", niter)

	ov1=QTAG_calc.overlap(qpas1,qpas1)
	ov2=QTAG_calc.overlap(qpas2,qpas2)
	ov12=QTAG_calc.overlap(qpas1,qpas2)
	
	H1=QTAG_ham.hamiltonian(qpas1,ov1,1,vcalc,pot)
	H2=QTAG_ham.hamiltonian(qpas2,ov2,2,vcalc,pot)
	Hcpl=QTAG_ham.coupling(qpas1,qpas2,ov12,3,cplg,pot)

	ovt=QTAG_ham.nonad_assemble("cplx",ntraj,ntraj,ov1,ov2,dummy)
	Ht=QTAG_ham.nonad_assemble("cplx",ntraj,ntraj,H1,H2,Hcpl)

	ct_new=QTAG_ham.basis_diag(2*ntraj,Ht,ovt,bt)
	pop_submatrix(ct_new,c1_new,pop_list,[0])
	pop_submatrix(ct_new,c2_new,[ntraj+i for i in pop_list],[0])

	norm1=QTAG_calc.norm(c1_new,ov1)
	norm2=QTAG_calc.norm(c2_new,ov2)
	energy=QTAG_calc.energy(ct_new,Ht)
	
	qpas1n,qpas2n,b1,b2=propagate(mom_calc,props,qpas1,c1_new,qpas2,c2_new,norm2)

	for vals in range(4):
		qpas1[vals]=QTAG_calc.update(ndof,qpas1[vals],qpas1n[vals])
		qpas2[vals]=QTAG_calc.update(ndof,qpas2[vals],qpas2n[vals])

	bt=QTAG_ham.nonad_assemble("cplx",ntraj,1,b1,b2,dummy)
#	ct_old=QTAG_calc.update(1,ct_old,ct_new)

	norm12=np.abs((c1_new.T().conj()*ov12*c2_new).get(0))**2

	if iter%n_data_out==0:
		out1=list(qpas1[0].get(i,j) for i in range(ndof) for j in range(ntraj))
		out2=list(qpas2[0].get(i,j) for i in range(ndof) for j in range(ntraj))
		print(t,*out1,sep=' ',end='\n',file=f_traj1)
		print(t,*out2,sep=' ',end='\n',file=f_traj2)
		print(t,norm1,norm2,energy,norm12,sep=' ',end='\n',file=f_obs)

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
			QTAG_calc.wf_print_1D(qpas1,c1_new,f_wf1)
			QTAG_calc.wf_print_1D(qpas2,c2_new,f_wf2)
		elif ndof == 2:
			QTAG_calc.wf_print_2D(qpas1,c1_new,f_wf1)
			QTAG_calc.wf_print_2D(qpas2,c2_new,f_wf2)
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
		QTAG_calc.wf_print_1D(qpas1,c1_new,f_wf1)
		QTAG_calc.wf_print_1D(qpas2,c2_new,f_wf2)
	elif ndof == 2:
		QTAG_calc.wf_print_2D(qpas1,c1_new,f_wf1)
		QTAG_calc.wf_print_2D(qpas2,c2_new,f_wf2)
	else:
		print("WF printing for ndof>2 not yet implemented!")

#Close output files...
f_traj1.close();f_traj2.close();f_obs.close();
f_mom1.close(),f_mom2.close()
f_wf1.close();f_wf2.close()
f_a1.close();f_a2.close()

"""
#Arrange data to print to screen at the end of simulations...
#Everything below can be commented if no perfunctory plot output is desired.
wf1=np.loadtxt('wft1.txt')
wf2=np.loadtxt('wft2.txt')
wf1x = wf1[:,0]; wf1y = wf1[:,1]
wf2x = wf2[:,0]; wf2y = wf2[:,1]

gbc1=np.loadtxt('gbc1.txt')
gbc2=np.loadtxt('gbc2.txt')
#gbc1x=gbc1[:,0]
obsrv=np.loadtxt('obs.txt')
obx = obsrv[:,0]; nrm1 = obsrv[:,1]; nrm2 = obsrv[:,2]; en = obsrv[:,3]

fig,axs=plt.subplots(2,2)
axs[0,0].plot(wf1x,wf1y,label='S$_0$')
axs[0,0].plot(wf2x,wf2y,label='S$_1$')
legend= axs[0,0].legend(loc='upper left')

for i in range(1,ntraj+1):
	axs[0,1].plot(gbc1[:,0],gbc1[:,i],'tab:blue',label='S$_0$')
	axs[0,1].plot(gbc2[:,0],gbc2[:,i],'tab:orange',label='S$_1$')

axs[1,0].plot(obx,nrm1)
axs[1,0].plot(obx,nrm2)
axs[1,1].plot(obx,(en-en[0])/en[0]*100,'tab:green')

plt.show()
"""
