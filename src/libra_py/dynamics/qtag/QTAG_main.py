"""
..module:: QTAG_main
  :platform: Unix, Windows
  :synopsis: This module is the main body of the QTAG method, responsible for running cycles of the BOT algorithm.
..moduleauthors :: Matthew Dutra
"""

import sys
import cmath
import math
import os
from liblibra_core import *
import util.libutil as comn
import time

from libra_py import data_outs
import matplotlib.pyplot as plt
import numpy as np

from QTAG_config import *
from QTAG_assembler import *
from QTAG_init import coeffs
from QTAG_ham import *

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

#Assemble the appropriate functions from the input file QTAG_config using QTAG_assembler...
initialize=assemble_init(traj0)
pot,cplg=assemble_pot(model)
if model['pottype'] == 'HO':
	model_params['d2']=model_params['d2']*np.sqrt(model_params['k1'])
mom_calc=assemble_mom(mom_params)
propagate=assemble_prop(mss)
props=assemble_basis(basis)

#Rename variables locally for convenience...
ntraj,niter,dt=univ['ntraj'],univ['niter'],univ['dt']
nout=univ['nout']

b1=CMATRIX(ntraj,1);b2=CMATRIX(ntraj,1)
c1_new=CMATRIX(ntraj,1);c2_new=CMATRIX(ntraj,1)
ct_old=CMATRIX(2*ntraj,1);ct_new=CMATRIX(2*ntraj,1)
dummy=CMATRIX(ntraj,ntraj)

pop_list=[]
for i in range(ntraj):
	pop_list.append(i)

#Initialize the basis parameters {q,p,a,s} and create initial b-matrix...
qpas1=initialize(ntraj,traj0,wf0)
qpas2=initialize(ntraj,traj0,wf0)

b1=coeffs(ntraj,wf0,qpas1,1)
b2=coeffs(ntraj,wf0,qpas2,2)
bt=nonad_assemble("cplx",ntraj,1,b1,b2,dummy)

#Begin the iterative process...
t=0.0
t_start=time.perf_counter()
for iter in range(niter):
	if iter%(niter/20)==0:
		print("Iteration = ",iter, " of ", niter)

	ov1=QTAG_calc.overlap(qpas1,qpas1)
	ov2=QTAG_calc.overlap(qpas2,qpas2)
	ov12=QTAG_calc.overlap(qpas1,qpas2)

	H1=hamiltonian(qpas1,ov1,1,pot)
	H2=hamiltonian(qpas2,ov2,2,pot)
	Hcpl=coupling(qpas1,qpas2,ov12,3,cplg,pot)

	ovt=nonad_assemble("cplx",ntraj,ntraj,ov1,ov2,dummy)
	Ht=nonad_assemble("cplx",ntraj,ntraj,H1,H2,Hcpl)

	ct_new=basis_diag(2*ntraj,Ht,ovt,bt)
	pop_submatrix(ct_new,c1_new,pop_list,[0])
	pop_submatrix(ct_new,c2_new,[ntraj+i for i in pop_list],[0])

	norm1=QTAG_calc.norm(c1_new,ov1).real
	norm2=QTAG_calc.norm(c2_new,ov2).real
	energy=QTAG_calc.energy(ct_new,Ht).real

	qpas1n,qpas2n,b1,b2=propagate(mom_calc,props,qpas1,c1_new,qpas2,c2_new,norm2)

	qpas1=QTAG_calc.update(4,qpas1,qpas1n)
	qpas2=QTAG_calc.update(4,qpas2,qpas2n)

	bt=nonad_assemble("cplx",ntraj,1,b1,b2,dummy)
	ct_old=QTAG_calc.update(1,ct_old,ct_new)

	norm12=np.abs((c1_new.T().conj()*ov12*c2_new).get(0))**2

	if iter%nout==0:
		out1=list(qpas1.get(i,0) for i in range(ntraj))
		out2=list(qpas2.get(i,0) for i in range(ntraj))
		print(t,*out1,sep=' ',end='\n',file=f_traj1)
		print(t,*out2,sep=' ',end='\n',file=f_traj2)
		print(t,norm1,norm2,energy,norm12,sep=' ',end='\n',file=f_obs)

		out1=list(qpas1.get(i,1) for i in range(ntraj))
		out2=list(qpas2.get(i,1) for i in range(ntraj))
		print(t,*out1,sep=' ',end='\n',file=f_mom1)
		print(t,*out2,sep=' ',end='\n',file=f_mom2)

		out1=list(qpas1.get(i,2) for i in range(ntraj))
		out2=list(qpas2.get(i,2) for i in range(ntraj))
		print(t,*out1,sep=' ',end='\n',file=f_a1)
		print(t,*out2,sep=' ',end='\n',file=f_a2)

#	if iter%10==1:
#		QTAG_calc.wf_print(qpas1,c1_new,f_wf1)
#		QTAG_calc.wf_print(qpas2,c2_new,f_wf2)
	t+=dt

t_end=time.perf_counter()
print("Avg iter time =", (t_end-t_start)/niter, "seconds")
print("Total time =", (t_end-t_start), "seconds")

#Print final WF output on each surface...
QTAG_calc.wf_print(qpas1,c1_new,f_wf1)
QTAG_calc.wf_print(qpas2,c2_new,f_wf2)

#Close output files...
f_traj1.close();f_traj2.close();f_obs.close();
f_mom1.close(),f_mom2.close()
f_wf1.close();f_wf2.close()
f_a1.close();f_a2.close()

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
