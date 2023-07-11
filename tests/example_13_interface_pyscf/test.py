from pyscf import gto, scf
from pyscf import ci


#mol = gto.M(atom="H 0 0 0; H 0 0 1.5", basis="631g")
#mol.atom = [['O',(0, 0, 0)], ['H',(0, 1, 0)], ['H',(0, 0, 1)]]
mol = gto.M(atom= [['H',(0, 0, 0)], ['H',(0, 0, 1)]], basis="631g")

mf = scf.RHF(mol)
mf.verbose=-1
res_mf = mf.kernel()
print "HF energy = ", res_mf
force = mf.nuc_grad_method().kernel()
print force

myci = ci.cisd.CISD(mf)
myci.verbose = -1
myci.nstates = 5
res_ci = myci.kernel()
print "CISD correlation energy = ", res_ci[0]
print "CISD vectors = ", res_ci[1]

for i in xrange(5):
    force = myci.nuc_grad_method().kernel(state=i)
    print force




#dm1 = myci.make_rdm1()

#print dm1

#ci = mol.apply("CISD").run()
#force = ci.nuc_grad_method().kernel()
#print force


