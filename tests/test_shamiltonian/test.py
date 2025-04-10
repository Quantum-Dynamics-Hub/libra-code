#import pdb
import sys
import torch
#pdb.set_trace()
print(torch.__version__)

#sys.path.append("../../_build_pytorch/src/") 
#sys.path.append("../../_build_pytorch/src/shamiltonian/")

#import liblibra_core
from liblibra_core import *
#from libshamiltonian import *

import torchHolstein

# Constructor with the default initialization
print("Constructor of the sHamiltonian object...")
s = sHamiltonian(2, 2, 2)
print( s.ovlp_dia) 
print( cpp2py(s.ovlp_dia) )

# Create Py tensor
print("Creating external Py tensor")
x = torch.ones(2,2,2)
address = x.data_ptr()
print(f"Memory address: {hex(address)}")
print(x)


# Bind the Py tensor with the object
print("sHamiltonian member after it is binded to the PyTorch object")
s.bind("ovlp_dia", x)
print(s.ovlp_dia)
print( cpp2py(s.ovlp_dia) )


print("Creating a copy sHamiltonian")
s_copy = sHamiltonian(s)
print(s_copy.ovlp_dia)
print( cpp2py(s_copy.ovlp_dia) )


# Change x in the Python side
print("Change the Python Tensor, let's see how the internal one changed")
x[0,0,1] = -1
print(x)
print(s.ovlp_dia)
print(cpp2py(s.ovlp_dia))

print("What about the copied object?")
print(s_copy.ovlp_dia)
print( cpp2py(s_copy.ovlp_dia) )


print("Showing transpose in place")
x.transpose_(1,2)
print(s.ovlp_dia)
print(cpp2py(s.ovlp_dia))


# Add calling external Python function:

def quadratic(q, params):
    return q*q;

y = torch.randn(2,2,2, dtype=complex)
print(y)

s.compute("ham_dia", quadratic,  y, {})

ham_dia = cpp2py(s.ham_dia)
print("Diabatic Ham")
print(ham_dia)

# Do the diagonalization
s.dia2adi()


#sys.exit()

#ham_adi = torch.diag_embed( cpp2py(s.ham_adi) )
ham_adi = cpp2py( s.ham_adi ) 
print("Adiabatic Ham")
print(ham_adi)
U = cpp2py(s.basis_transform)

print("Eigenvectors")
print(U)

#sys.exit()

print("Checking:  H_dia * U = U * H_adi")
print("LHS", ham_dia @ U  )
print("RHS", U @ ham_adi)

#print("LHS", ham_dia @ U.transpose(-1,-2)  )
#print("RHS", U.transpose(-1,-2) @ ham_adi)

#sys.exit()


s.compute_nacs_and_grads()
print("NACs")
print( cpp2py(s.dc1_adi) )
print( cpp2py(s.d1ham_adi))

print("Adiabatic forces")
f_adi = s.forces_adi()
#print(cpp2py(f_adi) )


#========
print("Now testing torch Holstein")
q = torch.randn(2, 1, dtype=float)  # nbeads, nnucl
#Holstein2(q, {})

s.compute( torchHolstein.Holstein2,  q, {"E_n":[0.0, 0.001], "x_n":[0.0, 1.0], "k_n":[0.001, 0.001]})

ham_dia = cpp2py(s.ham_dia)
print("Diabatic Ham")
print(ham_dia)


#sys.exit()

# Do the diagonalization
s.dia2adi()


#ham_adi = torch.diag_embed( cpp2py(s.ham_adi) )
ham_adi = cpp2py( s.ham_adi )
print("Adiabatic Ham")
print(ham_adi)
U = cpp2py(s.basis_transform)

print("Eigenvectors")
print(U)

print("Checking:  H_dia * U = U * H_adi")
print("LHS", ham_dia @ U)
print("RHS", U @ ham_adi)


s.compute_nacs_and_grads()
print("NACs")
print( cpp2py(s.dc1_adi) )
print( cpp2py(s.d1ham_adi))

print("Adiabatic forces")
f_adi = s.forces_adi()
#print(cpp2py(f_adi) )


d1ham_dia = cpp2py(s.d1ham_dia)
print("Diabatic Ham, derivatives")
print(d1ham_dia)
















#s.get_tensor(x)

#y = t2p(x)
#print(y)

#ss = s.get_ovlp_dia()

#print(ss)

##s.ovlp_dia[0,0,0] = 1.0


#h = nHamiltonian(1, 1, 2)


