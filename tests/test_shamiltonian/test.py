#import pdb
#import sys
import torch
#pdb.set_trace()
print(torch.__version__)

#sys.path.append("../../_build_pytorch/src/") 
#sys.path.append("../../_build_pytorch/src/shamiltonian/")

#import liblibra_core
from liblibra_core import *
#from libshamiltonian import *

#a = MATRIX(2,2)
s = sHamiltonian(2, 2, 2)

#print( s.ovlp_dia) 
#print(t2p(s.ovlp_dia))

x = torch.ones(2,2,2)
#x = torch.tensor([1.0, 2.0, 3.0])
#y = process_tensor(x)

#print(y)
#print(x)

s.set_tensor(x)

print( s.ovlp_dia)
print(t2p(s.ovlp_dia))

x[0,0,0] = -1
print(x)
print(t2p(s.ovlp_dia))

#s.get_tensor(x)

#y = t2p(x)
#print(y)

#ss = s.get_ovlp_dia()

#print(ss)

##s.ovlp_dia[0,0,0] = 1.0


#h = nHamiltonian(1, 1, 2)


