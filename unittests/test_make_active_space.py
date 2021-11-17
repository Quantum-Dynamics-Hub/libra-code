from liblibra_core import *
from libra_py.workflows.nbra.step3 import make_active_space

print('You are in testing the make_active_space function test.')
num_of_occ_orbital = input('number of occupied orbitals in the new active space:')
num_of_unocc_orbital = input('number of unoccupied orbitals in the new active space:')
ks_homo_index = input('the Kohn-Sham HOMO index in the raw data files:')
data_dim = input('data dimension of the raw data files:')
# The new active space and the KS HOMO index in that active space which again starts from 1
new_active_space, new_ks_homo_index = make_active_space(int(num_of_occ_orbital), int(num_of_unocc_orbital), int(data_dim), int(ks_homo_index))
print('The new active space is:', new_active_space)
print('The new KS HOMO index is:', new_ks_homo_index)

