import unittest
import numpy as np
from liblibra_core import *
from libra_py.workflows.nbra.step3 import make_active_space




class Test_make_active_space(unittest.TestCase):

    def test_full_active_pace(self):
        """Test full active space"""
        
        
        #sample_matrix = np.random.rand(8,8)
        # Testing the full active space
        num_of_occ_orbital = 2
        num_of_unocc_orbital = 2
        ks_homo_index = 2
        data_dim = 8
        # The new active space and the KS HOMO index in that active space which again starts from 1
        new_active_space, new_ks_homo_index = make_active_space(num_of_occ_orbital, num_of_unocc_orbital, data_dim, ks_homo_index)
        self.assertEqual(list(range(8)), new_active_space)
        self.assertEqual(ks_homo_index, new_ks_homo_index)
        #print('The new active space is:', new_active_space)
        #print('The new KS HOMO index is:', new_ks_homo_index)
        #print('The initial matrix is:\n',sample_matrix)
        #print('The matrix with the new active space:\n',sample_matrix[new_active_space,:][:,new_active_space])

    def test_partial_active_space(self):
        """Test partial active space"""


        #sample_matrix = np.random.rand(10,10)
        # Testing the full active space
        num_of_occ_orbital = 3
        num_of_unocc_orbital = 1
        ks_homo_index = 3
        data_dim = 10
        # The new active space and the KS HOMO index in that active space which again starts from 1
        new_active_space, new_ks_homo_index = make_active_space(num_of_occ_orbital, num_of_unocc_orbital, data_dim, ks_homo_index)
        self.assertEqual(___, new_active_space)
        self.assertEqual(___, new_ks_homo_index)
       # print('The new active space is:', new_active_space)
       # print('The new KS HOMO index is:', new_ks_homo_index)
       # print('The initial matrix is:\n',sample_matrix)
       # print('The matrix with the new active space:\n',sample_matrix[new_active_space,:][:,new_active_space])



if __name__=='__main__':
    unittest.main()

