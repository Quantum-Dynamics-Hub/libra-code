import unittest
import numpy as np
from liblibra_core import *
from libra_py.workflows.nbra.step3 import make_active_space




class Test_make_active_space(unittest.TestCase):

    def test_full_active_pace(self):
        """Test full active space"""

        full_active_space = list(range(8))
        # Testing the full active space
        num_of_occ_orbital = 2
        num_of_unocc_orbital = 2
        ks_homo_index = 2
        data_dim = 8
        # The new active space and the KS HOMO index in that active space which again starts from 1
        new_active_space, new_ks_homo_index = make_active_space(num_of_occ_orbital, num_of_unocc_orbital, data_dim, ks_homo_index)
        self.assertEqual(full_active_space, new_active_space)
        self.assertEqual(ks_homo_index, new_ks_homo_index)

    def test_partial_active_space(self):
        """Test partial active space"""

        partial_active_space = [0, 1, 2, 3, 5, 6, 7, 8]
        partial_active_space_homo_index = 3
        # Testing the full active space
        num_of_occ_orbital = 3
        num_of_unocc_orbital = 1
        ks_homo_index = 3
        data_dim = 10
        # The new active space and the KS HOMO index in that active space which again starts from 1
        new_active_space, new_ks_homo_index = make_active_space(num_of_occ_orbital, num_of_unocc_orbital, data_dim, ks_homo_index)
        self.assertEqual(partial_active_space, new_active_space)
        self.assertEqual(partial_active_space_homo_index, new_ks_homo_index)


if __name__=='__main__':
    unittest.main()

