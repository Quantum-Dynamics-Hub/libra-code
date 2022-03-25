import unittest
import numpy as np
from liblibra_core import *
from libra_py.workflows.nbra.step3 import make_active_space




class Test_make_active_space(unittest.TestCase):

    def test_full_active_pace(self):
        """Test full active space"""

        # Testing the full active space
        # The new active space and the KS HOMO index in that active space which again starts from 1
        new_active_space, new_ks_homo_index = make_active_space(2, 2, 8, 2)
        self.assertEqual(new_active_space, list(range(8)) )
        self.assertEqual(new_ks_homo_index, 2)

    def test_partial_active_space(self):
        """Test partial active space"""

        # Testing the partial active space
        # The new active space and the KS HOMO index in that active space which again starts from 1
        new_active_space, new_ks_homo_index = make_active_space(3, 1, 10, 3)
        self.assertEqual(new_active_space, [0, 1, 2, 3, 5, 6, 7, 8] )
        self.assertEqual(new_ks_homo_index, 3)


if __name__=='__main__':
    unittest.main()

