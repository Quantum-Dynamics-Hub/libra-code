#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
import cmath
import math
import os
import sys
import unittest


cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../_build/src/math_linalg")
sys.path.insert(1,cwd+"/../_build/src/math_random")
sys.path.insert(1,cwd+"/../_build/src/ivr")


# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cyglinalg import *
    from cygrandom import *
    from cygivr import *


elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from liblinalg import *
    from librandom import *
    from libivr import *



rnd = Random()


class Test_IVR(unittest.TestCase):
    """ Summary of the tests:
    """


    def test_1(self):
        """Run Husimi - 1 point """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()
        print "q,p = "; ivr_Husimi(q, p, W, rnd).show_matrix()


    def test_2(self):
        """Run Husimi - N points """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()

        res = ivr_Husimi(q, p, W, rnd, 10)
        print "q,p = "; res[0].show_matrix(); res[9].show_matrix();



    def test_3(self):
        """Run DHK - 1 point """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()
        print "q,p = "; ivr_DHK(q, p, W, rnd).show_matrix()


    def test_4(self):
        """Run DHK - N points """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()

        res = ivr_DHK(q, p, W, rnd, 10)
        print "q,p = "; res[0].show_matrix(); res[9].show_matrix();



    def test_5(self):
        """Run FB_MQC - 1 point """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        cQ = MATRIX(3,3);  tunq = 0.1;
        cQ.set(0,0, tunq);  cQ.set(1,1, tunq);  cQ.set(2,2, tunq);

        cP = MATRIX(3,3);  tunp = 0.1;
        cP.set(0,0, tunp);  cP.set(1,1, tunp);  cP.set(2,2, tunp);
                
        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()
        for flag in [0, 1, 2]:
            print "flag = ", flag
            res = ivr_FB_MQC(q, p, W, cQ, cP, flag, rnd)
            print "q,p = "; res.show_matrix()


    def test_6(self):
        """Run FB_MQC - N points """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        cQ = MATRIX(3,3);  tunq = 0.1;
        cQ.set(0,0, tunq);  cQ.set(1,1, tunq);  cQ.set(2,2, tunq);

        cP = MATRIX(3,3);  tunp = 0.1;
        cP.set(0,0, tunp);  cP.set(1,1, tunp);  cP.set(2,2, tunp);
                
        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()
        for flag in [0, 1, 2]:
            print "flag = ", flag
            res = ivr_FB_MQC(q, p, W, cQ, cP, flag, rnd, 10)
            print "q,p = "; res[0].show_matrix(); res[9].show_matrix()



    def test_7(self):
        """Run FF_MQC - 1 point """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        cQ = MATRIX(3,3);  tunq = 0.1;
        cQ.set(0,0, tunq);  cQ.set(1,1, tunq);  cQ.set(2,2, tunq);

        cP = MATRIX(3,3);  tunp = 0.1;
        cP.set(0,0, tunp);  cP.set(1,1, tunp);  cP.set(2,2, tunp);
                
        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()
        res = ivr_FF_MQC(q, p, W, cQ, cP, rnd)
        print "q,p = "; res.show_matrix()


    def test_8(self):
        """Run FF_MQC - N points """

        q,p = MATRIX(3,1), MATRIX(3,1)         
        q.set(0, 0.0);     p.set(0, 0.0)
        q.set(1, 2.0);     p.set(1, -1.0)
        q.set(2, -1.0);    p.set(2, 1.0)

        W = MATRIX(3,3)
        W.set(0,0, 0.1);  W.set(1,1, 0.1);  W.set(2,2, 0.1);

        cQ = MATRIX(3,3);  tunq = 0.1;
        cQ.set(0,0, tunq);  cQ.set(1,1, tunq);  cQ.set(2,2, tunq);

        cP = MATRIX(3,3);  tunp = 0.1;
        cP.set(0,0, tunp);  cP.set(1,1, tunp);  cP.set(2,2, tunp);
                
        print "qIn = "; q.show_matrix()
        print "pIn = "; p.show_matrix()
        res = ivr_FF_MQC(q, p, W, cQ, cP, rnd, 10)
        print "q,p = "; res[0].show_matrix(); res[9].show_matrix()




#        self.assertAlmostEqual( F[2].x, -10.0 );   # 0 + 0 - 10

        

               


if __name__=='__main__':
    unittest.main()




