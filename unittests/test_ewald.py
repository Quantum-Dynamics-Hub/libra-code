#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
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
sys.path.insert(1,cwd+"/../_build/src/pot")
sys.path.insert(1,cwd+"/../_build/src/converters")
sys.path.insert(1,cwd+"/../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cygpot import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libpot import *
    from liblinalg import *







class Test_Ewald(unittest.TestCase):
    """ Summary of the tests:
    """

    def test_1(self):
        """Computing the electrostatic sum for NaCl system
           with 4 NaCl molecules in the unit cell (8 atoms)
           The refrence energy value is close to the one reported 
           in the paper:
           Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for 
           Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327 
        """

        Angst = 1.889725989       # 1 Angstrom in atomic units
        hartree = 627.5094709     # 1 Ha = 627.5.. kcal/mol
        
        R = VECTORList()
        Q = doubleList()
        F = VECTORList()
        
        a = 5.63 * Angst   # in a.u.
        etha = 2.5 * Angst # in a.u.
        
        tv1 = a*VECTOR(1.0, 0.0, 0.0)
        tv2 = a*VECTOR(0.0, 1.0, 0.0)
        tv3 = a*VECTOR(0.0, 0.0, 1.0)
        
        box = MATRIX3x3(tv1, tv2, tv3)
        stress = MATRIX3x3()
        
        R.append( a*VECTOR(0.0, 0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na
        R.append( a*VECTOR(0.5, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na
        R.append( a*VECTOR(0.5, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na
        R.append( a*VECTOR(0.0, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na

        R.append( a*VECTOR(0.5, 0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
        R.append( a*VECTOR(0.0, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
        R.append( a*VECTOR(0.0, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
        R.append( a*VECTOR(0.5, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
                                                              
        
        pbc_deg = 2
        rec_deg = 1
        R_on = 10 * Angst
        R_off = 12 * Angst
        epsilon = 1.0  # dielectric constant

        # We'll get energies in kcal/mol        
        energy = Elec_Ewald3D(R, Q, box, epsilon, F, stress, rec_deg, pbc_deg, etha, R_on, R_off ) * hartree

        self.assertAlmostEqual( energy, -824.8591 , 4);


        # The total sum should be zero, no matter what
        Ftot = VECTOR(0.0, 0.0, 0.0)        
        for f in F:
            Ftot += f
            #print f.x, f.y, f.z
        
        self.assertAlmostEqual( Ftot.x, 0.0);
        self.assertAlmostEqual( Ftot.y, 0.0);
        self.assertAlmostEqual( Ftot.z, 0.0);

        # Stress is close to zero
        self.assertAlmostEqual( stress.xy, 0.0);
        self.assertAlmostEqual( stress.xz, 0.0);
        self.assertAlmostEqual( stress.yx, 0.0);
        self.assertAlmostEqual( stress.yz, 0.0);
        self.assertAlmostEqual( stress.zx, 0.0);
        self.assertAlmostEqual( stress.zy, 0.0);

        self.assertAlmostEqual( stress.xx, stress.yy);
        self.assertAlmostEqual( stress.yy, stress.zz);
       


    def test_2(self):
        """Computing the electrostatic sum for NaCl system
           with 4 NaCl molecules in the unit cell (8 atoms)
           The refrence energy value is close to the one reported 
           in the paper:
           Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for 
           Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327 

           Here, we distort the lattice to induce some forces and stress
        """

        Angst = 1.889725989       # 1 Angstrom in atomic units
        hartree = 627.5094709     # 1 Ha = 627.5.. kcal/mol
        
        R = VECTORList()
        Q = doubleList()
        F = VECTORList()
        
        a = 5.63 * Angst   # in a.u.
        etha = 2.5 * Angst # in a.u.
        
        tv1 = a*VECTOR(1.0, 0.0, 0.0)
        tv2 = a*VECTOR(0.0, 1.0, 0.0)
        tv3 = a*VECTOR(0.0, 0.0, 1.0)
        
        box = MATRIX3x3(tv1, tv2, tv3)
        stress = MATRIX3x3()
        
        R.append( a*VECTOR(0.0, 0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na
        R.append( a*VECTOR(0.5, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na
        R.append( a*VECTOR(0.5, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na
        R.append( a*VECTOR(0.0, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( 1.0)   # Na

        R.append( a*VECTOR(0.55,0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
        R.append( a*VECTOR(0.0, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
        R.append( a*VECTOR(0.0, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
        R.append( a*VECTOR(0.5, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append(-1.0)   # Cl
                                                              
        
        pbc_deg = 2
        rec_deg = 1
        R_on = 10 * Angst
        R_off = 12 * Angst
        epsilon = 1.0  # dielectric constant

        # We'll get energies in kcal/mol        
        energy = Elec_Ewald3D(R, Q, box, epsilon, F, stress, rec_deg, pbc_deg, etha, R_on, R_off ) * hartree

#        self.assertAlmostEqual( energy, -824.8591 , 4);


        # The total sum should be zero, no matter what
        Ftot = VECTOR(0.0, 0.0, 0.0)        
        for f in F:
            Ftot += f
            #print f.x, f.y, f.z
        
        self.assertAlmostEqual( Ftot.x, 0.0);
        self.assertAlmostEqual( Ftot.y, 0.0);
        self.assertAlmostEqual( Ftot.z, 0.0);

        """
        print stress.xx, stress.xy, stress.xz
        print stress.yx, stress.yy, stress.yz
        print stress.zx, stress.zy, stress.zz
        """

        # Stress is close to zero
        """
        self.assertAlmostEqual( stress.xy, 0.0);
        self.assertAlmostEqual( stress.xz, 0.0);
        self.assertAlmostEqual( stress.yx, 0.0);
        self.assertAlmostEqual( stress.yz, 0.0);
        self.assertAlmostEqual( stress.zx, 0.0);
        self.assertAlmostEqual( stress.zy, 0.0);

        self.assertAlmostEqual( stress.xx, stress.yy);
        self.assertAlmostEqual( stress.yy, stress.zz);       
        """
      
               


if __name__=='__main__':
    unittest.main()




