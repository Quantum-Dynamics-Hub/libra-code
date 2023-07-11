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







class Test_Ewald_VdW(unittest.TestCase):
    """ Summary of the tests:
    """


    def test_1(self):
        """Computing the vdW (dispersion) sum for NaCl system
           with 4 NaCl molecules in the unit cell (8 atoms)
           The refrence energy value is close to the one reported 
           in the paper:
           Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for 
           Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327 

           Here we take general form of the pairwise coefficients

        """

        Angst = 1.889725989       # 1 Angstrom in atomic units
        hartree = 627.5094709     # 1 Ha = 627.5.. kcal/mol
        
        R = VECTORList()
        F = VECTORList()
        types = intList()
        Bij = doubleList() 
        
        a = 5.63 * Angst   # in a.u.
        etha = 2.5 * Angst # in a.u.
        
        tv1 = a*VECTOR(1.0, 0.0, 0.0)
        tv2 = a*VECTOR(0.0, 1.0, 0.0)
        tv3 = a*VECTOR(0.0, 0.0, 1.0)
        
        box = MATRIX3x3(tv1, tv2, tv3)
        stress = MATRIX3x3()
        # 1 kcal/mol * A^6  = (1.0/hartree) * (1.0/Angstrom)**6
        B_NaNa =  24.180 * (1.0/hartree) * (Angst)**6
        B_ClCl =  1669.58 * (1.0/hartree) * (Angst)**6
        B_NaCl =  161.20 * (1.0/hartree) * (Angst)**6

        types.append(0)
        types.append(0)
        types.append(0)
        types.append(0)

        types.append(1)
        types.append(1)
        types.append(1)
        types.append(1)

        Bij.append(B_NaNa)
        Bij.append(B_NaCl)
        Bij.append(B_NaCl)
        Bij.append(B_ClCl)

        
        R.append( a*VECTOR(0.0, 0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));    # Na
        R.append( a*VECTOR(0.5, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));    # Na
        R.append( a*VECTOR(0.5, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));    # Na
        R.append( a*VECTOR(0.0, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));    # Na

        R.append( a*VECTOR(0.5, 0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));    # Cl
        R.append( a*VECTOR(0.0, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));    # Cl
        R.append( a*VECTOR(0.0, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));    # Cl
        R.append( a*VECTOR(0.5, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));    # Cl
                                                              
        
        pbc_deg = 2
        rec_deg = 2
        R_on = 20 * Angst
        R_off = 22 * Angst

        # We'll get energies in kcal/mol        
        energy = VdW_Ewald3D(R, types, 2, Bij, box, F, stress, rec_deg, pbc_deg, etha, R_on, R_off ) * hartree

        print energy
#        self.assertAlmostEqual( energy, -824.8591 , 4);


        # The total sum should be zero, no matter what
        Ftot = VECTOR(0.0, 0.0, 0.0)        
        for f in F:
            Ftot += f
            print f.x, f.y, f.z
        
#        self.assertAlmostEqual( Ftot.x, 0.0);
#        self.assertAlmostEqual( Ftot.y, 0.0);
#        self.assertAlmostEqual( Ftot.z, 0.0);

        # Stress is close to zero
#        self.assertAlmostEqual( stress.xy, 0.0);
#        self.assertAlmostEqual( stress.xz, 0.0);
#        self.assertAlmostEqual( stress.yx, 0.0);
#        self.assertAlmostEqual( stress.yz, 0.0);
#        self.assertAlmostEqual( stress.zx, 0.0);
#        self.assertAlmostEqual( stress.zy, 0.0);

#        self.assertAlmostEqual( stress.xx, stress.yy);
#        self.assertAlmostEqual( stress.yy, stress.zz);
       


    def test_2(self):
        """Computing the electrostatic sum for NaCl system
           with 4 NaCl molecules in the unit cell (8 atoms)
           The refrence energy value is close to the one reported 
           in the paper:
           Karasawa, N.; Goddard III,W. A. "Acceleration of Convergence for 
           Lattice Sums" J.Phys.Chem. 1989, 93,7320-7327 

           Here we assume geometric mean is applied to the coefficients
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
        # 1 kcal/mol * A^6  = (1.0/hartree) * (1.0/Angstrom)**6
        B_NaNa = math.sqrt( 24.180 * (1.0/hartree) * (Angst)**6)
        B_ClCl = math.sqrt( 1669.58 * (1.0/hartree) * (Angst)**6)

        print math.sqrt(24.180 * 1669.58)
        
        R.append( a*VECTOR(0.0, 0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_NaNa)   # Na
        R.append( a*VECTOR(0.5, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_NaNa)   # Na
        R.append( a*VECTOR(0.5, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_NaNa)   # Na
        R.append( a*VECTOR(0.0, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_NaNa)   # Na

        R.append( a*VECTOR(0.5, 0.0, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_ClCl)   # Cl
        R.append( a*VECTOR(0.0, 0.5, 0.0));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_ClCl)   # Cl
        R.append( a*VECTOR(0.0, 0.0, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_ClCl)   # Cl
        R.append( a*VECTOR(0.5, 0.5, 0.5));  F.append( VECTOR(0.0, 0.0, 0.0));  Q.append( B_ClCl)   # Cl
                                                              
        
        pbc_deg = 20
        rec_deg = 10
        R_on = 200 * Angst
        R_off = 220 * Angst

        # We'll get energies in kcal/mol        
        energy = VdW_Ewald3D(R, Q, box, F, stress, rec_deg, pbc_deg, etha, R_on, R_off ) * hartree

        print energy
#        self.assertAlmostEqual( energy, -824.8591 , 4);


        # The total sum should be zero, no matter what
        Ftot = VECTOR(0.0, 0.0, 0.0)        
        for f in F:
            Ftot += f
            print f.x, f.y, f.z
        
#        self.assertAlmostEqual( Ftot.x, 0.0);
#        self.assertAlmostEqual( Ftot.y, 0.0);
#        self.assertAlmostEqual( Ftot.z, 0.0);

        # Stress is close to zero
#        self.assertAlmostEqual( stress.xy, 0.0);
#        self.assertAlmostEqual( stress.xz, 0.0);
#        self.assertAlmostEqual( stress.yx, 0.0);
#        self.assertAlmostEqual( stress.yz, 0.0);
#        self.assertAlmostEqual( stress.zx, 0.0);
#        self.assertAlmostEqual( stress.zy, 0.0);

#        self.assertAlmostEqual( stress.xx, stress.yy);
#        self.assertAlmostEqual( stress.yy, stress.zz);
       


      
               


if __name__=='__main__':
    unittest.main()




