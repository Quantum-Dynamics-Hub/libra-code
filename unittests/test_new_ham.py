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
sys.path.insert(1,cwd+"/../_build/src/hamiltonian/nHamiltonian_Generic")
sys.path.insert(1,cwd+"/../_build/src/converters")
sys.path.insert(1,cwd+"/../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cygnhamiltonian_generic import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libnhamiltonian_generic import *
    from liblinalg import *




def model1(Hdia, Sdia, d1ham_dia, dc1_dia, x,x0,k,D,V):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    2.0*k*(0.5^2 + 2.0*0.5*dx + dx^2) +
    2.0*k*((-0.5)^2 - 2.0*0.5*dx + dx^2) 
    so forces at 0.5 should be zero!

    """

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );  Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    



class Test_New_Hamiltonian(unittest.TestCase):
    """ Summary of the tests:
    """

    def test_1(self):
        """Set by reference"""

        ham = nHamiltonian(2,2,1)  

        Hdia = CMATRIX(2,2)
        Hadi = CMATRIX(2,2)
        Sdia = CMATRIX(2,2)
        U = CMATRIX(2,2)
        Ddia = CMATRIX(2,2)
        Dadi = CMATRIX(2,2)



        ham.set_ham_dia_by_ref(Hdia);  
        ham.set_ham_adi_by_ref(Hadi);  
        ham.set_ovlp_dia_by_ref(Sdia); 
        ham.set_basis_transform_by_ref(U); 
        ham.set_den_mat_dia_by_ref(Ddia);  
        ham.set_den_mat_adi_by_ref(Dadi);  

        
        # This part is trivial
        for i in [0,1]:
            for j in [0,1]:
                self.assertAlmostEqual( ham.get_ham_dia().get(i,j) , Hdia.get(i,j) )
                self.assertAlmostEqual( ham.get_ham_adi().get(i,j) , Hadi.get(i,j) )
                self.assertAlmostEqual( ham.get_ovlp_dia().get(i,j) , Sdia.get(i,j) )
                self.assertAlmostEqual( ham.get_basis_transform().get(i,j) , U.get(i,j) )
                self.assertAlmostEqual( ham.get_den_mat_dia().get(i,j) , Ddia.get(i,j) )
                self.assertAlmostEqual( ham.get_den_mat_adi().get(i,j) , Dadi.get(i,j) )


        # Now, lest channge the outside matrices
        for i in [0,1]:
            for j in [0,1]:

                Hdia.set(i,j, 1.0+i,2.0+j)
                Hadi.set(i,j, i+j, i-2.0-j)
                Sdia.set(i,i, 1.0, 0.0)
                U.set(i,j, i*j, 2.0*i-3*j)
                Dadi.set(i,j, 10*(i+1.0), 20.0*(j-2))
                Ddia.set(i,j, 10*(i+1.0), 20.0*(j-2))


        # This is not so trivial anymore
        for i in [0,1]:
            for j in [0,1]:
                self.assertAlmostEqual( ham.get_ham_dia().get(i,j) , Hdia.get(i,j) )
                self.assertAlmostEqual( ham.get_ham_adi().get(i,j) , Hadi.get(i,j) )
                self.assertAlmostEqual( ham.get_ovlp_dia().get(i,j) , Sdia.get(i,j) )
                self.assertAlmostEqual( ham.get_basis_transform().get(i,j) , U.get(i,j) )
                self.assertAlmostEqual( ham.get_den_mat_dia().get(i,j) , Ddia.get(i,j) )
                self.assertAlmostEqual( ham.get_den_mat_adi().get(i,j) , Dadi.get(i,j) )




    def test_2(self):
        """Set by reference:  the 1-st order derivatives"""

        ham = nHamiltonian(2,2,1)  

        d1ham_dia = CMATRIXList()
        d1ham_adi = CMATRIXList()
        dc1_dia = CMATRIXList()
        dc1_adi = CMATRIXList()
        for i in [0]:
            d1ham_dia.append( CMATRIX(2,2) )
            d1ham_adi.append( CMATRIX(2,2) )
            dc1_dia.append( CMATRIX(2,2) )
            dc1_adi.append( CMATRIX(2,2) )

        ham.set_d1ham_dia_by_ref(d1ham_dia)
        ham.set_d1ham_adi_by_ref(d1ham_adi)
        ham.set_dc1_dia_by_ref(dc1_dia)
        ham.set_dc1_adi_by_ref(dc1_adi)

        
        # This part is trivial
        for i in [0,1]:
            for j in [0,1]:
                for k in [0]:
                    self.assertAlmostEqual( ham.get_d1ham_dia(k).get(i,j) , d1ham_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_d1ham_adi(k).get(i,j) , d1ham_adi[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_dc1_dia(k).get(i,j) , dc1_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_dc1_adi(k).get(i,j) , dc1_adi[k].get(i,j) )


        # Now, lest channge the outside matrices
        for i in [0,1]:
            for j in [0,1]:
                for k in [0]:
                    d1ham_dia[k].set(i,j, i+k*j, (j+1)*k-i)
                    d1ham_adi[k].set(i,j, i+k*j, (j+1)*k-i)
                    dc1_dia[k].set(i,j, i+k*j, (j+1)*k-i)
                    dc1_adi[k].set(i,j, i+k*j, (j+1)*k-i)


        # This is not so trivial anymore
        for i in [0,1]:
            for j in [0,1]:
                for k in [0]:
                    self.assertAlmostEqual( ham.get_d1ham_dia(k).get(i,j) , d1ham_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_d1ham_adi(k).get(i,j) , d1ham_adi[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_dc1_dia(k).get(i,j) , dc1_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_dc1_adi(k).get(i,j) , dc1_adi[k].get(i,j) )



    def test_3(self):
        """Set by reference:  2-nd order derivatives"""

        ham = nHamiltonian(2,2,2)  

        d2ham_dia = CMATRIXList()
        d2ham_adi = CMATRIXList()

        for i in xrange(2*2):
            d2ham_dia.append( CMATRIX(2,2) )
            d2ham_adi.append( CMATRIX(2,2) )

        ham.set_d2ham_dia_by_ref(d2ham_dia)
        ham.set_d2ham_adi_by_ref(d2ham_adi)

        
        # This part is trivial
        for i in [0,1]:
            for j in [0,1]:
                for k in xrange(2*2):
                    self.assertAlmostEqual( ham.get_d2ham_dia(k).get(i,j) , d2ham_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_d2ham_adi(k).get(i,j) , d2ham_adi[k].get(i,j) )

        # Now, lest channge the outside matrices
        for i in [0,1]:
            for j in [0,1]:
                for k in xrange(2*2):
                    d2ham_dia[k].set(i,j, i+k*j, (j+1)*k-i)
                    d2ham_adi[k].set(i,j, i+k*j, (j+1)*k-i)


        # This is not so trivial anymore
        for i in [0,1]:
            for j in [0,1]:
                for k in xrange(2*2):
                    self.assertAlmostEqual( ham.get_d2ham_dia(k).get(i,j) , d2ham_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_d2ham_adi(k).get(i,j) , d2ham_adi[k].get(i,j) )




    def test_4(self):
        """Set by value: matrices"""

        ham = nHamiltonian(2,2,1)  

        Hdia = CMATRIX(2,2)
        Hadi = CMATRIX(2,2)
        Sdia = CMATRIX(2,2)
        U = CMATRIX(2,2)
        Ddia = CMATRIX(2,2)
        Dadi = CMATRIX(2,2)

        ham.set_ham_dia_by_val(Hdia);  
        ham.set_ham_adi_by_val(Hadi);  
        ham.set_ovlp_dia_by_val(Sdia); 
        ham.set_basis_transform_by_val(U); 
        ham.set_den_mat_dia_by_val(Ddia);  
        ham.set_den_mat_adi_by_val(Dadi);  

        
        # This part is trivial
        for i in [0,1]:
            for j in [0,1]:
                self.assertAlmostEqual( ham.get_ham_dia().get(i,j) , Hdia.get(i,j) )
                self.assertAlmostEqual( ham.get_ham_adi().get(i,j) , Hadi.get(i,j) )
                self.assertAlmostEqual( ham.get_ovlp_dia().get(i,j) , Sdia.get(i,j) )
                self.assertAlmostEqual( ham.get_basis_transform().get(i,j) , U.get(i,j) )
                self.assertAlmostEqual( ham.get_den_mat_dia().get(i,j) , Ddia.get(i,j) )
                self.assertAlmostEqual( ham.get_den_mat_adi().get(i,j) , Dadi.get(i,j) )


        # Now, lest channge the outside matrices
        for i in [0,1]:
            for j in [0,1]:

                Hdia.set(i,j, 1.0+i,2.0+j)
                Hadi.set(i,j, i+j, i-2.0-j)
                Sdia.set(i,i, 1.0, 0.0)
                U.set(i,j, i*j, 2.0*i-3*j)
                Dadi.set(i,j, 10*(i+1.0), 20.0*(j-2))
                Ddia.set(i,j, 10*(i+1.0), 20.0*(j-2))


        # This is not so trivial anymore
        for i in [0,1]:
            for j in [0,1]:
                self.assertAlmostEqual( ham.get_ham_dia().get(i,j) , 0.0+0.0j )
                self.assertAlmostEqual( ham.get_ham_adi().get(i,j) , 0.0+0.0j )
                self.assertAlmostEqual( ham.get_ovlp_dia().get(i,j) , 0.0+0.0j )
                self.assertAlmostEqual( ham.get_basis_transform().get(i,j) , 0.0+0.0j )
                self.assertAlmostEqual( ham.get_den_mat_dia().get(i,j) , 0.0+0.0j )
                self.assertAlmostEqual( ham.get_den_mat_adi().get(i,j) , 0.0+0.0j )





    def test_5(self):
        """Set by value:  now also the derivatives"""

        ham = nHamiltonian(2,2,1)  

        d1ham_dia = CMATRIXList()
        d1ham_adi = CMATRIXList()
        dc1_dia = CMATRIXList()
        dc1_adi = CMATRIXList()
        for i in [0]:
            d1ham_dia.append( CMATRIX(2,2) )
            d1ham_adi.append( CMATRIX(2,2) )
            dc1_dia.append( CMATRIX(2,2) )
            dc1_adi.append( CMATRIX(2,2) )

        ham.set_d1ham_dia_by_val(d1ham_dia)
        ham.set_d1ham_adi_by_val(d1ham_adi)
        ham.set_dc1_dia_by_val(dc1_dia)
        ham.set_dc1_adi_by_val(dc1_adi)

        
        # This part is trivial
        for i in [0,1]:
            for j in [0,1]:
                for k in [0]:
                    self.assertAlmostEqual( ham.get_d1ham_dia(k).get(i,j) , d1ham_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_d1ham_adi(k).get(i,j) , d1ham_adi[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_dc1_dia(k).get(i,j) , dc1_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_dc1_adi(k).get(i,j) , dc1_adi[k].get(i,j) )


        # Now, lest channge the outside matrices
        for i in [0,1]:
            for j in [0,1]:
                for k in [0]:
                    d1ham_dia[k].set(i,j, i+k*j, (j+1)*k-i)
                    d1ham_adi[k].set(i,j, i+k*j, (j+1)*k-i)
                    dc1_dia[k].set(i,j, i+k*j, (j+1)*k-i)
                    dc1_adi[k].set(i,j, i+k*j, (j+1)*k-i)


        # This is not so trivial anymore
        for i in [0,1]:
            for j in [0,1]:
                for k in [0]:
                    self.assertAlmostEqual( ham.get_d1ham_dia(k).get(i,j) , 0.0+0.0j )
                    self.assertAlmostEqual( ham.get_d1ham_adi(k).get(i,j) , 0.0+0.0j )
                    self.assertAlmostEqual( ham.get_dc1_dia(k).get(i,j) , 0.0+0.0j )
                    self.assertAlmostEqual( ham.get_dc1_adi(k).get(i,j) , 0.0+0.0j )


    def test_6(self):
        """Set by value:  2-nd order derivatives"""

        ham = nHamiltonian(2,2,2)  

        d2ham_dia = CMATRIXList()
        d2ham_adi = CMATRIXList()

        for i in xrange(2*2):
            d2ham_dia.append( CMATRIX(2,2) )
            d2ham_adi.append( CMATRIX(2,2) )

        ham.set_d2ham_dia_by_val(d2ham_dia)
        ham.set_d2ham_adi_by_val(d2ham_adi)

        
        # This part is trivial
        for i in [0,1]:
            for j in [0,1]:
                for k in xrange(2*2):
                    self.assertAlmostEqual( ham.get_d2ham_dia(k).get(i,j) , d2ham_dia[k].get(i,j) )
                    self.assertAlmostEqual( ham.get_d2ham_adi(k).get(i,j) , d2ham_adi[k].get(i,j) )

        # Now, lest channge the outside matrices
        for i in [0,1]:
            for j in [0,1]:
                for k in xrange(2*2):
                    d2ham_dia[k].set(i,j, i+k*j, (j+1)*k-i)
                    d2ham_adi[k].set(i,j, i+k*j, (j+1)*k-i)


        # This is not so trivial anymore
        for i in [0,1]:
            for j in [0,1]:
                for k in xrange(2*2):
                    self.assertAlmostEqual( ham.get_d2ham_dia(k).get(i,j) , 0.0+0.0j )
                    self.assertAlmostEqual( ham.get_d2ham_adi(k).get(i,j) , 0.0+0.0j )



    def test_7(self):
        """ Diabatic-to-adiabatic : only matrices"""

        ham = nHamiltonian(2,2,1)  

        Hdia = CMATRIX(2,2)
        Hdia.set(0,0, 1.0+0.0j);  Hdia.set(0,1, 0.1+0.0j);
        Hdia.set(1,0, 0.1+0.0j);  Hdia.set(1,1, 1.0+0.0j);

        Sdia = CMATRIX(2,2)
        Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
        Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

        Hadi = CMATRIX(2,2)
        U = CMATRIX(2,2)

        ham.set_ham_dia_by_ref(Hdia);  
        ham.set_ham_adi_by_ref(Hadi);  
        ham.set_ovlp_dia_by_ref(Sdia); 
        ham.set_basis_transform_by_ref(U); 


        ham.compute_adiabatic(0);         
        self.assertAlmostEqual( Hadi.get(0,0) , 0.9+0.0j )
        self.assertAlmostEqual( Hadi.get(1,1) , 1.1+0.0j )


        Hdia.set(0,1, 0.5+0.0j)
        Hdia.set(1,0, 0.5+0.0j)

        ham.compute_adiabatic(0);         
        self.assertAlmostEqual( Hadi.get(0,0) , 0.5+0.0j )
        self.assertAlmostEqual( Hadi.get(1,1) , 1.5+0.0j )

        x = U.H() * Sdia * U
        self.assertAlmostEqual( x.get(0,0) , 1.0+0.0j )
        self.assertAlmostEqual( x.get(0,1) , 0.0+0.0j )
        self.assertAlmostEqual( x.get(1,0) , 0.0+0.0j )
        self.assertAlmostEqual( x.get(1,1) , 1.0+0.0j )



    def test_8(self):
        """ Diabatic-to-adiabatic : also the 1-st order derivatives 

        (-1 - x) * (1-x) - V^2 = 0
        (1+x) * (1-x) = -V^2
        1 + V^2 = x^2        
        """

        ham = nHamiltonian(2,2,1)  

        Hdia = CMATRIX(2,2)
        Hdia.set(0,0,-1.0+0.0j);  Hdia.set(0,1, 0.5+0.0j);
        Hdia.set(1,0, 0.5+0.0j);  Hdia.set(1,1, 1.0+0.0j);

        Sdia = CMATRIX(2,2)
        Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
        Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

        Hadi = CMATRIX(2,2)
        U = CMATRIX(2,2)

        ham.set_ham_dia_by_ref(Hdia);  
        ham.set_ham_adi_by_ref(Hadi);  
        ham.set_ovlp_dia_by_ref(Sdia); 
        ham.set_basis_transform_by_ref(U); 

        d1ham_dia = CMATRIXList()
        d1ham_adi = CMATRIXList()
        dc1_dia = CMATRIXList()
        dc1_adi = CMATRIXList()
        for i in [0]:
            d1ham_dia.append( CMATRIX(2,2) )
            d1ham_adi.append( CMATRIX(2,2) )
            dc1_dia.append( CMATRIX(2,2) )
            dc1_adi.append( CMATRIX(2,2) )

            #  d Hdia / dR_0
            d1ham_dia[i].set(0,0, 1.0+0.0j);   d1ham_dia[i].set(0,1, 0.1+0.0j);
            d1ham_dia[i].set(1,0, 0.1+0.0j);   d1ham_dia[i].set(1,1,-1.1+0.0j);

            #  <dia| d/dR_0| dia >
            dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.1+0.0j);
            dc1_dia[i].set(1,0,-0.1+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

       
        ham.set_d1ham_dia_by_ref(d1ham_dia)
        ham.set_d1ham_adi_by_ref(d1ham_adi)
        ham.set_dc1_dia_by_ref(dc1_dia)
        ham.set_dc1_adi_by_ref(dc1_adi)


        ham.compute_adiabatic(1);         
        self.assertAlmostEqual( Hadi.get(0,0) ,-math.sqrt(1.25)*(1.0+0.0j) )
        self.assertAlmostEqual( Hadi.get(1,1) , math.sqrt(1.25)*(1.0+0.0j) )

        for i in [0]:
            print "Adi dE/dR =  "; d1ham_adi[i].show_matrix()
            print "Adi D = "; dc1_adi[i].show_matrix()



        

    def test_9(self):
        """ Ehrenfest energies : only matrices

              |  1   0.1  |   (1)          (1.1)
        (1 1) |           |   (1)  = (1 1) (1.1) = 2.2    Norm: 2, so E = 1.1
              |  0.1   1  |

        |adi> = |dia> U

        |adi> * C_adi = |dia> * C_dia

        so: U * C_adi = C_dia =>  C_adi = U^-1 * C_dia = U^+ * S * C_dia

        """

        ham = nHamiltonian(2,2,1)  

        Hdia = CMATRIX(2,2)
        Hdia.set(0,0, 1.0+0.0j);  Hdia.set(0,1, 0.1+0.0j);
        Hdia.set(1,0, 0.1+0.0j);  Hdia.set(1,1, 1.0+0.0j);

        Sdia = CMATRIX(2,2)
        Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
        Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

        Hadi = CMATRIX(2,2)
        U = CMATRIX(2,2)

        ham.set_ham_dia_by_ref(Hdia);  
        ham.set_ham_adi_by_ref(Hadi);  
        ham.set_ovlp_dia_by_ref(Sdia); 
        ham.set_basis_transform_by_ref(U); 

        ham.compute_adiabatic(0);   


        Cdia = CMATRIX(2,1)
        Cdia.set(0, 1.0+0.0j); Cdia.set(1, 1.0+0.0j)
        ham.set_ampl_dia_by_ref(Cdia)
        
        self.assertAlmostEqual( ham.Ehrenfest_energy_dia(), 1.1+0.0j )

      
        Cadi = U.H() * Sdia * Cdia
        ham.set_ampl_adi_by_ref(Cadi)


        self.assertAlmostEqual( ham.Ehrenfest_energy_adi(), 1.1+0.0j )

        



    def test_10(self):
        """ Ehrenfest energies : also 1-st order derivatives

              |  1   0.1  |   (1)          (1.1)
        (1 1) |           |   (1)  = (1 1) (1.1) = 2.2    Norm: 2, so E = 1.1
              |  0.1   1  |

        |adi> = |dia> U

        |adi> * C_adi = |dia> * C_dia

        so: U * C_adi = C_dia =>  C_adi = U^-1 * C_dia = U^+ * S * C_dia

        """

        ham = nHamiltonian(2,2,1)  


        # Allocate memory
        Hdia = CMATRIX(2,2);   Sdia = CMATRIX(2,2);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
        Hadi = CMATRIX(2,2);   ham.set_ham_adi_by_ref(Hadi);  
        U = CMATRIX(2,2);      ham.set_basis_transform_by_ref(U); 
        Cdia = CMATRIX(2,1);   ham.set_ampl_dia_by_ref(Cdia)
        Cadi = CMATRIX(2,1);   ham.set_ampl_adi_by_ref(Cadi)


        d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
        dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
        for i in [0]:
            d1ham_dia.append( CMATRIX(2,2) ); d1ham_adi.append( CMATRIX(2,2) )
            dc1_dia.append( CMATRIX(2,2) );   dc1_adi.append( CMATRIX(2,2) )
       
        ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
        ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)


        #  Set up the models and compute internal variables
        Cadi.set(0, 1.0+0.0j); Cadi.set(1, 1.0+0.0j)        
        tmp = U * Cadi;      Cdia.set(0, tmp.get(0)); Cdia.set(1, tmp.get(1));
#        tmp = U.H() * Sdia * Cdia;   Cadi.set(0, tmp.get(0)); Cadi.set(1, tmp.get(1));

        # def model1(Hdia, Sdia, d1ham_dia, dc1_dia, x,x0,k,D,V):
        model1(Hdia, Sdia, d1ham_dia, dc1_dia, 0.5, 1.0, 1.0, -0.1, 0.1)
        ham.compute_adiabatic(1);   

        tmp = U * Cadi;      Cdia.set(0, tmp.get(0)); Cdia.set(1, tmp.get(1));

#        tmp = U.H() * Sdia * Cdia;   Cadi.set(0, tmp.get(0)); Cadi.set(1, tmp.get(1));

        Edia1 = ham.Ehrenfest_energy_dia()
        Eadi1 = ham.Ehrenfest_energy_adi()

        print "Hdia = "; Hdia.show_matrix()
        print "Hadi = "; Hadi.show_matrix()
        print "Cdia = "; Cdia.show_matrix()
        print "Cadi = "; Cadi.show_matrix()

        print "d1ham_dia = "; d1ham_dia[0].show_matrix()
        print "d1ham_adi = "; d1ham_adi[0].show_matrix()
        print "dc1_dia = "; dc1_dia[0].show_matrix()
        print "dc1_adi = "; dc1_adi[0].show_matrix()

        
        print "Edia1 = ", Edia1
        print "Eadi1 = ", Eadi1

        
        f_dia1 = ham.forces_tens_dia()
        f_adi1 = ham.forces_tens_adi()

        nrm_dia = (Cdia.H() * Sdia * Cdia).get(0)
        nrm_adi = (Cadi.H() * Cadi).get(0)

        print "Forces dia = "; (Cdia.H() * f_dia1[0] * Cdia / nrm_dia ).show_matrix()
        print "Forces adi = "; (Cadi.H() * f_adi1[0] * Cadi / nrm_adi ).show_matrix()



        model1(Hdia, Sdia, d1ham_dia, dc1_dia, 0.5+0.001, 1.0, 1.0, -0.1, 0.1)
        ham.compute_adiabatic(1);   

        tmp = U * Cadi;      Cdia.set(0, tmp.get(0)); Cdia.set(1, tmp.get(1));

#        tmp = U * Cdia
#        Cadi.set(0, tmp.get(0)); Cadi.set(1, tmp.get(1));
#        tmp = U.H() * Sdia * Cdia;   Cadi.set(0, tmp.get(0)); Cadi.set(1, tmp.get(1));

        Edia2 = ham.Ehrenfest_energy_dia()
        Eadi2 = ham.Ehrenfest_energy_adi()

        print "Hdia = "; Hdia.show_matrix()
        print "Hadi = "; Hadi.show_matrix()
        print "Cdia = "; Cdia.show_matrix()
        print "Cadi = "; Cadi.show_matrix()

        print "Edia2 = ", Edia2
        print "Eadi2 = ", Eadi2

        print "Numerical force (dia) = ", (Edia2 - Edia1)/0.001
        print "Numerical force (adi) = ", (Eadi2 - Eadi1)/0.001
        
      

      
         



                      


if __name__=='__main__':
    unittest.main()



