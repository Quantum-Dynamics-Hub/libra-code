#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import math

# Fisrt, we add the location of the library to test to the PYTHON path
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


print "\nTest 2: Constructor"
print "rb = RigidBody()"
rb = RigidBody()

print "\nTest 3: Showing the content"
print "rb.show_info()"
rb.show_info()

print "\nTest 4: RB initialization"
print "masses = [1.0,1.0,1.0]"
print "centers = [VECTOR(0.0,0.0,0.0),VECTOR(1.0,0.0,0.0),VECTOR(0.0,1.0,0.0)]"
print "rb.init(3,masses,centers)"
masses = [1.0,1.0,1.0]
centers = [VECTOR(0.0,0.0,0.0),VECTOR(1.0,0.0,0.0),VECTOR(0.0,1.0,0.0)]
rb.init(3,masses,centers)
print "Show info again:"
print "rb.show_info()"
rb.show_info()


print "\nTest 5: Print properties"
print "rb.rb_centers = ",rb.rb_centers
print "rb.rb_mass = ",rb.rb_mass
print "rb.rb_iM = ",rb.rb_iM
print "rb.rb_cm = ",rb.rb_cm, rb.rb_cm.x, rb.rb_cm.y, rb.rb_cm.z
print "rb.rb_p = ",rb.rb_p, rb.rb_p.x, rb.rb_p.y, rb.rb_p.z
print "rb.rb_v = ",rb.rb_v
print "rb.rb_force = ",rb.rb_force
print "rb.rb_I_I = ",rb.rb_I_I
print "rb.rb_I_e = ",rb.rb_I_e
print "rb.rb_invI_I = ",rb.rb_invI_I
print "rb.rb_invI_e = ",rb.rb_invI_e
print "rb.rb_A = ",rb.rb_A
print "rb.rb_B = ",rb.rb_B
print "rb.rb_C = ",rb.rb_C
print "rb.rb_A_I_to_e = ",rb.rb_A_I_to_e
print "rb.rb_A_I_to_e_T = ",rb.rb_A_I_to_e_T
print "rb.rb_L = ",rb.rb_L
print "rb.rb_p_r = ",rb.rb_p_r
print "rb.rb_l_e = ",rb.rb_l_e
print "rb.rb_w_e = ",rb.rb_w_e
print "rb.rb_torque_e = ",rb.rb_torque_e
print "rb.is_fixed_translation = ",rb.is_fixed_translation
print "rb.is_fixed_rotation = ",rb.is_fixed_rotation


print "==== Need to continue this tutorial ===="


