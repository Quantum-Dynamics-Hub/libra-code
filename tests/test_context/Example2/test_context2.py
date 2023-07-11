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

import os
import math
import sys
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



# Default (empty) context object
dctx = Context()

# Load an XML file to construct the context object
# and setup the path separator to be "/", because the tag names may contain
# "." in them - so the default path separator "." won't work

ctx = Context("x0.xml")      
ctx.set_path_separator("/")

print "path = ", ctx.get_path()   # what is the current path - the topmost tag name
print "size = ", ctx.size()       # how many children this node has


print "==========="
ctx1 = ctx.get_child("step", dctx)  # accessing the node's children - create a new Contxt object
print "path = ", ctx1.get_path()    # show the new object's path
print "size = ", ctx1.size()        # show the new object's number of children
ctx1.save_xml("_0.xml")             # save only that info into an XML file

ctx2 = ctx1.get_child("atomic_structure", dctx)    # we can repeat all the same as above, but at another level
print "path = ", ctx2.get_path()
print "size = ", ctx2.size()
ctx2.save_xml("_01.xml")



print "==========="
steps = ctx.get_children("step")                 # get all the children of the current node such that they all
                                                 # have the "step" tag. Create an array of Context objects
print "The number of <step> nodes =", len(steps)
print "path = ", steps[0].get_path()             
print "size = ", steps[0].size()
steps[0].save_xml("_1.xml")


ctx2 = steps[0].get_child("atomic_structure", dctx)   # access the elements of one of the objects in the created list
print "path = ", ctx2.get_path()
print "size = ", ctx2.size()
ctx2.save_xml("_11.xml")


print "==========="
print "path = ", ctx.get_path()   # what is the current path - the topmost tag name
all1 = ctx.get_children_all("qes:espresso")                     # get all the children of the current node 
print "The number of all nodes =", len(all1)

all1 = ctx.get_children_all()                     # get all the children of the current node 
print "The number of all nodes =", len(all1)


print "==========="
# Further processing:
# Note: so far, we can't access the nodes by giving complex path variables, like
# "atomic_structure/atomic_positions"
atoms = steps[0].get_children_all() 
atoms[0].show_children()

atoms = steps[0].get_child("atomic_structure",dctx).get_child("atomic_positions", dctx).get_children("atom")
print "The number of <atom> nodes =", len(atoms)

atoms[0].show_children()


# Now, access to the actual data (although in a form of a string)
nat = len(atoms)
for i in xrange(nat):
    print i, atoms[i].get("", ""), atoms[i].get("<xmlattr>/index", -1.0)
   

