
import sys
mods = [m.__name__ for m in sys.modules.values() if m]

print mods


if sys.platform=="cygwin":
    from cyglibra_core import *
    #import cyglibra_core
    #print dir(cyglibra_core)

elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


#mods = [m.__name__ for m in cyglibra_core.modules.values() if m]

#print mods
