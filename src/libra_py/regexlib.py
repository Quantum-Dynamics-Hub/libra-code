## \file regexlib.py 
# This module implements a collection of the general-purpose regular expressions

import re
import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


#------- Here are some basic patterns -------------
INT    = '([1-9]([0-9]*))'
NINT   = '([0-9]+)'
SP     = '\s+'    
DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
WORD   = '([a-zA-Z]+)'
ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
PHRASE = '"((\w|\W)+)"'
compINT = re.compile(INT)


#=== p - means "Pattern" =======================

pElement_name = '(?P<Element_name>'+WORD+')'+SP
pX_val = '(?P<X_val>'+DOUBLE+')'+SP
pY_val = '(?P<Y_val>'+DOUBLE+')'+SP
pZ_val = '(?P<Z_val>'+DOUBLE+')'+SP
