#*********************************************************************************                     
#* Copyright (C) 2017 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file parse_gamess.py
# This module implements reading a GAMESS format basis as obtained from EMSL
# or maybe from the GAMESS output. This file is a significantly re-written 
# version of the Basis/Tools.py file from the PyQuante library. See the
# license notice below.


"""\
This program inputs a GAMESS basis, e.g. from the EMSL basis set order
form, and formats it into PyQuante format

Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

import re,sys

def name2no(name):
    """Converts the name of an element to its atomic number. The input name can be upper, lower, or mixed case
    """
        
    all_names = [
    "dummy",
    "hydrogen", "helium",
    "lithium","beryllium","boron","carbon","nitrogen",
    "oxygen","fluorine","neon","sodium","magnesium",
    "aluminum","silicon","phosphorus","sulfur","chlorine",
    "argon","potassium", "calcium", "scandium", "titanium",
    "vanadium", "chromium", "manganese", "iron",
    "cobalt", "nickel", "copper", "zinc",
    "gallium", "germanium", "arsenic", "selenium", "bromine",
    "krypton", "rubidium", "strontium", "yttrium", "zirconium",
    "niobium", "molybdenum", "technetium", "ruthenium","rhodium",
    "palladium", "silver", "cadmium","indium", "tin", "antimony",
    "tellerium", "iodine", "xenon","cesium", "barium",
    "lanthanum","cerium","praseodymium","neodymium","promethium",
    "samarium","europium","gadolinium","terbium","dysprosium",
    "holmium","erbium","thulium","ytterbium","lutetium",
    "halfnium","tantalum","tungsten","rhenium","osmium","iridium",
    "platinum","gold","mercury","thallium","lead","bismuth",
    "polonium","astatine","radon"]

    res = -1
    lname = name.lower()

    if lname in all_names:
        res = all_names.index(lname) 

    else:
        print("Error in name2no function: The element ", name, " is not defined")
        print("Possible values are: ", all_names)
        sys.exit(0)

    return res



def parse_gamess_basis(filename):
    """Read a file <filename> in GAMESS format containing the info about basis
       for a number of elements. The function returns a dictionary containing 
       basis info for all read elements
    """

    # Read all lines
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    # Find the section which contains the info about basis functions
    # This is determined by the $DATA and $END blocks
    i1, i2, i = 0, 0, 0
    for line in lines:
        if re.search('\$DATA',line):
            i1 = i
        if re.search('\$END',line):
            i2 = i
        i = i + 1

    # Prepare to populate the basis for all elements
    basis = {} 
    atno = 0  # current element we are working on

    # Look in the correct info block!
    i = i1+1
    while i<i2:
         
        words = lines[i].split()
        sz = len(words)

        if sz==0:
            pass

        elif sz==1:
            el_string = words[0]
            atno = name2no(el_string)
            basis.setdefault(atno, [] )

        else:
            symb = words[0]
            nprim = int(words[1])
            if symb == "L":
                sprims = []
                pprims = []
                for j in range(0,nprim):
                    i = i + 1
                    line = lines[i]
                    words = line.split()
                    sprims.append((float(words[1]),float(words[2])))
                    pprims.append((float(words[1]),float(words[3])))
                basis[atno].append(("S",sprims))
                basis[atno].append(("P",pprims))
            else:
                prims = []
                for j in range(0,nprim):
                    i = i + 1
                    line = lines[i]
                    words = line.split()
                    prims.append((float(words[1]),float(words[2])))
                basis[atno].append((symb,prims))

        i = i + 1 
    return basis


#if __name__ == '__main__': 
#    # Eample usage
#    print name2no("Hydrogen")
#    bas = parse_gamess_basis("sto_2g.txt")
#    print bas


