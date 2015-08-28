###################################################################
# Tutorial: Grouping atoms and functional group determination
###################################################################

import lcccsObjects
from DataBase.ForceFields.UFF import * #import LoadUFF
from DataBase.ForceFields.DREIDING import * #import LoadDREIDING
from DataBase.ForceFields.GAFF import * #import LoadGAFF
from DataBase.ForceFields.MMFF94 import * #import LoadMMFF94

##############################################################
import os
import LoadMolecule

#i = 1
#while i<=12:
for i in [1]:
    print "=================== System ",i,"======================="
    # System creation section
    syst = lcccsObjects.System()

    # Create atoms, link them, define groups
    #inp_file = os.getcwd()+"/Input/GAFF_Tests/test"+str(i)+"a.pdb"
    inp_file = os.getcwd()+"/Input/Rings/test"+str(i)+".pdb"
    LoadMolecule.Load_Molecule(syst,inp_file,"pdb")

    # Create force field objects
    ff = lcccsObjects.ForceField()

    # Load parameters
    LoadUFF.Load_UFF(ff)
    #LoadGAFF.Load_GAFF(ff)

    # Set up functional forms
    ff.set_functionals({"bond":"Harmonic","angle":"Harmonic","vdw":"LJ12_6"})

    # Define regions for different levels of theory
    atlst1 = range(1,syst.Number_of_atoms+1)

    # Create interactions between atoms from atom list 1
    syst.set_interactions_for_atoms(atlst1,atlst1,ff)    # The simplest sintax
    syst.show_rings()
    syst.show_atoms()
    

    print "========================================================="
    i = i + 1


