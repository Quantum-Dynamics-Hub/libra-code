###################################################################
# Tutorial: Setting interactions
###################################################################

import lcccsObjects
from DataBase.AtomicData.PeriodicTable import *
from DataBase.ForceFields.UFF import * #import LoadUFF
from DataBase.ForceFields.DREIDING import * #import LoadDREIDING
from DataBase.ForceFields.GAFF import * #import LoadGAFF
from DataBase.ForceFields.MMFF94 import * #import LoadMMFF94

##############################################################
import os
import LoadMolecule

#======= Set up the Universe here ================

U = lcccsObjects.Universe()
LoadPT.Load_PT(U)


i = 1
while i<=1:
#for i in [1]:
    print "=================== System ",i,"======================="
    # System creation section
    syst = lcccsObjects.System()

    # Create atoms, link them, define groups
    inp_file = os.getcwd()+"/Input/GAFF_Tests/test"+str(i)+"a.pdb"
    #inp_file = os.getcwd()+"/Input/Rings/test"+str(i)+".pdb"
    LoadMolecule.Load_Molecule(U, syst,inp_file,"pdb")

    print "======================== GROUP ATOMS ======================"
    syst.GROUP_ATOMS([1,2,3,4],1)
    print "======================== GROUP ATOMS ======================"
    syst.GROUP_ATOMS([5,6,7],2)
    print "======================== GROUP ATOMS ======================"
    syst.GROUP_ATOMS([8,9,10,11],3)

    # Create force field objects
    ff = lcccsObjects.ForceField()

    # Load parameters
    #LoadUFF.Load_UFF(ff)
    LoadDREIDING.Load_DREIDING(ff)
    #LoadGAFF.Load_GAFF(ff)

    # Set up functional forms
#    ff.set_functionals({"bond":"Harmonic","angle":"Fourier","vdw":"LJ","dihedral":"General0"}) # This is UFF functional
    ff.set_functionals({"bond":"Harmonic","angle":"Harmonic_Cos","vdw":"LJ","dihedral":"General2"}) # This is DREIDING functional


    # Define regions for different levels of theory
    atlst1 = range(1,syst.Number_of_atoms+1)

    # Create interactions between atoms from atom list 1
    syst.set_interactions_for_atoms(atlst1,atlst1,ff)    # The simplest sintax

#    syst.show_interactions()
#    syst.show_pairs()
#    syst.show_frag_pairs()
#    syst.show_rings()
#    syst.show_atoms()

    syst.zero_atom_forces()
    print "system bond energy = ", syst.energy("bond")/627.5094709
    print "system angle energy = ", syst.energy("angle")/627.5094709
    print "system dihedral energy = ", syst.energy("dihedral")/627.5094709
    print "system vdw energy = ", syst.energy("vdw")/627.5094709

    print "system energy = ",syst.energy()/627.5094709

    syst.init_fragments()
    syst.show_fragments()
    

    print "========================================================="
    i = i + 1


