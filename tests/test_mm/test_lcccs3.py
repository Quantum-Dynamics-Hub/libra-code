###################################################################
# Tutorial: Grouping atoms and functional group determination
###################################################################

import lcccsObjects

##############################################################
# System creation section
syst = lcccsObjects.System()
# Create atoms, link them, define groups

syst.CREATE_ATOM( lcccsObjects.Atom({"Atom_element":"C","Atom_cm_x":0.0,"Atom_cm_y":0.0,"Atom_cm_z":0.0})  )

syst.CREATE_ATOM( lcccsObjects.Atom({"Atom_element":"H","Atom_cm_x":-1.0,"Atom_cm_y":-1.0,"Atom_cm_z":0.0}) )

at = lcccsObjects.Atom({"Atom_element":"H","Atom_cm_x":-1.0,"Atom_cm_y":1.0,"Atom_cm_z":0.0}) 
syst.CREATE_ATOM(at)

at = lcccsObjects.Atom({"Atom_element":"C","Atom_cm_x":2.0,"Atom_cm_y":0.0,"Atom_cm_z":0.0})    
syst.CREATE_ATOM(at)

at = lcccsObjects.Atom({"Atom_element":"H","Atom_cm_x":3.0,"Atom_cm_y":1.0,"Atom_cm_z":0.0})   
syst.CREATE_ATOM(at)

at = lcccsObjects.Atom({"Atom_element":"H","Atom_cm_x":3.0,"Atom_cm_y":-1.0,"Atom_cm_z":0.0}) 
syst.CREATE_ATOM(at)

at = lcccsObjects.Atom({"Atom_element":"Au","Atom_cm_x":0.0,"Atom_cm_y":0.0,"Atom_cm_z":-2.0}) 
syst.CREATE_ATOM(at)

at = lcccsObjects.Atom({"Atom_element":"Au","Atom_cm_x":1.0,"Atom_cm_y":0.0,"Atom_cm_z":-2.0}) 
syst.CREATE_ATOM(at)

at = lcccsObjects.Atom({"Atom_element":"Au","Atom_cm_x":0.0,"Atom_cm_y":1.0,"Atom_cm_z":-2.0}) 
syst.CREATE_ATOM(at)

at = lcccsObjects.Atom({"Atom_element":"Au","Atom_cm_x":1.0,"Atom_cm_y":1.0,"Atom_cm_z":-2.0}) 
syst.CREATE_ATOM(at)

syst.show_atoms()


syst.LINK_ATOMS(1,2)
syst.LINK_ATOMS(1,3)
syst.LINK_ATOMS(1,4)
syst.LINK_ATOMS(4,5)
syst.LINK_ATOMS(4,6)

syst.show_bonds()
syst.show_fragments()
syst.show_molecules()


syst.GROUP_ATOMS([1,2,3],1)
syst.GROUP_ATOMS([4,5,6],2)
syst.GROUP_ATOMS([7,8,9,10],3)

syst.show_fragments()

syst.determine_functional_groups()

syst.show_atoms()


###############################################################

from DataBase.ForceFields.UFF import * #import LoadUFF


##############################################################
# Create force field objects
uff = lcccsObjects.ForceField()

# Load parameters
LoadUFF.Load_UFF(uff)

# Set up functional forms
uff.set_functionals({"bond":"Harmonic","angle":"Harmonic","vdw":"LJ12_6"})


# Define regions for different levels of theory
atlst1 = [1,2,3,4,5,6]
atlst2 = [7,8,9,10]
grlst1 = [1,2]
grlst2 = [3]

# Create interactions between atoms from atom list 1
#syst.set_interactions_for_atoms(atlst1,atlst1,uff)    # The simplest sintax

# Set up interactions in another region
#syst.set_interactions_for_atoms(atlst2,atlst2,uff) 

# Interactions between the regions is set up with following sintax
#syst.set_interactions_for_atoms(atlst1,atlst2,uff) 

