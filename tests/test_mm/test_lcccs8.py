###################################################################
# Tutorial: Creating State and setting its functioning
# Add barostat!
# This is for developing barostat things
###################################################################

import lcccsObjects
from DataBase.AtomicData.PeriodicTable import *
from DataBase.ForceFields.UFF import * #import LoadUFF
from DataBase.ForceFields.DREIDING import * #import LoadDREIDING
from DataBase.ForceFields.TRIPOS import * # import LoadTRIPOS
from DataBase.ForceFields.GAFF import * #import LoadGAFF
from DataBase.ForceFields.MMFF94 import * #import LoadMMFF94
from DataBase.ForceFields.TIP3P import * #import LoadTIP3P

##############################################################
import os
import LoadMolecule

#======= Set up objects here ==========================
# Universe
U = lcccsObjects.Universe()
LoadPT.Load_PT(U)

# MD
md = lcccsObjects.MD({"max_step":1000,"ensemble":"NVT","integrator":"KLN","terec_exp_size":10})
md.show_info()

# Thermostat
#therm = lcccsObjects.Thermostat({"Temperature":250,"Q":100.0})
therm = lcccsObjects.Thermostat({"Temperature":300,"Q":100.0,"thermostat_type":"Nose-Hoover","nu_therm":0.1,"NHC_size":5})
therm.show_info()

# Barostat
baro = lcccsObjects.Barostat({"W":100.0,"Pressure":1.0,"nu_baro":1.1})
baro.show_info()

# Force fields
#uff = lcccsObjects.ForceField({"bond_functional":"Harmonic","vdw_functional":"LJ","angle_functional":"Fourier"})
#uff = lcccsObjects.ForceField({"bond_functional":"Harmonic","angle_functional":"Fourier","vdw_functional":"LJ","dihedral_functional":"General0"})
uff = lcccsObjects.ForceField({"bond_functional":"Harmonic","angle_functional":"Fourier","vdw_functional":"LJ","dihedral_functional":"General0","R_vdw_on":6.0,"R_vdw_off":7.0})
LoadUFF.Load_UFF(uff)
uff.show_info()

dreiding = lcccsObjects.ForceField({"bond_functional":"Harmonic","angle_functional":"Harmonic_Cos","vdw_functional":"LJ","dihedral_functional":"General2"})
LoadDREIDING.Load_DREIDING(dreiding)
dreiding.show_info()

tripos = lcccsObjects.ForceField({})
LoadTRIPOS.Load_TRIPOS(tripos)
tripos.show_info()

gaff = lcccsObjects.ForceField({})
LoadGAFF.Load_GAFF(gaff)
gaff.show_info()

mmff94 = lcccsObjects.ForceField({})
LoadMMFF94.Load_MMFF94(mmff94)
mmff94.show_info()

#tip3p = lcccsObjects.ForceField({"bond_functional":"Harmonic","angle_functional":"Harmonic","vdw_functional":"LJ","elec_functional":"Coulomb","R_vdw_on":6.0,"R_vdw_off":7.0,"R_elec_on":6.0,"R_elec_off":7.0})
tip3p = lcccsObjects.ForceField({"bond_functional":"Harmonic","angle_functional":"Harmonic","vdw_functional":"LJ","mb_functional":"Ewald_3D","R_vdw_on":6.0,"R_vdw_off":7.0})
#tip3p = lcccsObjects.ForceField({"bond_functional":"Harmonic","angle_functional":"Harmonic","vdw_functional":"LJ","R_vdw_on":6.0,"R_vdw_off":7.0})
#tip3p = lcccsObjects.ForceField({"bond_functional":"Harmonic","angle_functional":"Harmonic","vdw_functional":"LJ","elec_functional":"Coulomb"})
LoadTIP3P.Load_TIP3P(tip3p)
tip3p.show_info()







testff = lcccsObjects.ForceField({"bond_functional":"Harmonic","vdw_functional":"LJ","angle_functional":"Fourier","R_vdw_on":6.0,"R_vdw_off":7.0})
#testff = lcccsObjects.ForceField({"bond_functional":"Harmonic","vdw_functional":"LJ","angle_functional":"Fourier"})
LoadUFF.Load_UFF(testff)
testff.show_info()




#========== Create the system to simulate ===============

i = 1
while i<=1:
#for i in [1]:
    print "=================== System ",i,"======================="
    # System creation section
    syst = lcccsObjects.System()

    # Create atoms, link them, define groups
    #inp_file = os.getcwd()+"/Input/GAFF_Tests/test"+str(i)+"a.pdb"
#    inp_file = os.getcwd()+"/Input/23waters_aa.ent"
#    inp_file = os.getcwd()+"/Input/23waters.ent"
#    inp_file = os.getcwd()+"/Input/Cub_flex.ent"
    inp_file = os.getcwd()+"/Input/Rotors/rotor1_flex.ent"
    #inp_file = os.getcwd()+"/Input/Rings/test"+str(i)+".pdb"
    LoadMolecule.Load_Molecule(U, syst,inp_file,"pdb")

    syst.init_fragments()
    syst.init_molecules()

    syst.show_info()
    syst.show_fragments()
    syst.show_molecules()
    
 
#    j = 0
#    while j<2:
#        syst.GROUP_ATOMS([3*j+1,3*j+2,3*j+3],j+1)
#        j = j + 1
#    syst.show_fragments()

    syst.print_ent(os.getcwd()+"/Run/original.ent")
#    syst.init_box(20.0,20.0,20.0)
    syst.show_info()
    syst.print_ent(os.getcwd()+"/Run/after_box.ent")
 
    #syst.GROUP_ATOMS([1,2,3,4],1)
    #syst.GROUP_ATOMS([5,6,7],2)
    #syst.GROUP_ATOMS([8,9,10,11],3)


    # Define regions for different levels of theory
    atlst1 = range(1,syst.Number_of_atoms+1)

    # Create interactions between atoms from atom list 1
    t1 = os.times()
    syst.set_interactions_for_atoms(atlst1,atlst1,uff)    # The simplest sintax
    t2 = os.times()
    print "Setting force field parameters takes ", t2[0]-t1[0], " s"
#    exit(0)
    syst.show_interactions_statistics()
#    exit(0)
 
#    syst.show_interactions("angle")
#    syst.show_pairs()
#    syst.show_frag_pairs()
#    syst.show_rings()
#    syst.show_atoms()

    syst.zero_atom_forces()
#    exit(0)
    print "system bond energy = ", syst.energy("bond")/627.5094709
    print "system angle energy = ", syst.energy("angle")/627.5094709
#    exit(0)
    print "system dihedral energy = ", syst.energy("dihedral")/627.5094709
    print "system vdw energy = ", syst.energy("vdw")/627.5094709
#    exit(0)
    print "system elec energy = ", syst.energy("elec")/627.5094709
#    exit(0)
    print "system mb energy = ", syst.energy("mb")/627.5094709

    print "system energy = ",syst.energy()/627.5094709

#    exit(0)

    syst.zero_atom_forces()
#    syst.init_fragments()
#    syst.show_fragments()
#    syst.show_info()

 
#    st = lcccsObjects.State(syst)         # Good for NVE ensemble
    st = lcccsObjects.State(syst,therm)   # Good for NVT or NVE ensemble
#    st = lcccsObjects.State(syst,baro)     # Good for NPH or NVE ensemble
#    st = lcccsObjects.State(syst,therm,baro) # Good for NPT ensemble
    st.set_md(md)
    f = open(os.getcwd()+"/Run/output.txt","w")
    f.close()

    j=0
    st.init_md()
#    exit(0)
    t1 = os.times()
    while j<5000:
#        st.init_md()
#        print j
        st.run_md()
        line = " step "+str(j+1)+"  "+str(st.E_kin)+"  "+str(st.E_pot)+"  "+str(st.E_tot)+" curr_T= "+str(st.curr_T)+" H_NP= "+str(st.H_NP)+\
               " s_var = "+str(therm.s_var) + " Ps = "+str(therm.Ps)+" volume = "+str(st.curr_V)+" pressure =  "+str(st.curr_P)+"\n"
        f = open(os.getcwd()+"/Run/output.txt","a")
        f.write(line)
        f.close()
        syst.print_ent(os.getcwd()+"/Run/"+str(j)+".ent")
#        st.cool()
        j = j + 1
    t2 = os.times()
    print "MD simulation takes ", t2[0]-t1[0], " s"
 
    exit(0)
    
#    f = open(os.getcwd()+"/Run/output.txt","w")
#    f.close()

    st.init_md()
    j=0
    while j<5000:
        st.run_md()
        line = " step "+str(j+1)+"  "+str(st.E_kin)+"  "+str(st.E_pot)+"  "+str(st.E_tot)+" curr_T= "+str(st.curr_T)+" H_NP= "+str(st.H_NP)+\
               " s_var = "+str(therm.s_var) + " Ps = "+str(therm.Ps) + "\n"
        f = open(os.getcwd()+"/Run/output.txt","a")
        f.write(line)
        syst.print_ent(os.getcwd()+"/Run/"+str(int((j+1)/10.0))+".ent")
        f.close()
        j = j + 1

    print "========================================================="
    i = i + 1


