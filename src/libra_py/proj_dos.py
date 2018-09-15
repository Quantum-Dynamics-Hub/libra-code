import math
import os
import sys

Ha_to_eV = 1.0/0.036749309


def fermi_integral(bnds, ef, degen):
#// Compute integral(sum):
#// sum [degen / (1 + exp(e-ef))]  = N
#//  i
#// where N - is a number of electrons (valence), 2 - is because we use closed shell formulation - each orbital represents 2 spin-orbitals
#// so 2 - is a degeneracy of the energy level

    kT = 0.00095; # kT in a.u. for T = 300 K
  
    Norb = len(bnds)
  
    summ = 0.0
    for i in range(0,Norb):
        argg = ((bnds[i]-ef)/kT)
        pop = degen
        if(argg>50.0):
            pass
        elif(argg<-50.0):
            summ += degen
        else:
          pop = degen/(1.0 + math.exp(argg))
          summ += pop
  
    return summ



def fermi_energy(bnds, Nel, degen):
#// Computes Fermi energy

  tol = 1e-5
  err = 2.0*tol

  Norb = len(bnds)

  ef_l = bnds[0] - 10.0
  ef_r = bnds[Norb-1] + 10.0

  i_l = fermi_integral(bnds,ef_l,degen) - Nel
  i_r = fermi_integral(bnds,ef_r,degen) - Nel

#  ef_l,ef_m,ef_r = 0.0,0.0,0.0
#  i_l,i_m,i_r = 0.0,0.0,0.0

  state = 1

  while state==1:
      ef_m = 0.5*(ef_l + ef_r)
      i_m = fermi_integral(bnds,ef_m,degen) - Nel
      
      if(i_m*i_r<=0.0 and i_m*i_l>=0.0):
          i_l = i_m
          ef_l = ef_m
      elif(i_m*i_r>=0.0 and i_m*i_l<=0.0):
          i_r = i_m
          ef_r = ef_m
      else:
          print "Error in fermi_energy\n"
          exit(0) 
      
      err = 0.5*(i_r - i_l)
      
      if err<tol:
          state = 0
      
 
  return 0.5*(ef_l+ef_r)



def pdos(emin, emax, de, projections, prefix, outfile, Nel):
# projections - groups of atoms and types of projections
# e.g. projections = [["s",[1,2,3]], ["p",[1,2,3]], ... ]
# prefix - points to the files containing projection information
# Nel - number of electrons - to compute Fermi energy

#    f = open(prefix+str()"alp","r")
#    A = f.readlines()
#    f.close()


    nproj = len(projections)  # number of projections
    
    # Determine dimensionality and prepare arrays
    en0 = []
    en = []
    dos = []  # dos[i][proj] - dos for level i projected on projection proj

    N = int(math.floor((emax - emin)/de))
    for i in range(0,N+1):
        en.append(emin + i*de)
        ds = []
        for proj in projections:
            ds.append(0.0) 
        dos.append(ds)


    
    for proj in projections:  # loop over all projection
        ang_mom = proj[0]
        atoms = proj[1]

        proj_indx = projections.index(proj)

        for a in atoms: # open files for all atoms in given group
            fa = open(prefix+str(a),"r")
            B = fa.readlines()
            fa.close()

            for lin in B[1:-4]:  # read all lines
                tmp = lin.split()
                 
                e = float(tmp[0])
                if a==0:
                    en0.append(e)  # for Fermi energy calculations

                e = e*Ha_to_eV  # convert to eV

                x = 0.0
                if ang_mom=="s":
                    x = float(tmp[2])
                elif ang_mom=="p":
                    x = float(tmp[3])
                elif ang_mom=="d":
                    x = float(tmp[4])
                elif ang_mom=="tot":
                    x = float(tmp[1])
                else:
                    x = 0.0


                if e<emin or e>emax:
                    pass
                else:
                    state_indx = int(math.floor((e - emin)/de))

#                    if x>0.0:
#                        print state_indx, proj_indx, x

                    dos[state_indx][proj_indx] = dos[state_indx][proj_indx] + x


    e_f = fermi_energy(en0,Nel,2.0)*Ha_to_eV



    f2 = open(outfile,"w")
    f2.write("Ef = %5.3f eV\n" % e_f)
    f2.close()

    # Now compute projections
    for i in range(0,N+1):  # loop over orbitals

        line = str(en[i])+"   "

        tot = 0.0
        for j in range(0,nproj):
            tot = tot + dos[i][j]
            line = line + str(dos[i][j])+"   "

        line = line + str(tot)+"\n"


        f2 = open(outfile,"a")
        f2.write(line)
        f2.close()





# Possible projections        
#proj = [["s",range(0,360)],["p",range(0,360)],["d",range(0,360)]]
#proj = [["s",range(0,1)],["p",range(0,1)],["d",range(0,1)]]
#proj = [["tot",range(0,112)]]


# Example of call - for Si QD
# Si
#main(-35.0, 35.0, 0.1,[["tot",range(0,103)]],"_alpha_wfc_atom","dos_proj.txt",238)    






def convolve(X, dx0, dx, var):

    mult = int(dx0/dx)  # making grid mult times bigger
    print "multiplication factor is = ", mult
#    var = var0/float(mult)
    print "original grid spacing = ", dx0
    print "new grid spacing = ", dx
    print "gaussian variance = ", var
       
 
    # Prepare arrays
    T = len(X)     # how many original grid points
    TM = T*mult    # how many new grid points
    S = len(X[0])  # how many components

    R = []         # result goes here
    t = 0
    while t<TM:
        r = []
        s = 0
        while s<S:
            r.append(0.0)
            s = s + 1
        R.append(r)
        t = t + 1

    # Now do the calculations for all components
    
    prefac = 1.0/(var*math.sqrt(2.0*math.pi)*float(mult))
    alp = 0.5/(var*var)

    t = 0
    while t<TM:
        R[t][0] = X[0][0] + dx*t
        t = t + 1
   
    s = 1
    while s<S:
        t = 0
        while t<TM:
            t1 = 0
            while t1<T:
                w = prefac*math.exp(-alp*(dx*t-dx0*t1)*(dx*t-dx0*t1))  # normalized Gaussian
                R[t][s] = R[t][s] +  X[t1][s] * w    # weight all points of the grid
                t1 = t1 + 1
            t = t + 1
        s = s + 1

    return TM,list(R)


def printout(out,TM,X,R):
# Now print out the sum to file

    f = open(out,'w')
    t = 0
    while t<TM:
        line = ""
        x = 0
        while x<X:
            line = line + str(R[t][x])+"  "
            x = x + 1

        if X==3:  # spin-polarized case
            line = line + str(R[t][1]+R[t][2])+"   "+str(R[t][1]-R[t][2])        
        line = line + "\n"
        f.write(line)
        t = t + 1
    
    f.close()



def sum_pdos(prefix,lst,nat,symb_lst,out,E_f,nspin,do_convolve,dx0,dx,var,PT):
# prefix - a common part of the files containg raw pdos data (including the directory name)
# lst - list of the atoms to sum
# nat - total number of atoms in the system
# symb_lst - symbolyc meaning of the wfc number
# out - file for the output 
# Fermi energy (or HOMO) - center of coordinate system
# nspin - spin-polarization, =1 - non-polarized, =2 - polarized
# var - variance of the Gaussian to convolve with
# dx0  - initial grid spacing
# dx - final grid spacing
# PT - periodic table - a list of the atom types to look for

# Determine the size of the files (number and coordinates of x grid points)
    stat = 0
    XP = []
    scl = len(lst)/float(nat) # fraction of the atoms in given subset from the total number of atoms in the system

    for j in lst:
        for symb in symb_lst:
            for wfc in range(0,5): # Specify max wfc type index - usually no more than 3, 5 - should be more than enough  
                for Elt in PT:

                    filename = prefix+str(j)+"("+Elt+")_wfc#"+str(wfc)+"("+symb+")"  # first file
            
                    if(os.path.exists(filename) and stat==0):
                        print "First file is ", filename
            
                        f = open(filename,'r')
                        A = f.readlines()
                        f.close()

                        start_line = 5
                        if symb=="s":
                            start_line = nspin*1 + 3
                        elif symb=="p":
                            start_line = nspin*3 + 3
                        elif symb=="d":
                            start_line = nspin*5 + 3
                    
                        for a in A[start_line:]:
                            tmp = a.split()
                            line = []
                            i = 0
                            while i<len(tmp):
                                if i==0:
                                    line.append( float(tmp[0]) - E_f)
                                else:
                                    line.append(0.0)
                                i = i + 1
            
                            XP.append( line )
            
                        stat = 1

    # Dimension of the file matrix
    T = len(XP)
    X = nspin + 1  # Only total (for unpolarized) and total up + total down (for spin-polarized) 

    print "Dimensions are: ", T, "by", X


# Read all files from the list
    count = 0.0

    for i in lst:
        for symb in symb_lst:
            for wfc in range(0,5): # Specify max wfc type index - usually no more than 3, 5 - should be more than enough  
                for Elt in PT:

                    filename = prefix+str(i)+"("+Elt+")_wfc#"+str(wfc)+"("+symb+")"  # file
                    
                    if os.path.exists(filename):

                        print "using file", filename
            
                        f = open(filename,'r')
                        A = f.readlines()
                        f.close()

                        start_line = 5
                        if symb=="s":
                            start_line = nspin*1 + 3
                        elif symb=="p":
                            start_line = nspin*3 + 3
                        elif symb=="d":
                            start_line = nspin*5 + 3

            
                        # Summ up all densities
                        t = 0
                        while t<T:
                            x = 1
                            tmp = A[t+start_line].split()
                            while x<X:
                                XP[t][x] = XP[t][x] + scl*float(tmp[x])
                                x = x + 1                           
                                
                            t = t + 1
                    
                        count = count + 1.0



    if do_convolve==1:
        # Convolve XP with gaussian
        TM,R = convolve(XP,dx0,dx,var)
        printout(out,TM,X,R)

    else:
        printout(out,T,X,XP)



# Example of usage:
# PT = ["H","C"]
#E_f = 4.5446  # Fermi energy 
#P1 = range(1,37)
#P2 = range(37,73)
#nat = 72
#sum_pdos('pdos/x0.pdos_atm#',P1,nat,["s","p"],"P1", E_f,1, 1, 0.1, 0.02, 0.1,PT)
#sum_pdos('pdos/x0.pdos_atm#',P2,nat,["s","p"],"P2", E_f,1, 1, 0.1, 0.02, 0.1,PT)

