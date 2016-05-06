import math

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



