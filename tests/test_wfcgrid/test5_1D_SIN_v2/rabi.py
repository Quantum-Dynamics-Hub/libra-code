import math

E0 =-0.001
E1 = 0.001
V01= 0.001

#E0 = -0.010
#E1 = 0.010
#V01 = 0.005


print "E_adi_1 = ", 0.5*(E0 + E1) - math.sqrt( 0.25*(E0 - E1)**2 + V01**2)
print "E_adi_2 = ", 0.5*(E0 + E1) + math.sqrt( 0.25*(E0 - E1)**2 + V01**2)


A = 0.5/math.sqrt(1.0 + 0.25*((E1-E0)/V01)**2)
OMEGA = math.sqrt(V01**2 + 0.25*(E1-E0)**2)

nsnap = 200
nstep = 250

f = open("relax.txt","w")
f.close()

dt = 0.1
for snap in range(0,nsnap):
    t = (snap+1)*nstep*dt
    pop0_ex = (2.0*A*math.sin(OMEGA*t))**2
    pop1_ex = 1.0 - pop0_ex

    f = open("relax.txt","a")
    f.write("%i  %10.8f  %10.8f\n" % (snap+1, pop0_ex, pop1_ex))    
    f.close()
