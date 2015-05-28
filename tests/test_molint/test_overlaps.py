import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/qchem")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygqchem import *


print "\nTest 2: 1D overlaps: Unnormalized gaussians"
f = open("1D_overlaps.txt","w")
for i in range(0,50):
    x = 0.1 * i

    #gaussian_overlap(nxa, alp_a, Xa, nxb, alp_b, Xb )
    ss = gaussian_overlap(0, 1.3, 0.0, 0, 1.3, x, 0 )  # s(H)-s(H)
    sp = gaussian_overlap(0, 1.3, 0.0, 1, 1.3, x, 0 )  # s(H)-p(H)
    sd = gaussian_overlap(0, 1.3, 0.0, 2, 1.3, x, 0 )  # s(H)-d(H)
    sf = gaussian_overlap(0, 1.3, 0.0, 3, 1.3, x, 0 )  # s(H)-f(H)
    pp = gaussian_overlap(1, 1.3, 0.0, 1, 1.3, x, 0 )  # p(H)-p(H)
    pd = gaussian_overlap(1, 1.3, 0.0, 2, 1.3, x, 0 )  # p(H)-d(H)
    pf = gaussian_overlap(1, 1.3, 0.0, 3, 1.3, x, 0 )  # p(H)-f(H)
    dd = gaussian_overlap(2, 1.3, 0.0, 2, 1.3, x, 0 )  # d(H)-d(H)
    df = gaussian_overlap(2, 1.3, 0.0, 3, 1.3, x, 0 )  # d(H)-f(H)
    ff = gaussian_overlap(3, 1.3, 0.0, 3, 1.3, x, 0 )  # f(H)-f(H)

    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )

f.close()


print "\nTest 3: 1D overlaps: Normalized gaussians"
print "These are actually normalization coefficients"
s_norm = gaussian_norm(0, 1.3)
p_norm = gaussian_norm(1, 1.3)
d_norm = gaussian_norm(2, 1.3)
f_norm = gaussian_norm(3, 1.3)
print "s_norm = ", s_norm, "|<s|s>|^2 = ", gaussian_overlap(0, 1.3, 0.0,   0, 1.3, 0.0, 0 )*(s_norm**2)
print "p_norm = ", p_norm, "|<p|p>|^2 = ", gaussian_overlap(1, 1.3, 0.0,   1, 1.3, 0.0, 0 )*(p_norm**2)
print "d_norm = ", d_norm, "|<d|d>|^2 = ", gaussian_overlap(2, 1.3, 0.0,   2, 1.3, 0.0, 0 )*(d_norm**2)
print "f_norm = ", f_norm, "|<f|f>|^2 = ", gaussian_overlap(3, 1.3, 0.0,   3, 1.3, 0.0, 0 )*(f_norm**2)

f = open("1D_overlaps_norm.txt","w")
for i in range(0,50):
    x = 0.1 * i

    ss = gaussian_overlap(0, 1.3, 0.0, 0, 1.3, x, 0 )  # s(H)-s(H)
    sp = gaussian_overlap(0, 1.3, 0.0, 1, 1.3, x, 0 )  # s(H)-p(H)
    sd = gaussian_overlap(0, 1.3, 0.0, 2, 1.3, x, 0 )  # s(H)-d(H)
    sf = gaussian_overlap(0, 1.3, 0.0, 3, 1.3, x, 0 )  # s(H)-f(H)
    pp = gaussian_overlap(1, 1.3, 0.0, 1, 1.3, x, 0 )  # p(H)-p(H)
    pd = gaussian_overlap(1, 1.3, 0.0, 2, 1.3, x, 0 )  # p(H)-d(H)
    pf = gaussian_overlap(1, 1.3, 0.0, 3, 1.3, x, 0 )  # p(H)-f(H)
    dd = gaussian_overlap(2, 1.3, 0.0, 2, 1.3, x, 0 )  # d(H)-d(H)
    df = gaussian_overlap(2, 1.3, 0.0, 3, 1.3, x, 0 )  # d(H)-f(H)
    ff = gaussian_overlap(3, 1.3, 0.0, 3, 1.3, x, 0 )  # f(H)-f(H)

    # Normalization:
    ss = s_norm * s_norm * ss
    sp = s_norm * p_norm * sp
    sd = s_norm * d_norm * sd
    sf = s_norm * f_norm * sf
    pp = p_norm * p_norm * pp
    pd = p_norm * d_norm * pd
    pf = p_norm * f_norm * pf
    dd = d_norm * d_norm * dd
    df = d_norm * f_norm * df
    ff = f_norm * f_norm * ff

    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )


f.close()



print "\nTest 4: 3D Gaussians (normalized)"
f = open("3D_overlaps_norm.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(1.0, 0.0, 0.0)

for i in range(0,50):
    Rb.x = 0.1 * i

    ss  = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb )  # s(H)-s(H)
    spx = gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # s(H)-px(H)
    spy = gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # s(H)-py(H)
    spz = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rb )  # s(H)-pz(H)
    pxpx = gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # px(H)-px(H) 
    pxpy = gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # px(H)-py(H) = 0
    pypz = gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rb )  # py(H)-pz(H) = 0
    sdz2 = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rb )  # s(H)-dz2(H) 
    sdxy = gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rb )  # s(H)-dxy(H) 
    pxdxy = gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rb )  # s(H)-dxy(H) 


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (Rb.x, ss, spx, spy, spz, pxpx, pxpy, pypz, sdz2, sdxy, pxdxy) )


f.close()


print "\nTest 5: 1D moments (normalized)"
print "Normalization coefficients"
s_norm = gaussian_norm(0, 1.3)
p_norm = gaussian_norm(1, 1.3)
d_norm = gaussian_norm(2, 1.3)
f_norm = gaussian_norm(3, 1.3)
print "s_norm = ", s_norm, "|<s|s>|^2 = ", gaussian_overlap(0, 1.3, 0.0,   0, 1.3, 0.0, 0 )*(s_norm**2)
print "p_norm = ", p_norm, "|<p|p>|^2 = ", gaussian_overlap(1, 1.3, 0.0,   1, 1.3, 0.0, 0 )*(p_norm**2)
print "d_norm = ", d_norm, "|<d|d>|^2 = ", gaussian_overlap(2, 1.3, 0.0,   2, 1.3, 0.0, 0 )*(d_norm**2)
print "f_norm = ", f_norm, "|<f|f>|^2 = ", gaussian_overlap(3, 1.3, 0.0,   3, 1.3, 0.0, 0 )*(f_norm**2)

f = open("1D_moments_norm.txt","w")
for i in range(0,50):
    x = 0.1 * i

    ss1 = gaussian_moment(0, 0.0, 0.0,  0, 1.3, 0.0, 0, 1.3, x )  # <s(A)| 1 | s(B)>
    ss2 = gaussian_moment(1, 0.0, 0.0,  0, 1.3, 0.0, 0, 1.3, x )  # <s(A)| (x-x(A))| s(B) >
    ss3 = gaussian_moment(1, 0.0, 0.5*x,0, 1.3, 0.0, 0, 1.3, x )  # <s(A)| (x-(x(A)+x(B))/2) | s(B)>
    ss4 = gaussian_moment(2, 0.0, 0.0,  0, 1.3, 0.0, 0, 1.3, x )  # <s(A)| (x-X(A))^2 | s(B)>

    pp1 = gaussian_moment(0, 0.0, 0.0,  1, 1.3, 0.0, 1, 1.3, x )  # <p(A)| 1 | p(B)>
    pp2 = gaussian_moment(1, 0.0, 0.0,  1, 1.3, 0.0, 1, 1.3, x )  # <p(A)| (x-x(A))| s(B) >
    pp3 = gaussian_moment(1, 0.0, 0.5*x,1, 1.3, 0.0, 1, 1.3, x )  # <p(A)| (x-(x(A)+x(B))/2) | p(B)>
    pp4 = gaussian_moment(2, 0.0, 0.0,  1, 1.3, 0.0, 1, 1.3, x )  # <p(A)| (x-X(A))^2 | p(B)>

    sp1 = gaussian_moment(0, 0.0, 0.0,  0, 1.3, 0.0, 1, 1.3, x )  # <s(A)| 1 | p(B)>
    sp2 = gaussian_moment(1, 0.0, 0.0,  0, 1.3, 0.0, 1, 1.3, x )  # <s(A)| (x-x(A))| s(B) >
    sp3 = gaussian_moment(1, 0.0, 0.5*x,0, 1.3, 0.0, 1, 1.3, x )  # <s(A)| (x-(x(A)+x(B))/2) | p(B)>
    sp4 = gaussian_moment(2, 0.0, 0.0,  0, 1.3, 0.0, 1, 1.3, x )  # <s(A)| (x-X(A))^2 | p(B)>



    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss1, ss2, ss3, ss4, pp1, pp2, pp3, pp4, sp1, sp2, sp3, sp4) )


f.close()



print "\nTest 6: Pseudopotentials with 3D Gaussians (normalized)"
f = open("3D_pseudopot02.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(1.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

for i in range(0,50):
    Rb.x = 0.1 * i

    ss_pp   = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb )  # s(H)-s(H)
    spx_pp  = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # s(H)-px(H)
    spy_pp  = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # s(H)-py(H)
    pxpx_pp = pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # px(H)-px(H)
    pxpy_pp = pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # px(H)-py(H)


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (Rb.x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )


f.close()




print "\nTest 7: 1D overlaps: normalized STOs"  # similar to Test 3
print "These are actually normalization coefficients"
Rcut = 12.0
s_norm = sto_norm(0, 1.3)
p_norm = sto_norm(1, 1.3)
d_norm = sto_norm(2, 1.3)
f_norm = sto_norm(3, 1.3)
print "sto_overlap  returns overlaps assuming STOs are normalized already!"
print "s_norm = ", s_norm, "|<s|s>|^2 = ", sto_overlap(0, 0, 0, 1.3,  0, 0, 0, 1.3,  0.0, Rcut )#*(s_norm**2)
print "p_norm = ", p_norm, "|<p|p>|^2 = ", sto_overlap(1, 0, 0, 1.3,  1, 0, 0, 1.3,  0.0, Rcut )#*(p_norm**2)
print "d_norm = ", d_norm, "|<d|d>|^2 = ", sto_overlap(2, 0, 0, 1.3,  2, 0, 0, 1.3,  0.0, Rcut )#*(d_norm**2)
print "f_norm = ", f_norm, "|<f|f>|^2 = ", sto_overlap(3, 0, 0, 1.3,  3, 0, 0, 1.3,  0.0, Rcut )#*(f_norm**2)

f = open("1D_sto_overlaps_norm.txt","w")
for i in range(0,50):
    x = 0.1 * i

    ss = sto_overlap(0, 0, 0, 1.3,  0, 0, 0, 1.3,  x, Rcut )
    sp = sto_overlap(0, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    sd = sto_overlap(0, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    sf = sto_overlap(0, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    pp = sto_overlap(1, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    pd = sto_overlap(1, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    pf = sto_overlap(1, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    dd = sto_overlap(2, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    df = sto_overlap(2, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    ff = sto_overlap(3, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )


f.close()


f = open("1D_sto_overlaps_norm1.txt","w")
for i in range(0,50):
    x = 0.1 * i

    ss = sto_overlap_fast(0, 0, 0, 1.3,  0, 0, 0, 1.3,  x, Rcut )
    sp = sto_overlap_fast(0, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    sd = sto_overlap_fast(0, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    sf = sto_overlap_fast(0, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    pp = sto_overlap_fast(1, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    pd = sto_overlap_fast(1, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    pf = sto_overlap_fast(1, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    dd = sto_overlap_fast(2, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    df = sto_overlap_fast(2, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    ff = sto_overlap_fast(3, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )


f.close()





