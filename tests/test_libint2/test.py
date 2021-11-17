from liblibra_core import *
from libra_py import data_outs
import math

x = libint2_Shell()



expts = Py2Cpp_double([1.0, 0.5, 0.02])
coeff = Py2Cpp_double([1.0, -0.5, 1.2])
r = VECTOR(0.0, 0.0, 0.0)

expts1 = Py2Cpp_double([10.0, 0.15, 0.2])
coeff1 = Py2Cpp_double([0.1, -0.5, 0.2])
r1 = VECTOR(0.1, 0.0, 0.0)



print( type(expts))
print( type(coeff))

my_shells = initialize_shell(  0,  True,  expts, coeff,  r )
print( type(my_shells))

print_shell(my_shells)

add_to_shell( my_shells,    0,  True,  expts1, coeff1,  r1 )

print_shell(my_shells)
                      
nthreads = 4    
S = compute_overlaps(my_shells, my_shells, nthreads)

data_outs.print_matrix(S)


ao1 = AO()
ao2 = AO()
for i in range(3):
    ao1.add_primitive(coeff[i],  PrimitiveG(0,0,0, expts[i], r)   )
    ao2.add_primitive(coeff1[i], PrimitiveG(0,0,0, expts1[i], r1) )

print ( gaussian_overlap(ao1, ao1))
print ( gaussian_overlap(ao1, ao2))
print ( gaussian_overlap(ao2, ao1))
print ( gaussian_overlap(ao2, ao2))


print("==================== L = 1 =====================\n")

my_shells = initialize_shell(  1,  True,  expts, coeff,  r )
add_to_shell( my_shells,    1,  True,  expts1, coeff1,  r1 )

nthreads = 4    
S = compute_overlaps(my_shells, my_shells, nthreads)
print( "AO overlaps computed with the wrapped Libint2 codes\n")
data_outs.print_matrix(S)



ao = AOList()
for i in range(6):
    ao.append( AO() )

for i in range(3):
    ao[0].add_primitive(coeff[i],  PrimitiveG(1,0,0, expts[i], r)   )
    ao[1].add_primitive(coeff[i],  PrimitiveG(0,1,0, expts[i], r)   )
    ao[2].add_primitive(coeff[i],  PrimitiveG(0,0,1, expts[i], r)   )
    ao[3].add_primitive(coeff1[i], PrimitiveG(1,0,0, expts1[i], r1)   )
    ao[4].add_primitive(coeff1[i], PrimitiveG(0,1,0, expts1[i], r1)   )
    ao[5].add_primitive(coeff1[i], PrimitiveG(0,0,1, expts1[i], r1)   )

S = MATRIX(6,6)
for i in range(6):
    for j in range(6):
        S.set(i,j, gaussian_overlap(ao[i], ao[j]) )
print( "AO overlaps computed with internal Libra codes\n")
data_outs.print_matrix(S)



print("==================== L = 2 =====================\n")

my_shells = initialize_shell(  2,  True,  expts, coeff,  r )
add_to_shell( my_shells,    2,  True,  expts1, coeff1,  r1 )

nthreads = 4    
S = compute_overlaps(my_shells, my_shells, nthreads)
print( "AO overlaps computed with the wrapped Libint2 codes\n")
data_outs.print_matrix(S)




ao = AOList()
for i in range(10):
    ao.append( AO() )

sq = math.sqrt(2.0)
for i in range(3):
    ao[0].add_primitive(2.0*coeff[i],  PrimitiveG(0,0,2, expts[i], r)   )  # 2*z^2
    ao[0].add_primitive(-coeff[i],  PrimitiveG(2,0,0, expts[i], r)   )  # -x^2
    ao[0].add_primitive(-coeff[i],  PrimitiveG(0,2,0, expts[i], r)   )  # -y^2
    ao[1].add_primitive( coeff[i],  PrimitiveG(1,1,0, expts[i], r)   )  # xy
    ao[2].add_primitive( coeff[i],  PrimitiveG(1,0,1, expts[i], r)   )  # xz
    ao[3].add_primitive( coeff[i],  PrimitiveG(2,0,0, expts[i], r)   )  # x^2
    ao[3].add_primitive(-coeff[i],  PrimitiveG(0,2,0, expts[i], r)   )  # -y^2
    ao[4].add_primitive( coeff[i],  PrimitiveG(0,1,1, expts[i], r)   )  # yz

    ao[5].add_primitive(2.0*coeff1[i],  PrimitiveG(0,0,2, expts1[i], r1)   )  # 2*z^2
    ao[5].add_primitive(-coeff1[i],  PrimitiveG(2,0,0, expts1[i], r1)   )  # -x^2
    ao[5].add_primitive(-coeff1[i],  PrimitiveG(0,2,0, expts1[i], r1)   )  # -y^2
    ao[6].add_primitive( coeff1[i],  PrimitiveG(1,1,0, expts1[i], r1)   )  # xy
    ao[7].add_primitive( coeff1[i],  PrimitiveG(1,0,1, expts1[i], r1)   )  # xz
    ao[8].add_primitive( coeff1[i],  PrimitiveG(2,0,0, expts1[i], r1)   )  # x^2
    ao[8].add_primitive(-coeff1[i],  PrimitiveG(0,2,0, expts1[i], r1)   )  # -y^2
    ao[9].add_primitive( coeff1[i],  PrimitiveG(0,1,1, expts1[i], r1)   )  # yz



S = MATRIX(10,10)
for i in range(10):
    for j in range(10):
        S.set(i,j, gaussian_overlap(ao[i], ao[j]) )
print( "AO overlaps computed with internal Libra codes\n")
data_outs.print_matrix(S)

