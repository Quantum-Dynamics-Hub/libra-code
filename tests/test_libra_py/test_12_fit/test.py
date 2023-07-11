import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


X, Y = fit.get_data_from_file("relax.txt", 0, 1)
Ypred, A, B = fit.fit_exp(X,Y, 0.0)  # Y = A * exp(-B*X)
#Ypred, A, B = fit.fit_gau(X,Y, 0.0)  # Y = A * exp(-B*X^2)


f = open("relax_and_fit.txt","w")
for i in xrange(len(X)):
    f.write("%8.5f  %8.5f  %8.5f\n" % ( X[i],Y[i],Ypred[i]) )
f.close()

