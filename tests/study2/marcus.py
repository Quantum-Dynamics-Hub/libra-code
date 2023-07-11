from math import sqrt
from math import exp

def rates():
    Er = 2.39e-2
    omega = 3.5e-4
    kT = 9.5e-4
    V = 5.0e-5

    dE = (3.0e-2 - 1.5e-2)/16.0

    for i in range(0,17):
        eps = 1.5e-2 + i * dE

        k = (2.0*3.1415926*V*V/sqrt(4.0*3.1415926*Er*kT))*exp(-(Er-eps)**2/(4.0*Er*kT))

        print "%12.8f  %12.8E "  %(eps, k)

rates()
