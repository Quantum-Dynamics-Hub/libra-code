import math

def Regression(X,Y,sz):
# Finds y = a + b*x
    x,y,x2,y2,xy = 0.0, 0.0, 0.0, 0.0, 0.0

    i = 0
    while i<sz:
        x = x + X[i]
        y = y + Y[i]
        x2 = x2 + X[i]*X[i]
        y2 = y2 + Y[i]*Y[i]
        xy = xy + X[i]*Y[i]
        i = i + 1
    N = float(sz)
    b = (N*xy - x*y)/(N*x2 - x*x)
    a = (y - b*x)/N

    b = y/x
    a = 0.0

    return [a,b]


def main(Tcut, Toff, val_indices, filename):
# val_indices - the numberation starts with 0

    f = open(filename,"r")
    A = f.readlines()
    f.close()

    T, P = [], []
    sz = len(A)
    act_sz = 0
    i = 0
    while i<sz:
        tmp = A[i].split()

        t = float(tmp[1])
        val = 0.0
        for val_indx in val_indices:
            val = val + float(tmp[val_indx])

        if val>0.0001 and t>Tcut and t<Toff:
            T.append((t - Tcut))
            P.append(math.log(val) ) 
            act_sz = act_sz + 1
        i = i + 1

    a,b = Regression(T,P,act_sz)

    tau = -1.0/b
    A = math.exp(a)
    print "tau = ",tau, "sqrt(tau)= ", math.sqrt(tau), "k = ", -b, "A = ", A

    f = open(filename+"_model.txt","w")

    for t in T:
        p = A * math.exp(-t/tau)

        line = str(t + Tcut)+"   "+str(p)+"\n"
        f.write(line)
    f.close()



main(0.0, 1500.0, [9,10], "_populations.txt")
main(0.0, 1500.0, [9,10], "_populations_sh.txt")

    
