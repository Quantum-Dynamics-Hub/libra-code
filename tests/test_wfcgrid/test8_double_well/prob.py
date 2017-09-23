import sys

N = 4000
c = open("prob.txt","w")

for i in range(N):
    sum_all = 0.0
    sum_right = 0.0
    prob_right = 0.0
    a = open("res/wfc.state0.frame"+str(i),"r")
    for line in a:
        b = line.strip().split()
        sum_all += float(b[1]) 
        if float(b[0]) > 0.0:
            sum_right += float(b[1])
    prob_right = sum_right / sum_all
    c.write( str(i*250) + " " +  str(prob_right) + "\n" )
  
