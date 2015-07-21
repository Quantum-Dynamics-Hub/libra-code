#




set output sprintf("relax%d.png",counter)

#set fit logfile sprintf("relax%d.log",counter)
#set xrange [0:2000000]
#f1(x) = exp(-A1*x) + B1       # define the function to be fit
#A1 = 0.00001; B1 = 0.0;       # initial guess for a1 and b1
#f1(x) = -a*x + b
#a = 0.0001;  b = 0.00;
#fit f1(x) sprintf("relax%d.txt", counter) using $1:$2 via a, b

plot sprintf("relax%d.txt",counter) using  1:2   w l  ls 31  lw 3   t "Population",\
     sprintf("model%d",counter) using  1:2   w l  ls 11  lw 3   t "fit"


unset output
print sprintf('frame%d.png written.',counter)

counter = counter + 1

if(counter<17) reread
