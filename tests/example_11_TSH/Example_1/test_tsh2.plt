#
set terminal pngcairo font "arial,24" size 800, 600 enhanced rounded truecolor

set lmargin at screen 0.17
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.95

#set xtics 1000.0
#set xrange [0.0:4.5]
#set yrange [-40:40]
#set key spacing 1.0 font ",24"
#set key top horizontal

# color definitions
set style line 11 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 5 # --- red
set style line 12 lc rgb '#FF4500' pt 6 ps 1 lt 1 lw 5 # --- orangered
set style line 13 lc rgb '#B22222' pt 6 ps 1 lt 1 lw 5 # --- firebrick
set style line 14 lc rgb '#DC143C' pt 6 ps 1 lt 1 lw 5 # --- crimson

set style line 21 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 5 # --- green
set style line 22 lc rgb '#006400' pt 6 ps 1 lt 1 lw 5 # --- darkgreen
set style line 23 lc rgb '#228B22' pt 6 ps 1 lt 1 lw 5 # --- forestgreen
set style line 24 lc rgb '#808000' pt 6 ps 1 lt 1 lw 5 # --- olive

set style line 31 lc rgb '#8A2BE2' pt 6 ps 1 lt 1 lw 5 # --- blueviolet
set style line 32 lc rgb '#00008B' pt 6 ps 1 lt 1 lw 5 # --- darkblue

set style line 41 lc rgb '#2F4F4F' pt 6 ps 1 lt 1 lw 5 # --- darkslategray




#set xrange [0:5000]
#set xtics 1000
#set yrange [-0.06:0.16]
#set ytics 0.04

set xlabel "q" offset 0.0, 0.5
set ylabel "p" offset 1.5, 0.0 
set output "_q-p.png"
plot "_0_new.txt" using 2:3   w lp  ls 11  lw 5  t "q-p, dia rep"


set xlabel "t" offset 0.0, 0.5
set ylabel "q" offset 1.5, 0.0 
set output "_t-q.png"
plot "_0_new.txt" using 1:2   w l   ls 11  lw 5  t "t-q, dia rep"

set xlabel "t" offset 0.0, 0.5
set ylabel "p" offset 1.5, 0.0 
set output "_t-p.png"
plot "_0_new.txt" using 1:3   w l   ls 11  lw 5  t "t-p, dia rep"



set xlabel "t" offset 0.0, 0.5
set ylabel "Energy" offset 1.5, 0.0 

set output "_energy.png"
plot "_0_new.txt" using 1:4   w l  ls 11  lw 5  t "Ekin",\
     "_0_new.txt" using 1:5   w l  ls 21  lw 5  t "Epot",\
     "_0_new.txt" using 1:6   w l  ls 31  lw 5  t "Etot"


set xlabel "t" offset 0.0, 0.5
set ylabel "Population" offset 1.5, 0.0 

set output "_pop_adi.png"
plot "_0_new.txt" using 1:10   w l  ls 11  lw 5  t "P(0, adi)"


set output "_pop_sh.png"
plot "_0_new.txt" using 1:14   w p  ls 11  lw 5  t "SH, P(0)",\
     "_0_new.txt" using 1:15   w p  ls 21  lw 5  t "SH, P(1)",\
     "_0_new.txt" using 1:10   w p  ls 31  lw 5  t "SE, P(0)",\
     "_0_new.txt" using 1:11   w p  ls 41  lw 5  t "SE, P(1)"

set output "_ordering.png"
plot "_0_new.txt" using 1:16   w p  ls 11  lw 5  t "SH, P(0)"




#set output "_en_dia_new.png"
#plot "_0_new.txt" using 1:10   w l  ls 11  lw 5  t "Hdia_0, dia rep",\
#     "_0_new.txt" using 1:11   w l  ls 12  lw 5  t "Hdia_1, dia rep",\
#     "_1_new.txt" using 1:10   w l  ls 21  lw 5  t "Hdia_0, adi rep",\
#     "_1_new.txt" using 1:11   w l  ls 22  lw 5  t "Hdia_1, adi rep"

#set output "_en_adi_new.png"
#plot "_0_new.txt" using 1:12   w l  ls 11  lw 5  t "Hadi_0, dia rep",\
#     "_0_new.txt" using 1:13   w l  ls 12  lw 5  t "Hadi_1, dia rep",\
#     "_1_new.txt" using 1:12   w l  ls 21  lw 5  t "Hadi_0, adi rep",\
#     "_1_new.txt" using 1:13   w l  ls 22  lw 5  t "Hadi_1, adi rep"

#set output "_en_nac_new.png"
#plot "_0_new.txt" using 1:14   w l  ls 11  lw 5  t "NAC_01, dia rep",\
#     "_1_new.txt" using 1:14   w l  ls 21  lw 5  t "NAC_01, adi rep",\







