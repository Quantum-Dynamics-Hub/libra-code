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




#set xrange [-0.1:0.1]
#set xtics 1000
#set yrange [-0.06:0.16]
#set ytics 0.04

set xlabel "q" offset 0.0, 0.5
set ylabel "p" offset 1.5, 0.0 
set output "_pop.png"
plot "o" using 1:2   w l  ls 11  lw 3  t "P(0)",\
     "o" using 1:3   w l  ls 21  lw 3  t "P(1)",\
     "o" using 1:4   w l  ls 31  lw 3  t "S",\
     "o" using 1:6   w l  ls 41  lw 3  t "State",\

set xlabel "q" offset 0.0, 0.5
set ylabel "p" offset 1.5, 0.0 
set output "_en.png"
plot "o" using 1:7   w l  ls 11  lw 3  t "Eadi_00",\
     "o" using 1:8   w l  ls 21  lw 3  t "Eadi_11",\
     "o" using 1:10  w l  ls 31  lw 3  t "NAC_01"










