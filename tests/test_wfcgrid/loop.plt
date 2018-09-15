#

counter = counter + 1

#set size ratio -1

set terminal pngcairo font "arial,24" size 800, 600 enhanced rounded truecolor

set lmargin at screen 0.17
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.95

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



#set xlabel "Coordinate, a.u."
#set ylabel "Population density"

set xlabel "Coordinate, a.u." offset 0.0, 0.5
set ylabel "Population density" offset 1.5, 0.0 

#set zlabel "Amplitude"

set xrange [-20:20]
set xtics 5
set ytics 0.1

#set nokey
#set palette rgbformulae 22,13,-31

# 3D pictures, colored by z-value 
#set auto
#set parametric
#set pm3d implicit at s

# 2D pictures, colored by z-value
#set contour
#set pm3d map
#set palette

#plot sprintf('wfc.state1.frame%d',counter) t "frame"  w l  lt 1  lw 3  # for semiclassical

set yrange [0.0: 0.5]
#set ytics 0.025
set output sprintf('frame%d.png',counter)

plot sprintf('wfc.state0.frame%d',counter) using 1:2  t "state 0"  w l  ls 12  lw 5,\
     sprintf('wfc.state1.frame%d',counter) using 1:2  t "state 1"  w l  ls 32  lw 5



#set yrange [0.0: 0.5]
#set ytics 0.025
set output sprintf('reci_frame%d.png',counter)

plot sprintf('reci_wfc.state0.frame%d',counter) using 1:2  t "state 0"  w l  lt 1  lw 3,\
     sprintf('reci_wfc.state1.frame%d',counter) using 1:2  t "state 1"  w l  lt 2  lw 3






unset output
print sprintf('frame%d.png written.',counter)

if(counter<200) reread
