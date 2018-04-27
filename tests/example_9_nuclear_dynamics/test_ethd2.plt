#
set terminal pngcairo font "arial,24" size 800, 600 enhanced rounded truecolor

set lmargin at screen 0.17
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.95

# color definitions
set style line 11 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 5 # --- red
set style line 12 lc rgb '#FF0000' pt 1 ps 1 lt 1 lw 5 # --- bright red
set style line 13 lc rgb '#FF4500' pt 6 ps 1 lt 1 lw 5 # --- orangered
set style line 14 lc rgb '#B22222' pt 6 ps 1 lt 1 lw 5 # --- firebrick
set style line 15 lc rgb '#DC143C' pt 6 ps 1 lt 1 lw 5 # --- crimson
set style line 16 lc rgb '#FF7000' pt 6 ps 1 lt 1 lw 5 # --- dark orange

set style line 21 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 5 # --- green
set style line 22 lc rgb '#006400' pt 6 ps 1 lt 1 lw 5 # --- darkgreen
set style line 23 lc rgb '#228B22' pt 6 ps 1 lt 1 lw 5 # --- forestgreen
set style line 24 lc rgb '#808000' pt 6 ps 1 lt 1 lw 5 # --- olive
set style line 25 lc rgb '#00FF10' pt 6 ps 1 lt 1 lw 5 # --- lime
set style line 26 lc rgb '#20B2AA' pt 6 ps 1 lt 1 lw 5 # --- light sea green

set style line 31 lc rgb '#8A2BE2' pt 6 ps 1 lt 1 lw 5 # --- blueviolet
set style line 32 lc rgb '#000E8B' pt 6 ps 1 lt 1 lw 5 # --- royalblue
set style line 33 lc rgb '#00008B' pt 6 ps 1 lt 1 lw 5 # --- darkblue
set style line 34 lc rgb '#800080' pt 6 ps 1 lt 1 lw 5 # --- purple

set style line 41 lc rgb '#2F4F4F' pt 6 ps 1 lt 1 lw 5 # --- darkslategray
set style line 42 lc rgb '#C0C0C0' pt 6 ps 1 lt 1 lw 5 # --- silver
set style line 43 lc rgb '#D2691E' pt 6 ps 1 lt 1 lw 5 # --- chocolate

set output "_tunnel.png"
set xlabel "Time (a.u)"  offset 0.0, 0.5
set ylabel "Population on initial State" offset 1.3, 0.0
plot  "_output.txt" u ($1):5 w l ls 22 lw 10 t "ETHD"


set output "_energy.png"
set xlabel "Time (a.u)" offset 0.0, 0.5
set ylabel "Energy (a.u)" offset 1.7, 0.0
set ytics 0.1
plot  "_output.txt" u ($1):3 w l lt 7 lw 5 t "Potential",\
      "_output.txt" u ($1):2 w l lt 6 lw 5 t "Kinetic",\
      "_output.txt" u ($1):4 w l lt 8 lw 5 t "Total"


set output "_t-q.png"
set xlabel "Time (a.u)" offset 0.0, 0.5
set ylabel "Position (a.u)" offset 1.3, 0.0
set yrange[-2.0:2.0]
set ytics 1.0
fName = '_pos_space.txt'
stat fName nooutput
N = STATS_columns #number of columns found in file
plot for [i=2:N] fName u ($1):i w l lt 8 lw 1 t "",\


set output "_q-p.png"
set xlabel "Position (a.u)" offset 0.0, 0.5
set ylabel "Momentum (a.u)" offset 1.3, 0.0
set xtics 0.1
set xrange[-2.0:2.0]
set ytics 10.0
set yrange[-50.0:50.0]
fName = '_phase_space.txt'
stat fName nooutput
N = STATS_columns #number of columns found in file
plot for [i=2:N/2] fName u 2*i-1:2*i w l lt 8 lw 1 t "",\



