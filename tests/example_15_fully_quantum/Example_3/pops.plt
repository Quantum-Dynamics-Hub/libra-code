#
set terminal pngcairo font "arial,24" size 800, 600 enhanced rounded truecolor

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

set lmargin at screen 0.20

set output "_pops.png"
set xlabel "Time, fs" offset 0.0, 0.5
set ylabel "Transmission probability" offset 2.2, 0.0
plot  "_pops.txt" u ($1/41.0):5 w l ls 22 lw 5 t "Quantum"


set output "_energy.png"
set xlabel "Time, fs" offset 0.0, 0.5
set ylabel "Energy, a.u." offset 2.2, 0.0
plot  "_pops.txt" u ($1/41.0):2 w l ls 21 lw 5 t "Quantum, E_{kin}",\
      "_pops.txt" u ($1/41.0):3 w l ls 24 lw 5 t "Quantum, E_{pot}",\
      "_pops.txt" u ($1/41.0):4 w l ls 25 lw 5 t "Quantum, E_{tot}"

set output "_q-t.png"
set xlabel "Time, fs" offset 0.0, 0.5
set ylabel "q, a.u." offset 2.2, 0.0
plot  "_pops.txt" u ($1/41.0):6 w l ls 21 lw 5 t "Position"

set output "_p-t.png"
set xlabel "Time, fs" offset 0.0, 0.5
set ylabel "p, a.u." offset 2.2, 0.0
plot  "_pops.txt" u ($1/41.0):7 w l ls 21 lw 5 t "Momentum"

set output "_q-p.png"
set xlabel "q, a.u." offset 0.0, 0.5
set ylabel "p, a.u." offset 2.2, 0.0
plot  "_pops.txt" u 6:7 w l ls 21 lw 5 t "Phase portrait"








