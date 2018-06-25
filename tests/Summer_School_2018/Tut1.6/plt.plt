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


set output "_acf_vel.png"
plot "_acf_vel.txt" using 1:2 w l ls 11 lw 5 t "velocity ACF"


set lmargin at screen 0.20
set output "_rho_q.png"
set xlabel "Q (a.u)" offset 0.0, 0.5
set ylabel "\Rho_q" offset 2.2, 0.0
plot  "_density_q.txt" u ($1):2 w l lt 7 lw 5 notitle

set output "_rho_q-first.png"
set xlabel "Q (a.u)" offset 0.0, 0.5
set ylabel "\Rho_q" offset 2.2, 0.0
plot  "_density_q-first.txt" u ($1):2 w l lt 7 lw 5 notitle


set lmargin at screen 0.20
set output "_rho_p.png"
set xlabel "P (a.u)" offset 0.0, 0.5
set ylabel "\Rho_p" offset 2.2, 0.0
plot  "_density_p.txt" u ($1):2 w l lt 7 lw 5 notitle

set output "_rho_p-first.png"
set xlabel "P (a.u)" offset 0.0, 0.5
set ylabel "\Rho_p" offset 2.2, 0.0
plot  "_density_p-first.txt" u ($1):2 w l lt 7 lw 5 notitle



set lmargin at screen 0.20

set output "_1D_verlet_ens_energy.png"
set xlabel "Time (a.u)" offset 0.0, 0.5
set ylabel "Energy (a.u)" offset 2.2, 0.0
plot  "_output.txt" u ($1):3 w l lt 7 lw 5 t "Potential",\
      "_output.txt" u ($1):2 w l lt 6 lw 5 t "Kinetic",\
      "_output.txt" u ($1):4 w l lt 8 lw 5 t "Total",\
      "_output.txt" u ($1):6 w l lt 9 lw 5 t "Extended"

set output "_1D_verlet_ens_T.png"
set xlabel "Time (a.u)" offset 0.0, 0.5
set ylabel "Temperature, K" offset 2.2, 0.0
plot  "_output.txt" u ($1):5 w l lt 7 lw 5 notitle





set lmargin at screen 0.17

set output "_1D_verlet_ens_t-q.png"
set xlabel "Time (a.u)" offset 0.0, 0.5
set ylabel "Position (a.u)" offset 1.3, 0.0
set yrange[-0.75:0.75]
set ytics 0.5
fName = '_pos_space.txt'
stat fName nooutput
N = STATS_columns #number of columns found in file
plot for [i=2:N] fName u ($1):i w l lt 8 lw 1 t "",\


set output "_1D_verlet_ens_t-q-first.png"
set xlabel "Time (a.u)" offset 0.0, 0.5
set ylabel "Position (a.u)" offset 1.3, 0.0
set yrange[-0.75:0.75]
set ytics 0.5
fName = '_pos_space.txt'
plot fName u ($1):2 w l lt 8 lw 1 t ""





set output "_1D_verlet_ens_q-p.png"
set xlabel "Position (a.u)" offset 0.0, 0.5
set ylabel "Momentum (a.u)" offset 1.3, 0.0
set xtics 0.1
set xrange[-0.5:0.5]
set yrange[-1.5:1.5]
set ytics 1.0
fName = '_phase_space.txt'
stat fName nooutput
N = STATS_columns #number of columns found in file
plot for [i=2:N/2] fName u 2*i-1:2*i w l lt 8 lw 1 t "",\


set output "_1D_verlet_ens_q-p-first.png"
set xlabel "Position (a.u)" offset 0.0, 0.5
set ylabel "Momentum (a.u)" offset 1.3, 0.0
set xtics 0.1
set xrange[-0.5:0.5]
set yrange[-1.5:1.5]
set ytics 1.0
fName = '_phase_space.txt'
plot fName u 3:4 w l lt 8 lw 1 t ""





set xrange [100:5000]
set yrange [0:100000]
#set yrange [0:10000]
set xtics 1000
#set ytics 2000
set xlabel "Frequency, cm^{-1}" offset 0.0, 0.5
set ylabel "Intencity" offset 1.5, 0.0

set output "_spectrum_vel.png"
plot "_spectrum_vel.txt" using ($1):($2*$2) w l ls 21 lw 5 t "velocity ACF spectrum"

   


