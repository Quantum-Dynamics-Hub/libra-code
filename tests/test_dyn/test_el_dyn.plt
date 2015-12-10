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





set xlabel "x" offset 0.0, 0.5
set ylabel "Population" offset 1.5, 0.0 

# SAC, adiabatic Hamiltonian
set xrange [-10:10]
set xtics 2
set yrange [-0.2:1.5]
set ytics 0.4

set output "dyn_sac.png"
plot "dyn_sac.txt" using 2:3  w l  ls 11  lw 5  t "P(0)",\
     "dyn_sac.txt" using 2:5  w l  ls 21  lw 5  t "Re(rho(0,1))",\
     "dyn_sac.txt" using 2:6  w l  ls 31  lw 5  t "Im(rho(0,1))"


# Rabi2, diabatic Hamiltonian
set xrange [-10:90]
set xtics 20
set yrange [-0.6:1.5]
set ytics 0.4

set output "dyn_rabi2.png"
plot "dyn_rabi2.txt" using 2:3  w l  ls 11  lw 5  t "P(0)",\
     "dyn_rabi2.txt" using 2:5  w l  ls 21  lw 5  t "Re(rho(0,1))",\
     "dyn_rabi2.txt" using 2:6  w l  ls 31  lw 5  t "Im(rho(0,1))"

