#
set terminal pngcairo font "Arial,24" size 800, 600 enhanced

set lmargin at screen 0.17
set rmargin at screen 0.95
set bmargin at screen 0.15
set tmargin at screen 0.95


set xlabel "Time, a.u."      offset 0.0, 0.5
set ylabel "Population"    offset 1.5, 0.0


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


#set xrange [0.015 : 0.03]
#set yrange [0.0   : 1e-6]
#set xtics 0.005
#set ytics 2

set output "energy.png"
counter = 0

plot sprintf("relax%d.txt",counter)  using ($1):($8)         w l   ls  11   lw 3  t "e_{kin}",\
     sprintf("relax%d.txt",counter)  using ($1):($9+0.005)   w l   ls  21   lw 3  t "e_{pot}",\
     sprintf("relax%d.txt",counter)  using ($1):($10+0.005)   w l   ls  31   lw 3  t "e_{tot}",\
     sprintf("relax%d.txt",counter)  using ($1):($11+0.005)   w l   ls  41   lw 3  t "e_{ext}"


set output "position.png"
plot sprintf("relax%d.txt",counter)  using ($1):($4)   w l   ls  11   lw 3  t "t-x",\
     sprintf("relax%d.txt",counter)  using ($1):($5)   w l   ls  21   lw 3  t "t-y"

set output "momentum.png"
plot sprintf("relax%d.txt",counter)  using ($1):($14)   w l   ls  11   lw 3  t "t-p_x",\
     sprintf("relax%d.txt",counter)  using ($1):($15)   w l   ls  21   lw 3  t "t-p_y"

set output "phase_portrait.png"
plot sprintf("relax%d.txt",counter)  using ($4):($14)   w l   ls  11   lw 3  t "x-p_x",\
     sprintf("relax%d.txt",counter)  using ($5):($15)   w l   ls  21   lw 3  t "y-p_y"


set output "coherence.png"
plot sprintf("relax%d.txt",counter)  using ($1):($18)   w l   ls  11   lw 3  t "Re(rho_{01})",\
     sprintf("relax%d.txt",counter)  using ($1):($19)   w l   ls  21   lw 3  t "Im(rho_{01})"


#set xrange [:10]
set yrange [0.0:1.1]
set key top
set xlabel "Time, a.u."      offset 0.0, 0.5
set ylabel "Population on lower state"    offset 1.5, 0.0
set xtics 1250

set output "populations.png"
plot sprintf("relax%d.txt",counter)  using ($1):($2)  w p   ls  21   lw 3  t "FSSH/GFSH/MSSH",\
     "exact.txt"  using ($2*25):($4)   w l   ls  11   lw 3  t "exact",\
     sprintf("relax%d.txt",counter)  using ($1):($20)  w l   ls  31   lw 3  t "analytic solution(unperturbed)"
     


