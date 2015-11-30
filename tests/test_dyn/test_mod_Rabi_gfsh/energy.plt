#
set terminal pngcairo font "Arial,24" size 800, 600 enhanced

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




counter = 0


set output "energy.png"
set key bottom
set xlabel "Time, a.u."      offset 0.0, 0.5
set ylabel "Energy, a.u."    offset 2.5, 0.0
set xtics 1250
set ytics 0.005
plot sprintf("relax%d.txt",counter)  using ($1):($6)   w l   ls  11   lw 3  t "E_{kin}",\
     sprintf("relax%d.txt",counter)  using ($1):($7)   w l   ls  21   lw 3  t "E_{pot}",\
     sprintf("relax%d.txt",counter)  using ($1):($8)   w l   ls  31   lw 3  t "E_{tot}"



set output "phase_portrait.png"
set key top
set xlabel "q, a.u."      offset 0.0, 0.5
set ylabel "p, a.u."    offset 1.5, 0.0
set xtics 2
set ytics 0.5
#set yrange [-0.001:0.01]
plot sprintf("relax%d.txt",counter)  using ($4):($12)   w l   ls  11   lw 3  t "q-p"



set output "position.png"
set key top             
set xlabel "t, a.u."    offset 0.0, 0.5
set ylabel "q, a.u."    offset 1.5, 0.0
set xtics 1250
set ytics 2
plot sprintf("relax%d.txt",counter)  using ($1):($4)   w l   ls  11   lw 3  t "t-q"



set output "momentum.png"
set key top             
set xlabel "t, a.u."    offset 0.0, 0.5
set ylabel "p, a.u."    offset 1.5, 0.0
set xtics 1250
set ytics 0.5
plot sprintf("relax%d.txt",counter)  using ($1):($12)   w l   ls  11   lw 3  t "t-p"



set output "populations.png"
set yrange [0.0:0.7]
set key top
set xlabel "Time, a.u."                   offset 0.0, 0.5
set ylabel "Population on lower state"    offset 1.5, 0.0
set xtics 1250
set ytics 0.1
plot sprintf("relax%d.txt",counter)  using ($1):(1.0-$11)  w p   ls  21   lw 3  t "FSSH/GFSH/MSSH",\
     "exact.txt"  using ($2*25):($4)   w l   ls  11   lw 3  t "exact",\
     sprintf("relax%d.txt",counter)  using ($1):($16)  w l   ls  31   lw 3  t "analytic solution(unperturbed)"


     
set output "coherence.png"
set yrange [-0.4:0.7]
set key top
set xlabel "Time, a.u."      offset 0.0, 0.5
set ylabel "Coherence"       offset 1.5, 0.0
set xtics 1250
set ytics 0.2
plot sprintf("relax%d.txt",counter)  using ($1):($14)   w l   ls  11   lw 3  t "Real component",\
     sprintf("relax%d.txt",counter)  using ($1):($15)   w l   ls  21   lw 3  t "Imaginary component"


