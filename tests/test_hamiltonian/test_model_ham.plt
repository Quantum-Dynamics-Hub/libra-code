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
set ylabel "SAC Ham" offset 1.5, 0.0
set output "sac_ham.png"
plot "sac_ham.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "sac_ham.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "sac_ham.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "sac_ham1.png"
plot "sac_ham1.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "sac_ham1.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "sac_ham1.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "sac_ham2.png"
plot "sac_ham2.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "sac_ham2.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "sac_ham2.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"



set xlabel "x" offset 0.0, 0.5
set ylabel "m-SAC Ham" offset 1.5, 0.0
set output "sac_ham-m.png"
plot "m-sac_ham.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "m-sac_ham.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "m-sac_ham.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "sac_ham1-m.png"
plot "m-sac_ham1.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "m-sac_ham1.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "m-sac_ham1.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "sac_ham2-m.png"
plot "m-sac_ham2.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "m-sac_ham2.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "m-sac_ham2.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"




set xlabel "x" offset 0.0, 0.5
set ylabel "DAC Ham" offset 1.5, 0.0
set output "dac_ham.png"
plot "dac_ham.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "dac_ham.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "dac_ham.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "dac_ham1.png"
plot "dac_ham1.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "dac_ham1.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "dac_ham1.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "dac_ham2.png"
plot "dac_ham2.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "dac_ham2.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "dac_ham2.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"




set xlabel "x" offset 0.0, 0.5
set ylabel "ECWR Ham" offset 1.5, 0.0
set output "ecwr_ham.png"
plot "ecwr_ham.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "ecwr_ham.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "ecwr_ham.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "ecwr_ham1.png"
plot "ecwr_ham1.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "ecwr_ham1.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "ecwr_ham1.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "ecwr_ham2.png"
plot "ecwr_ham2.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "ecwr_ham2.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "ecwr_ham2.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"




set xlabel "x" offset 0.0, 0.5
set ylabel "Marcus Ham" offset 1.5, 0.0
set output "marcus_ham.png"
plot "marcus_ham.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "marcus_ham.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "marcus_ham.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "marcus_ham1.png"
plot "marcus_ham1.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "marcus_ham1.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "marcus_ham1.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"

set output "marcus_ham2.png"
plot "marcus_ham2.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "marcus_ham2.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "marcus_ham2.txt" using 1:4  w l  ls 31  lw 5  t "H(0,1)"




set xlabel "x" offset 0.0, 0.5
set ylabel "SEXCH Ham" offset 1.5, 0.0
set output "sexch_ham.png"
plot "sexch_ham.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "sexch_ham.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "sexch_ham.txt" using 1:4  w l  ls 31  lw 5  t "H(2,2)",\
     "sexch_ham.txt" using 1:5  w l  ls 12  lw 5  t "H(0,1)",\
     "sexch_ham.txt" using 1:6  w l  ls 22  lw 5  t "H(0,2)",\
     "sexch_ham.txt" using 1:7  w l  ls 32  lw 5  t "H(1,2)"

set output "sexch_ham1.png"
plot "sexch_ham1.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "sexch_ham1.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "sexch_ham1.txt" using 1:4  w l  ls 31  lw 5  t "H(2,2)",\
     "sexch_ham1.txt" using 1:5  w l  ls 12  lw 5  t "H(0,1)",\
     "sexch_ham1.txt" using 1:6  w l  ls 22  lw 5  t "H(0,2)",\
     "sexch_ham1.txt" using 1:7  w l  ls 32  lw 5  t "H(1,2)"


set output "sexch_ham2.png"
plot "sexch_ham2.txt" using 1:2  w l  ls 11  lw 5  t "H(0,0)",\
     "sexch_ham2.txt" using 1:3  w l  ls 21  lw 5  t "H(1,1)",\
     "sexch_ham2.txt" using 1:4  w l  ls 31  lw 5  t "H(2,2)",\
     "sexch_ham2.txt" using 1:5  w l  ls 12  lw 5  t "H(0,1)",\
     "sexch_ham2.txt" using 1:6  w l  ls 22  lw 5  t "H(0,2)",\
     "sexch_ham2.txt" using 1:7  w l  ls 32  lw 5  t "H(1,2)"




