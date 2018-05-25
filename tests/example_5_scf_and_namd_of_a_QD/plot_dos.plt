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



set output "DOS.png"
set xrange [-4.0:4.0]
set yrange [0:5]
set key top left
set xlabel "E-E_{F}, eV"     offset 0.0, 0.5
set ylabel "DOS"             offset 1.5, 0.0
set xtics 1.0
set ytics 2.5

plot "pdos.txt" using ($1-(-22.734)):($2)  w l   ls 42  lw 3  t "total DOS, EHT",\
     "pdos.txt" using ($1-(-22.734)):($3)  w l   ls 12  lw 3  t "s, EHT",\
     "pdos.txt" using ($1-(-22.734)):($4)  w l   ls 22  lw 3  t "p, EHT",\
     "pdos.txt" using ($1-(-22.734)):($5)  w l   ls 32  lw 3  t "d, EHT"

