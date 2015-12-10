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





set xlabel "t" offset 0.0, 0.5
set ylabel "Energy" offset 1.5, 0.0 
set xrange [0:500]
set xtics 100
set yrange [0:8]
set ytics 1.0
set output "_nucl_en.png"
plot "_nucl.txt" using 2:8   w l  ls 11  lw 5  t "Ekin",\
     "_nucl.txt" using 2:10  w l  ls 21  lw 5  t "Epot",\
     "_nucl.txt" using 2:12  w l  ls 32  lw 5  t "Etot"


set xlabel "t" offset 0.0, 0.5
set ylabel "Energy" offset 1.5, 0.0 
set xrange [0:500]
set xtics 100
set yrange [-5:5]
set ytics 1.0
set output "_nucl2_en.png"
plot "_nucl2.txt" using 2:8   w l  ls 11  lw 5  t "Ekin",\
     "_nucl2.txt" using 2:10  w l  ls 21  lw 5  t "Epot",\
     "_nucl2.txt" using 2:12  w l  ls 32  lw 5  t "Etot"



set xlabel "q" offset 0.0, 0.5
set ylabel "p" offset 1.5, 0.0 
set xrange [-1.5:1.5]
set xtics 0.5
set yrange [-150:150]
set ytics 50
set output "_nucl_phase.png"
plot "_nucl.txt" using 4:6  w l  ls 11  lw 5  notitle "Phase space"


set xlabel "q" offset 0.0, 0.5
set ylabel "p" offset 1.5, 0.0 
set xrange [-0.6:1.2]
set xtics 0.4
set yrange [-100:100]
set ytics 40
set output "_nucl2_phase.png"
plot "_nucl2.txt" using 4:6  w l  ls 11  lw 5  notitle "Phase space"

