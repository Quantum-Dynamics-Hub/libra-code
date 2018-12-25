#
set terminal pngcairo font "arial,26" size 800, 600 enhanced rounded truecolor

# color definitions
set style line 11 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 5 # --- red
set style line 12 lc rgb '#FF0000' pt 1 ps 1 lt 1 lw 5 # --- bright red
set style line 13 lc rgb '#FF4500' pt 6 ps 1 lt 1 lw 5 # --- orangered
set style line 14 lc rgb '#B22222' pt 6 ps 1 lt 1 lw 5 # --- firebrick
set style line 15 lc rgb '#DC143C' pt 6 ps 1 lt 1 lw 5 # --- crimson
set style line 16 lc rgb '#FF7000' pt 6 ps 1 lt 1 lw 5 # --- dark orange
set style line 17 lc rgb '#FF00FF' pt 6 ps 1 lt 1 lw 5 # --- fuchsia
                                                             

set style line 21 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 5 # --- green
set style line 22 lc rgb '#006400' pt 6 ps 1 lt 1 lw 5 # --- darkgreen
set style line 23 lc rgb '#228B22' pt 6 ps 1 lt 1 lw 5 # --- forestgreen
set style line 24 lc rgb '#808000' pt 6 ps 1 lt 1 lw 5 # --- olive
set style line 25 lc rgb '#00FF10' pt 6 ps 1 lt 1 lw 5 # --- lime
set style line 26 lc rgb '#20B2AA' pt 6 ps 1 lt 1 lw 5 # --- light sea green
set style line 27 lc rgb '#00FF7F' pt 6 ps 1 lt 1 lw 5 # --- spring green
set style line 28 lc rgb '#008080' pt 6 ps 1 lt 1 lw 5 # --- teal

set style line 31 lc rgb '#8A2BE2' pt 6 ps 1 lt 1 lw 5 # --- blueviolet
set style line 32 lc rgb '#000E8B' pt 6 ps 1 lt 1 lw 5 # --- royalblue
set style line 33 lc rgb '#00008B' pt 6 ps 1 lt 1 lw 5 # --- darkblue
set style line 34 lc rgb '#800080' pt 6 ps 1 lt 1 lw 5 # --- purple

set style line 41 lc rgb '#2F4F4F' pt 6 ps 1 lt 1 lw 5 # --- darkslategray
set style line 42 lc rgb '#C0C0C0' pt 6 ps 1 lt 1 lw 5 # --- silver
set style line 43 lc rgb '#D2691E' pt 6 ps 1 lt 1 lw 5 # --- chocolate
set style line 44 lc rgb '#B8860B' pt 6 ps 1 lt 1 lw 5 # --- gold
set style line 45 lc rgb '#800000' pt 6 ps 1 lt 1 lw 5 # --- maroon
set style line 46 lc rgb '#8B0000' pt 6 ps 1 lt 1 lw 5 # --- darkred


set lmargin at screen 0.20
set rmargin at screen 0.93
set bmargin at screen 0.20
set tmargin at screen 0.95

set border 31 lw 4
set key font ",16"

set output "_acf.png"
set ylabel "ACF" offset 1.5, 0.0
set xlabel "Time (fs)" offset 0.0, 0.5

#set xtics 250

plot "_spectr__re__acf_0_0.txt" using ($1):2 w l ls 11 lw 8 t "ACF(0,0)",\
     "_spectr__re__acf_1_1.txt" using ($1):2 w l ls 12 lw 8 t "ACF(1,1)",\
     "_spectr__im__acf_0_1.txt" using ($1):2 w l ls 21 lw 8 t "ACF(0,1)"


set output "_IS.png"
set ylabel "IS" offset 1.5, 0.0
set xlabel "Time (fs)" offset 0.0, 0.5

plot "_spectr__re__spectrum_0_0.txt" using ($1):2 w l ls 11 lw 8 t "(0,0)",\
     "_spectr__re__spectrum_1_1.txt" using ($1):2 w l ls 12 lw 8 t "(1,1)",\
     "_spectr__im__spectrum_0_1.txt" using ($1):2 w l ls 21 lw 8 t "(0,1)"





                                  