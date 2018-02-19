#
set terminal pngcairo font "arial,24" size 800, 600 enhanced rounded truecolor

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



set xlabel "time, fs" offset 0.0, 0.5
set ylabel "Energy, eV" offset 1.5, 0.0
#set xtics 1000.0
#set xrange [0.0:4.5]
#set yrange [-40:40]
set key spacing 1.0 font ",24"
set key top horizontal



set output "energy.png"
plot "_energy-adi.txt" using 2:4 w l ls 11 lw 5 t "SE energy"


set output "_pop_adi_se.png"
plot "_pop_adi_se.txt" using 2:5 w l ls 11 lw 5 t "P(0)",\
     "_pop_adi_se.txt" using 2:8 w l ls 21 lw 5 t "P(1)",\
     "_pop_adi_se.txt" using 2:11 w l ls 31 lw 5 t "P(2)"



set output "spectrum.png"
plot "spectrum.txt" using ($1*27.3):2 w l ls 11 lw 5 t "SE energy"



set nokey
#set palette rgbformulae 22,13,-31
set palette rgbformulae 7,5,15

# 3D pictures, colored by z-value 
#set auto
#set parametric
#set pm3d implicit at s

# 2D pictures, colored by z-value
#set contour
set pm3d map
set pm3d interpolate 0,0  # first - y axis, second - x axis
set pm3d explicit at b

set output "NAC.png"
splot "_Hvib_im.dat" using ($1):($2):($3*1000)


set output "TMO.png"
splot "Tmo.dat" using ($1):($2):($3*1000)
