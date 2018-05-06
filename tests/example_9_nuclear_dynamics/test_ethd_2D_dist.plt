#
set terminal pngcairo font "arial,20" size 800, 600 enhanced rounded truecolor

#set lmargin at screen 0.17
#set rmargin at screen 0.95
#set bmargin at screen 0.15
#set tmargin at screen 0.95

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


do for [i=0:1]    {
    set output "_initial dist_".i.".png"
    set xlabel "DOF #".i." axis" offset 0.0, 0.5
    set ylabel "Probabilitiy" offset 1.7, 0.0
    set xrange [-4.0:4.0]
    plot  "_q_distrib_dof_".i.".txt" w l ls 32 lw 5 t "q - DOF #".i."",\
          "_p_distrib_dof_".i.".txt" w l ls 12 lw 5 t "p - DOF #".i."",\
}

set nokey
set palette model CMY rgbformulae 7,5,15
# 2D pictures, colored by z-value
set contour
set pm3d map
set pm3d interpolate 0,0  # first - y axis, second - x axis
set pm3d explicit at b
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model CMY 
set cbrange [0:0.1]  
set colorbox user origin 0.05,0.25
set ylabel "DOF 2 (a.u)" offset -1.75, 5.5
set view 0, 0, 1, 1
set xrange[-4:4]
set yrange[-5:3]
set zrange[0.000005:]
unset ztics
### Enter bsed on your specific simulation ###
Nsteps = 20
dt = 1

do for [i=0:149]{
    set output "_2D_dist_output/Disttribution_".i.".png"
   
    set xlabel "DOF 1 (a.u), t = ".i*(dt*Nsteps)." (a.u) "  offset 0.0, -0.5
    splot "_pes_scan.txt" u ($1-70)*0.1:($2-50)*0.1:($3) w pm3d,\
          "_2D_dist/_2D_distrib_".i."_.txt" u ($1-70)*0.1:($2-50)*0.1:($3) w p nocontour pt 7 ps 1 lc "black",\
}

