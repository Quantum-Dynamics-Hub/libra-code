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





set xlabel "t, fs" offset 0.0, 0.5
set ylabel "Energy, a.u." offset 1.5, 0.0 
set xrange [0:5000]
set xtics 500
set yrange [-0.02:0.16]
set ytics 0.04
set output "_tsh_en.png"
plot "_tsh_extern.txt" using ($2/41.0):8   w l  ls 11  lw 5  t "Ekin",\
     "_tsh_extern.txt" using ($2/41.0):10  w l  ls 21  lw 5  t "Epot",\
     "_tsh_extern.txt" using ($2/41.0):12  w l  ls 32  lw 5  t "Etot",\
     "_tsh_extern.txt" using ($2/41.0):14  w l  ls 41  lw 5  t "Eext"


set xlabel "q, a.u." offset 0.0, 0.5
set ylabel "p, a.u." offset 1.5, 0.0 
set xrange [-40:30]
set xtics 5
set yrange [-10:10]
set ytics 5

set output "_tsh_phase.png"
plot "_tsh_extern.txt" using 4:6  w l  ls 11  lw 5  t "Phase space"


set xlabel "t, fs" offset 0.0, 0.5
set ylabel "Density matrix" offset 1.5, 0.0 
set xrange [0:41000]
set xtics 500
set yrange [-0.6:1.1]
set ytics 0.4

set output "_tsh_pop.png"
plot "_tsh_extern.txt" using ($2/41.0):16  w l  ls 11  lw 5  t "|c0|^2",\
     "_tsh_extern.txt" using ($2/41.0):18  w l  ls 21  lw 5  t "|c1|^2",\
     "_tsh_extern.txt" using ($2/41.0):20  w l  ls 32  lw 5  t "Re(c01)",\
     "_tsh_extern.txt" using ($2/41.0):24  w l  ls 11  lw 3  t "SH, 0",\
     "_tsh_extern.txt" using ($2/41.0):26  w l  ls 21  lw 3  t "SH, 1"


set xlabel "t, fs" offset 0.0, 0.5
set ylabel "Energy, a.u." offset 1.5, 0.0 
set xrange [0:41000]
set xtics 500
set yrange [-0.002:0.002]
set ytics 0.0005

set output "_tsh_Ham1.png"
plot "_tsh_extern.txt" using ($2/41.0):28  w l  ls 11  lw 5  t "E_0",\
     "_tsh_extern.txt" using ($2/41.0):30  w l  ls 21  lw 5  t "E_1"


set xlabel "t, fs" offset 0.0, 0.5
set ylabel "NAC, a.u." offset 1.5, 0.0 
set xrange [0:41000]
set xtics 500
set yrange [-0.001:0.001]
set ytics 0.0005

set output "_tsh_Ham2.png"
plot "_tsh_extern.txt" using ($2/41.0):32  w l  ls 11  lw 5  t "D_{01}"




set xlabel "t" offset 0.0, 0.5
set ylabel "istate" offset 1.5, 0.0 
set xrange [0:5000]
set xtics 500
set yrange [-0.1:1.1]
set ytics 0.2

set output "_tsh_istate.png"
plot "_tsh_extern.txt" using ($2/41.0):22   w p  ls 11  lw 3  t "istate"





##======= ACF ==========
set xlabel "t, fs" offset 0.0, 0.5
set ylabel "ACF" offset 1.5, 0.0 
set xrange [0:1000]
set xtics 250
set yrange [-1.1:1.1]
set ytics 0.2

set output "ACF.png"
plot "H_vib_im_D01__acf0_1.txt" using ($1):2   w p  ls 11  lw 3  t "ACF im 0,1",\
     "H_vib_re_E0__acf0_0.txt" using ($1):2   w p  ls 21  lw 3  t "ACF re 0,0",\
     "H_vib_re_E1__acf1_1.txt" using ($1):2   w p  ls 31  lw 3  t "ACF re 1,1"


##======= Spectrum ==========
set xlabel "Frequency, cm^{-1}" offset 0.0, 0.5
set ylabel "J" offset 1.5, 0.0 
set xrange [0:3000]
set xtics 500
set yrange [-1000: 50000000]
set ytics 5000000

set output "Spectrum.png"
plot "H_vib_im_D01__spectrum0_1.txt" using ($1):($2*$2)   w l  ls 11  lw 3  t "J im 0,1",\
     "H_vib_re_E0__spectrum0_0.txt" using ($1):($2*$2)   w l  ls 21  lw 3  t "J re 0,0",\
     "H_vib_re_E1__spectrum1_1.txt" using ($1):($2*$2)   w l  ls 31  lw 3  t "J re 1,1"




