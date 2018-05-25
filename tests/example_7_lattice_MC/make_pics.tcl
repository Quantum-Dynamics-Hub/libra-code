for {set i 0} { $i<100 } {incr i } { 
mol load pdb step$i
set indx [ molinfo top ]

mol addrep $indx
mol modstyle 1 $indx cpk

render snapshot pic$i.png

mol delrep 1 $indx
mol delete $indx
}