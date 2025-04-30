clear 

set xlabel "x" font ",18"
set ylabel "y" font ",18"

set size ratio -1
set cblabel "level-set" font ",18"

plot [0:1][0:1] '../lvl8_myc_ei/0-snapshot-0' u 1:2:8 w image not
