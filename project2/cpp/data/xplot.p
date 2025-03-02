
set xtics font ",12"
set ytics font ",12"

set title "Re = 1000" font ",14"

set xlabel "Y" font ",16"
set ylabel "u" font ",16"

set key font ",12" top center

plot [0:1]            '1000-128-x_L2' u 2:3 t "N=128" w l lw 2, \
            'x_ref_1000' u 1:2 t "Gina et al." ps 2 pt 5

# plot [0:1] '400-32-x_L2' u 2:3 t "N=32" w l lw 2, \
#         '400-64-x_L2' u 2:3 t "N=64" w l lw 2, \
#            '400-128-x_L2' u 2:3 t "N=128" w l lw 2, \
#            'x_ref_400' u 1:2 t "Gina et al." ps 2 pt 5

#plot [0:1] '100-32-x_L2' u 2:3 t "N=32" w l lw 2, \
#        '100-64-x_L2' u 2:3 t "N=64" w l lw 2, \
#            '100-128-x_L2' u 2:3 t "N=128" w l lw 2, \
#            'x_ref_100' u 1:2 t "Gina et al." ps 2 pt 5
