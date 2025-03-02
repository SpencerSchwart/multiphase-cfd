
set xtics font ",12"
set ytics font ",12"

set title "Re = 1000" font ",14"

set xlabel "X" font ",16"
set ylabel "v" font ",16"

set key font ",12" top right

plot [0:1]            '1000-128-y_L2' u 1:4 t "N=128" w l lw 2, \
            'y_ref_1000' u 1:2 t "Gina et al." ps 2 pt 5

# plot [0:1] '400-32-y_L2' u 1:4 t "N=32" w l lw 2, \
#         '400-64-y_L2' u 1:4 t "N=64" w l lw 2, \
#          '400-128-y_L2' u 1:4 t "N=128" w l lw 2, \
#            'y_ref_400' u 1:2 t "Gina et al." ps 2 pt 5

#plot [0:1] '100-32-y_L2' u 1:4 t "N=32" w l lw 2, \
#        '100-64-y_L2' u 1:4 t "N=64" w l lw 2, \
#            '100-128-y_L2' u 1:4 t "N=128" w l lw 2, \
#            'y_ref_100' u 1:2 t "Gina et al." ps 2 pt 5
