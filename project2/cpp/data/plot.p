
set xlabel "x" font ",14"
set ylabel "v" font ",14"

set xrange[0:1]
set yrange[0:1]

#set size square

#set cbrange [0:1]

set title "RE = 100" font ",14"

set key top right font ",12"

plot 'ref/y_ref_100' u 1:2 lw 2 lc rgb "black" t "Ghia et al. 1982", \
        're100_central/final-ycut' u 1:4 w l lc rgb "dark-green" t "cpp central", \
            're100_upwind/final-ycut' u 1:4 w l lc rgb "dark-green" dt 2 t "cpp upwind", \
                '../../basilisk/data/re100_central/final-ycut' u 1:4 w l lc rgb "red" t "basilisk central", \
                '../../basilisk/data/re100_upwind/final-ycut' u 1:4 w l lc rgb "red" dt 2 t "basilisk upwind"

