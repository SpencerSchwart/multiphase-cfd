clear

load '~/basilisk/reconstructed-ibm/figures/linestyles.p'

set xlabel "X" font ",18"
set ylabel "Y" font ",18"

set xtics 0.1 font ",11"
set ytics 0.1 font ",11"

set mxtics 2
set mytics 2

set size ratio -1

set key top left font ",15"
unset key

unset ylabel

t=9
i=t*2000

set xrange [0.4:1]
set yrange [0.4:1]

plot '../lvl8_myc_ei/'.i.'-interface-'.t.'' u 1:2 w l ls 2 t "EI w/CCM", \
     '../lvl8_yfd_ei/'.i.'-interface-'.t.'' u 1:2 w l ls 3 t "EI w/YFD", \
     '../lvl8_myc_bas/'.i.'-interface-'.t.'' u 1:2 w l ls 1 dt 3 t "BAS w/CCM", \
     '../lvl8_yfd_bas/'.i.'-interface-'.t.'' u 1:2 w l ls 4 dt 3 t "BAS w/YFD"

plot '../lvl8_myc_ei/final-interface' u 1:2 w l ls 2 t "EI w/CCM", \
     '../lvl8_yfd_ei/final-interface' u 1:2 w l ls 3 t "EI w/YFD", \
     '../lvl8_myc_bas/final-interface' u 1:2 w l ls 1 dt 3 t "BAS w/CCM", \
     '../lvl8_yfd_bas/final-interface' u 1:2 w l ls 4 dt 3 t "BAS w/YFD"
