clear

load '~/basilisk/reconstructed-ibm/figures/linestyles.p'

set xlabel "tU_0/L" font ",18"
set ylabel "Volume" font ",18"

set xtics 2 font ",11"
set ytics 0.002 font ",11"

set mxtics 2
set mytics 2

set key center right font ",15"

plot '../lvl8_ccm_ei/log' u 2:8 w l ls 2 t "EI w/CCM", \
     '../lvl8_yfd_ei/log' u 2:8 w l ls 3 t "EI w/YFD", \
     '../lvl8_myc_bas/log' u 2:8 w l ls 1 dt 3 t "BAS w/CCM", \
     '../lvl8_yfd_bas/log' u 2:8 w l ls 4 dt 3 t "BAS w/YFD"
