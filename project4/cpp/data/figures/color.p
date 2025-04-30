clear

load '~/basilisk/reconstructed-ibm/figures/linestyles.p'

set xlabel "X" font ",18"
set ylabel "Y" font ",18"

set xtics 0.25 font ",11"
set ytics 0.25 font ",11"

set mxtics 2
set mytics 2

set size ratio -1

set key top left font ",15"
unset key
unset ylabel

t=8
i=t*2000

set xrange [0:1]
set yrange [0:1]

plot '../lvl8_yfd_ei/'.i.'-snapshot-'.t.'' u 1:2:7 w image not

plot '../lvl8_yfd_ei/final-snapshot' u 1:2:7 w image not

