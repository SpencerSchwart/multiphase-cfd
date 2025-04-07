
set xlabel "Time" font ",18"
set ylabel "Volume" font ",18"

set key font ",16"

plot 'no_reinit/log' u 2:8 w l lw 2 t "w/o reinit | lvl 7", \
     'reinit/log' u 2:8 w l lw 2 dt 2  t "w/reinit | lvl 7", \
     'no_reinit_256/log' u 2:8 w l lw 2 t "w/reinit | lvl 8", \
     'reinit_256/log' u 2:8 w l lw 2 dt 2 t "w/reinit | lvl 8"
