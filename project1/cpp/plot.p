
set xlabel "x" font ",14"
set ylabel "y" font ",14"

set xrange[0:1]
set yrange[0:1]

set size square

set cbrange [0:0.5]

folder = "upwind-mu"
timestep = "780000"

set title "Tracer Field @i=".timestep." | ".folder."" font ",14"

plot for [f in system("ls data/".folder."/".timestep."-tracer-*")] f using 1:2:5 w image



