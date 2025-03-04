clear

set xlabel "X" font ",14"
set ylabel "Y" font ",14"

set size square

set title "Re = 1000 | t = 8.25" font ",14"

# set cbrange [0:1]

plot [0:1][0:1] 're1000_central/5500-snapshot-8.25' u 1:2:7 w image
# plot [0:1][0:1] 're100_central/final-snapshot' u 1:2:7 w image
