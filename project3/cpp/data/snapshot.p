
set xlabel "X" font ",18"
set ylabel "Y" font ",18"

set cblabel "H" font ",18"
#set cblabel "{/Symbol f}" font ",18"

set size ratio -1

t = 8
i = t*2000  # *1000 for lvl 7, *2000 for lvl 8

a = 8       # 7 = phi (level-set func.), 8 = h (heavy side func.)

#c = "no_reinit"
#set title "w/o reinitialization | lvl 7" font ",18"

#c = "reinit"
#set title "w/ reinitialization | lvl 7" font ",18"

#c = "no_reinit_256"
#set title "w/o reinitialization | lvl 8" font ",18"

c = "reinit_256"
set title "w/ reinitialization | lvl 8" font ",18"


#plot [0:1][0:1] "".c."/".i."-snapshot-".t."" u 1:2:a w image

plot [0:1][0:1] "".c."/final-snapshot" u 1:2:a w image
