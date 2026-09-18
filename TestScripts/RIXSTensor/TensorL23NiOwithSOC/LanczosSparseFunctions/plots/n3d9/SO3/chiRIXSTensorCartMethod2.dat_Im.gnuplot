set autoscale

unset xtics

unset border
unset colorbox

set palette defined (-0.25 "blue", 0 "white", 0.25 "red")
set cbrange [-0.25:0.25]

titles = "x x x y x z y x y y y z z x z y z z"

set terminal postscript portrait enhanced color  "Verdana" 8
set output 'plots/n3d9/SO3/chiRIXSTensorCartMethod2.dat_Im.ps'

set size 0.85, 1
set multiplot layout 9, 9

set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=3;k<=23;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=25;k<=45;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=47;k<=67;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=69;k<=89;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=91;k<=111;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=113;k<=133;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=135;k<=155;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=157;k<=177;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=179;k<=199;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=201;k<=221;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=223;k<=243;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=245;k<=265;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=267;k<=287;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=289;k<=309;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=311;k<=331;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=333;k<=353;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=355;k<=375;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=377;k<=397;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=399;k<=419;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=421;k<=441;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=443;k<=463;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=465;k<=485;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=487;k<=507;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=509;k<=529;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=531;k<=551;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=553;k<=573;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=575;k<=595;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=597;k<=617;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=619;k<=639;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=641;k<=661;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=663;k<=683;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=685;k<=705;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=707;k<=727;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=729;k<=749;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=751;k<=771;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=773;k<=793;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=795;k<=815;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=817;k<=837;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=839;k<=859;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=861;k<=881;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=883;k<=903;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=905;k<=925;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=927;k<=947;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=949;k<=969;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=971;k<=991;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=993;k<=1013;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1015;k<=1035;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1037;k<=1057;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1059;k<=1079;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1081;k<=1101;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1103;k<=1123;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1125;k<=1145;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1147;k<=1167;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1169;k<=1189;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1191;k<=1211;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1213;k<=1233;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1235;k<=1255;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1257;k<=1277;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1279;k<=1299;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1301;k<=1321;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1323;k<=1343;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1345;k<=1365;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1367;k<=1387;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1389;k<=1409;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1411;k<=1431;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1433;k<=1453;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1455;k<=1475;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1477;k<=1497;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1499;k<=1519;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1521;k<=1541;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1543;k<=1563;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1565;k<=1585;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1587;k<=1607;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1609;k<=1629;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1631;k<=1651;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1653;k<=1673;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1675;k<=1695;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1697;k<=1717;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1719;k<=1739;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1741;k<=1761;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
set colorbox user origin 1.01,0.05 size 0.02,0.9
set cbtics
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1763;k<=1783;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border



unset multiplot
