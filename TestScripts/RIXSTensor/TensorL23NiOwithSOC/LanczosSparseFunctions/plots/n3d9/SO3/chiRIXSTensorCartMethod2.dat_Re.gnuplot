set autoscale

unset xtics

unset border
unset colorbox

set palette defined (-0.25 "blue", 0 "white", 0.25 "red")
set cbrange [-0.25:0.25]

titles = "x x x y x z y x y y y z z x z y z z"

set terminal postscript portrait enhanced color  "Verdana" 8
set output 'plots/n3d9/SO3/chiRIXSTensorCartMethod2.dat_Re.ps'

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
plot "< awk '(NR>5){for(k=2;k<=22;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=24;k<=44;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=46;k<=66;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=68;k<=88;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=90;k<=110;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=112;k<=132;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=134;k<=154;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=156;k<=176;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*1 - 1).word(titles, 2*1) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=178;k<=198;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=200;k<=220;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=222;k<=242;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=244;k<=264;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=266;k<=286;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=288;k<=308;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=310;k<=330;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=332;k<=352;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=354;k<=374;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*2 - 1).word(titles, 2*2) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=376;k<=396;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=398;k<=418;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=420;k<=440;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=442;k<=462;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=464;k<=484;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=486;k<=506;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=508;k<=528;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=530;k<=550;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=552;k<=572;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*3 - 1).word(titles, 2*3) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=574;k<=594;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=596;k<=616;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=618;k<=638;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=640;k<=660;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=662;k<=682;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=684;k<=704;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=706;k<=726;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=728;k<=748;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=750;k<=770;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*4 - 1).word(titles, 2*4) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=772;k<=792;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=794;k<=814;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=816;k<=836;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=838;k<=858;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=860;k<=880;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=882;k<=902;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=904;k<=924;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=926;k<=946;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=948;k<=968;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*5 - 1).word(titles, 2*5) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=970;k<=990;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=992;k<=1012;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1014;k<=1034;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1036;k<=1056;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1058;k<=1078;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1080;k<=1100;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1102;k<=1122;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1124;k<=1144;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1146;k<=1166;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*6 - 1).word(titles, 2*6) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1168;k<=1188;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1190;k<=1210;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1212;k<=1232;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1234;k<=1254;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1256;k<=1276;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1278;k<=1298;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1300;k<=1320;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1322;k<=1342;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1344;k<=1364;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*7 - 1).word(titles, 2*7) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1366;k<=1386;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1388;k<=1408;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1410;k<=1430;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1432;k<=1452;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1454;k<=1474;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1476;k<=1496;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1498;k<=1518;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1520;k<=1540;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1542;k<=1562;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*9 - 1).word(titles, 2*9)." | ".word(titles, 2*8 - 1).word(titles, 2*8) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1564;k<=1584;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y "%g"
unset colorbox
#set xlabel word(titles, 2*1 - 1).word(titles, 2*1)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1586;k<=1606;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*2 - 1).word(titles, 2*2)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1608;k<=1628;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*3 - 1).word(titles, 2*3)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1630;k<=1650;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*4 - 1).word(titles, 2*4)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1652;k<=1672;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*5 - 1).word(titles, 2*5)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1674;k<=1694;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*6 - 1).word(titles, 2*6)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1696;k<=1716;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*7 - 1).word(titles, 2*7)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1718;k<=1738;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border

set border 15
set ytics
set format y ""
unset colorbox
#set xlabel word(titles, 2*8 - 1).word(titles, 2*8)." | ".word(titles, 2*9 - 1).word(titles, 2*9) font "Verdana,8"
set xlabel offset 0, 10
set xrange [0:1]
set yrange [0:4.02]
plot "< awk '(NR>5){for(k=1740;k<=1760;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
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
plot "< awk '(NR>5){for(k=1762;k<=1782;k+=2){printf \"%s \", $k}{printf \"%s\", RS}}' data/n3d9/SO3/chiRIXSTensorCartMethod2.dat" matrix using ($1*0.1+0):($2*0.02+0):3 with image notitle
unset border



unset multiplot
