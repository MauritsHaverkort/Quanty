set autoscale
set xtic auto
set ytic auto

set xlabel "zeta in units of 10Dq" font "Times,12"
set ylabel "Energy (eV)" font "Times,12"

set out 'EnergyLevelDiagram.ps'
set size 1.0, 0.625
set terminal postscript portrait enhanced color  "Times" 12

# set this to the range of your variable which you want to color-encode
# or leave it out
set cbrange [0.5:1.5]

# define the palette to your liking
set palette defined ( 0.5 "#FF0000", 0.75 "#00FF00", 1.25 "#00FF00", 1.5 "#0000FF" )

# in this example, column 3 is mapped to the colors of the palette
plot for [i=0:9] "EnergyLevelDiagram" using 1:2*i+2:2*i+3 with lines lc palette z notitle
