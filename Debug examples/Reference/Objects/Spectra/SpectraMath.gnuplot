set autoscale                          # scale axes automatically
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set style line  1 lt 1 lw 1 lc rgb "#000000"
set style line  2 lt 1 lw 1 lc rgb "#FF0000"
set style line  3 lt 1 lw 1 lc rgb "#00FF00"
set style line  4 lt 1 lw 1 lc rgb "#0000FF"
set style line  5 lt 1 lw 3 lc rgb "#808080"
set style line  6 lt 1 lw 3 lc rgb "#00FFFF"
set style line  7 lt 1 lw 3 lc rgb "#FF00FF"
set style line  8 lt 1 lw 3 lc rgb "#FFFF00"

set xlabel "Energy (arb. units)" font "Times-Roman,12"
set ylabel "Intensity (arb. units)" font "Times-Roman,12"

set size 1.0, 0.3
set terminal postscript portrait enhanced color  "Times-Roman" 8
     
set output "SpectraMath.ps"

plot "SpectraMath.dat"  using 1:(-$3  ) title 'G ' with lines ls  1,\
     "SpectraMath1.dat" using 1:(-$3  ) title 'G1' with lines ls  2,\
     "SpectraMath2.dat" using 1:(-$3  ) title 'G2' with lines ls  3
