#!/usr/bin/gnuplot

datafile = 'neldermead.log'

set multiplot layout 3,1 margins 0.12,0.95,0.12,0.95
set border lw 2
unset xlabel
set logscale y
do for [c in 'delta energy'] {
    set ylabel ''.c font ',14'
    plot datafile u 1:c w lines lw 2 lc rgb 'red'
    unset logscale y
}
set key autotitle columnhead
set xlabel 'evaluation' font ',14'
set ylabel 'parameters' font ',14'
stats datafile skip 1
plot for [i=0:(STATS_columns-5)] '' u 1:(column('p'.i)) w lines lw 2
unset multiplot
pause -1

