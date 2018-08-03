#!/bin/gnuplot
set term wxt enhanced dashed "arial,16"
set datafile separator ","

set title "Edge Standard Deviation vs Performance"
set key below
set xrange [0.5:5.5]
set yrange [0:70]
set xlabel 'Standard Deviation of Contigent Edges (sec)'
set ylabel 'Empirical Success Rate (%)'


v0='dyn'
v1='early'
v2='stat'
v3='dc'
v4='n_stat'

plot  v0 with lines lt 1 lc 1 title 'DREA', \
    '' with errorbars notitle lc 1 ps 2, \
    v1 with lines lt 2 lc 2 title 'Early Execution', \
    '' with errorbars notitle lc 2 ps 2, \
    v2 with lines lt 3 lc 3 title 'SREA', \
    '' with errorbars notitle lc 3 ps 2, \
    v3 with lines lt 4 lc -1 title 'Dyn. Contr.', \
    '' with errorbars notitle lc -1 ps 2, \
    v4 with lines lt 5 lc 4 title 'N. SREA', \
    '' with errorbars notitle lc 4 ps 2,

set term pngcairo dashed size 1024,768
set output "plot.png"
replot
set term x11
