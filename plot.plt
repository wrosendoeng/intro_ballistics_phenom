# gnuplot -c arg1 arg2 plot.plt
set terminal pngcairo size 960,640 font "Arial, 14" enhanced
set termoption enhanced
set encoding iso_8859_1

set output "picture_blasius.png"

set xrange [0:10]
set yrange [0:2]
set xlabel '{/Symbol h}'; set ylabel 'y'
set key top right
set title "Problema de Blasius ({/Symbol Dh} = 0.01)" font "Arial,16"
plot 'results.txt' u 1:2 w lines dt 2 lc rgb "red" lw 1.5 title 'f({/Symbol h})' ,\
     'results.txt' u 1:3 w lines dt 6 lc rgb "black" lw 1.5 title "f'({/Symbol h})",\
     'results.txt' u 1:4 w lines dt 10 lc rgb "blue" title "f''({/Symbol h})"