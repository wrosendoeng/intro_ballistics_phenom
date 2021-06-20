#!/bin/bash

./executecode.sh

frame ( )
  {
    echo "set terminal pngcairo size 960,640 font 'Arial, 14' enhanced"
    echo "set termoption enhanced"
    echo "set encoding iso_8859_1"
    echo "set output 'picture_blasius.png'"
    echo "set xrange [0:10] "
    echo "set yrange [0:2]"
    echo "set ylabel 'y'"
    echo "set title 'Problema de Blasius ({/Symbol Dh} = 0.01)'"
    echo "plot 'results.txt' u 1:2 w lines dt 2 lc rgb 'red' lw 1.5 title 'f({/Symbol h})',\
    'results.txt' u 1:3 w lines dt 6 lc rgb 'black' lw 1.5 title "f\'""({/Symbol h}),\"
    'results.txt' u 1:4 w lines dt 10 lc rgb 'blue' title "f\''""({/Symbol h})"
    return
  }

frame | gnuplot
