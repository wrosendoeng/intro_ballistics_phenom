#!/bin/bash

./executecode.sh

frame ( )
  {
    echo "set term pdfcairo size 5,3 font 'Arial,14'"
    echo "set output 'rk4_questao04plot.pdf'"
    echo "set xlabel 'time (s)'"
    echo -e "set format x \"%g\" "
    echo "set xrange [100:200]"
    echo "set ylabel 'velocity (m/s)'"
    echo "set title 'Particles movement - RK4 method'"
    echo -e "set format y \"%g\" "
    echo "plot 'resultados_novos.txt' using 1:5 w lines lt 1 title 'corpo1', 'resultados_novos.txt' using 1:6 w lines lt 2 title 'corpo2', 'resultados_novos.txt' using 1:7 w lines lt 8 title 'corpo3' "
    return
  }

frame | gnuplot
