# gnuplot -c arg1 arg2 plot.plt

set terminal pdfcairo size 5,3 font "Arial, 14"
filename1 = arg1
filename2 = arg2

set output "aprox_vs_interp.pdf"

set xrange [0:14]
set yrange [-250:250]
set xlabel 'x'; set ylabel 'y'
set key top left
set title 'Aproximacao Trigonometrica'
plot filename1 w lines lt -1 title "trig mmq" ,\
     filename2 w lines lt 3  title "cubic spline"