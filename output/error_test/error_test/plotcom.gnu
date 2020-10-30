set terminal pngcairo

# Graphic of error
set output 'output/error_test/error_test/error.png'

set title "Erreur en norme L1"
set xlabel "N"
set ylabel "error"

set logscale x 10
plot 'output/error_test/error_test/plot.dat' using 1:2 title "error" w lp

# Graphic of time
set output 'output/error_test/error_test/time.png'

set title "Duree"
set xlabel "N"
set ylabel "time (s)"

set logscale x 10
set yrange [-1:10]

plot 'output/error_test/error_test/plot.dat' using 1:3 title "time" w lp