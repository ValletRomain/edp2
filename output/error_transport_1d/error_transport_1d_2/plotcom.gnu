set terminal pngcairo

# Graphic of error
set output 'output/error_transport_1d/error_transport_1d_2/error.png'

set title "Erreur en norme L1"
set xlabel "N"
set ylabel "error"

set logscale x 10
plot 'output/error_transport_1d/error_transport_1d_2/error.dat' using 1:2 title "error"

# Graphic of time
set output 'output/error_transport_1d/error_transport_1d_2/time.png'

set title "Duree"
set xlabel "N"
set ylabel "time (s)"

set logscale x 10
set yrange [-1:10]

plot 'output/error_transport_1d/error_transport_1d_2/time.dat' using 1:2 title "time"