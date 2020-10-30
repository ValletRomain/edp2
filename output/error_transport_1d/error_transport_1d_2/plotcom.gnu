set terminal pngcairo

# Graphic of error
set output 'output/error_transport_1d/error_transport_1d_2/error.png'

set title "Erreur du schéma de Godunov pour transport_1d_1 en norm_L2"
set xlabel "N"
set ylabel "error"

set logscale x 10
stats 'output/error_transport_1d/error_transport_1d_2/plot.dat' using 1:2 nooutput
set xrange [STATS_min_x:STATS_max_x]
set yrange [0: STATS_max_y + 0.100000 * (STATS_max_y-STATS_min_y)]

plot 'output/error_transport_1d/error_transport_1d_2/plot.dat' using 1:2 title "error" w lp

reset

# Graphic of time
set output 'output/error_transport_1d/error_transport_1d_2/time.png'

set title "Duree"
set xlabel "N"
set ylabel "time (s)"

set autoscale y
stats 'output/error_transport_1d/error_transport_1d_2/plot.dat' using 1:3 nooutput
set xrange [STATS_min_x:STATS_max_x]
plot 'output/error_transport_1d/error_transport_1d_2/plot.dat' using 1:3 title "time" w lp