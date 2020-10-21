set terminal pngcairo
set output 'output/solve_transport_1d/solve_transport_1d_2/graphe.png'

set title "Resolution de l'equation de transport"
set xlabel "x"
set ylabel "u"

set yrange [0:1.2]

plot 'output/solve_transport_1d/solve_transport_1d_2/plot.dat' using 1:2 title "solution numerique", 'output/solve_transport_1d/solve_transport_1d_2/plot.dat' using 1:3 title "soluton exacte"