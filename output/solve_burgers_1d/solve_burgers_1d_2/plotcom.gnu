set terminal pngcairo
set output 'output/solve_burgers_1d/solve_burgers_1d_2/graphe.png'

set title "Resolution de burgers_1d_1 tmax=1.000000"
set xlabel "x"
set ylabel "u"

set yrange [0:1.2]

plot 'output/solve_burgers_1d/solve_burgers_1d_2/plot.dat' using 1:2 title "solution numerique" w lp, 'output/solve_burgers_1d/solve_burgers_1d_2/plot.dat' using 1:3 title "soluton exacte" w lp