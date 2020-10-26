set terminal pngcairo
set output 'output/solve_burgers_1d/solve_burgers_1d_3/graphe.png'

set title "Resolution de burgers_1d_1 tmax=2.000000"
set xlabel "x"
set ylabel "u"

plot 'output/solve_burgers_1d/solve_burgers_1d_3/plot.dat' using 1:2 title "solution numerique" w lp, 'output/solve_burgers_1d/solve_burgers_1d_3/plot.dat' using 1:3 title "soluton exacte" w lp