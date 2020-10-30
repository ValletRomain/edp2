set terminal pngcairo
set output 'output/solve_burgers_1d/solve_burgers_1d_3/graphe.png'

set title "Resolution de burgers_1d_1 tmax=2.000000"
set xlabel "x"
set ylabel "u"

stats 'output/solve_burgers_1d/solve_burgers_1d_3/plot.dat' using 1:2 nooutput
set xrange [STATS_min_x:STATS_max_x]
set yrange [STATS_min_y - 0.100000 * (STATS_max_y-STATS_min_y): STATS_max_y + 0.100000 * (STATS_max_y-STATS_min_y)]

plot 'output/solve_burgers_1d/solve_burgers_1d_3/plot.dat' using 1:2 title "solution numerique" w lp pt 0, 'output/solve_burgers_1d/solve_burgers_1d_3/plot.dat' using 1:3 title "soluton exacte" w lp pt 0