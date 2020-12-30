set terminal pngcairo
set output 'graphe.png'

set title "Riemann Saint Venant"
set xlabel "x"

#plot 'output/solve_test/plot.dat' using 1:2 title "godunov" w lp pt 0 lc rgb "blue"
plot 'riem.dat' using 1:2 title "u" w lp pt 0 lc rgb "blue"
plot 'riem.dat' using 1:3 title "h" w lp pt 0 lc rgb "red"