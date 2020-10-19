set terminal pngcairo
set output 'output/solve_test/graphe.png'

set title "Resolution de l'equation de transport"
set xlabel "x"
set ylabel "u"

set yrange [0:1.2]

plot 'output/solve_test/plot.dat' using 1:2 title "solution exacte", 'output/solve_test/plot.dat' using 1:3 title "soluton numerique"