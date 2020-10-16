set terminal pngcairo
set output 'output/benchmarks2/graphe.png'

set title "Resolution de l'equation de transport"
set xlabel "x"
set ylabel "u"

set yrange [0:1.2]

plot 'output/benchmarks2/plot.dat' using 1:2 title "solution exacte", \
'output/benchmarks2/plot.dat' using 1:3 title "soluton numerique"