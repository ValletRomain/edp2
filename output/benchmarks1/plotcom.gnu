set terminal pngcairo
set output 'output/benchmarks1/graphe.png'

set title "Resolution de l'equation de transport"
set xlabel "x"
set ylabel "u"

set yrange [0:1.2]

plot 'output/benchmarks1/plot.dat' using 1:2 title "solution exacte", \
'output/benchmarks1/plot.dat' using 1:3 title "soluton numerique"