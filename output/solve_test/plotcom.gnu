set terminal pngcairo

stats 'output/solve_test/plots/plot1.dat' using 1:2 nooutput

Xmin = STATS_min_x
Xmax = STATS_max_x
Ymin = STATS_min_y - 0.100000 * (STATS_max_y - STATS_min_y)
Ymax = STATS_max_y + 0.100000 * (STATS_max_y - STATS_min_y)

do for [i=0:49] {
	set output sprintf('output/solve_test/animation/graphe%d.png',i) 

	set title sprintf("Animation de burgers_1d_1 len_U=%d",i) 
	set key off

	set xlabel "x" 
	set ylabel "u" 

	set xrange [Xmin:Xmax] 
	set yrange [Ymin:Ymax] 

	plot sprintf('output/solve_test/plots/plot%d.dat', i) using 1:2 w lp pt 0 
}