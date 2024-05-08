# This script is used to plot the population plots found in
# figure 3 (c and d) of the paper.

set terminal postscript enhanced colour font "Times-Roman,32"
set output "populations.ps"

set xlabel "{/Times-Italic t}/au"
set ylabel "{/Symbol-Oblique r}({/Times-Italic t})_{33}"
set key top right

set xrange [0:50000]
set xtics 0,15000,45000

set yrange [0:1.0]
set key bottom right
set ylabel "{/Symbol-Oblique r}({/Times-Italic t})_{{/Times-Italic ii}}"
plot "nonsecular_relax_optimal_to_es1.txt" u 1:2 with lines lw 8 lc rgb "#f4a58s" title "{/Times-Italic i}=1", "nonsecular_relax_optimal_to_es1.txt" u 1:3 with lines lw 8 lc rgb "#92c5de" title "{/Times-Italic i}=2", "nonsecular_relax_optimal_to_es1.txt" u 1:4 with lines lw 8 lc rgb "#b2182b" title "{/Times-Italic i}=3"
unset key
plot "nonsecular_relax_optimal_to_es3.txt" u 1:2 with lines lw 8 lc rgb "#f4a58s" title "{/Times-Italic i}=1", "nonsecular_relax_optimal_to_es3.txt" u 1:3 with lines lw 8 lc rgb "#92c5de" title "{/Times-Italic i}=2", "nonsecular_relax_optimal_to_es3.txt" u 1:4 with lines lw 8 lc rgb "#b2182b" title "{/Times-Italic i}=3"
