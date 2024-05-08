# Plots the population evolution plots from several different
# starting points and with the three different propagation 
# techniques.

set terminal postscript enhanced colour font "Times-Roman,32"
set output "figure.ps"

set xlabel "{/Times-Italic t}/au"
set ylabel "{/Symbol-Oblique r}({/Times-Italic t})_{33}"
set key top right

set xrange [0:5000]
set xtics 0,1500,4500

set yrange [0:1.0]

set ylabel "{/Symbol-Oblique r}({/Times-Italic t})_{11}"
plot "wf_-0.87_-0.17_-0.31.txt" u 1:2 with lines lw 8 lc rgb "#f4a58s" title "{/Symbol-Oblique r}(0)={/Symbol-Oblique r}_1", "wf_0.87_0.383_-0.31.txt" u 1:2 with lines lw 8 lc rgb "#92c5de" title "{/Symbol-Oblique r}(0)={/Symbol-Oblique r}_2", "wf_0.87_0.17_0.46.txt" u 1:2 with lines lw 8 lc rgb "#b2182b" title "{/Symbol-Oblique r}(0)={/Symbol-Oblique r}_3" #, "above_nscp_on_scp.txt" u 1:6 with lines lw 8 lc rgb "#f46d43" title "{/Symbol-Oblique r}(0)={/Symbol-Oblique r}_1 Coherence 3,4"



set ylabel "{/Symbol-Oblique r}({/Times-Italic t})_{11}" #, {/Symbol-Oblique r}(0)={/Symbol-Oblique r}_2"
plot "wf_0.87_0.383_-0.31.txt" u 1:2 with lines lw 8 lc rgb "#f4a58s" title "Partial Secular", "wf_0.87_0.383_-0.31_coherence_killed.txt" u 1:2 with lines lw 8 lc rgb "#92c5de" title "Zeroed Coherence", "wf_0.87_0.383_-0.31_secular.txt" u 1:2 with lines lw 8 lc rgb "#b2182b" title "Secular"
