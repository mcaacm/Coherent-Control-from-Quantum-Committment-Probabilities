set terminal postscript enhanced colour font "Times-Roman,28"

set output "perturbations.ps"

set yrange [0.93:1.005]
set xrange [0.945:1.055]
set ytics 0.94,0.02,1.0
set xlabel "{/Symbol-Oblique w}{/Times-Italic _t/}{/Symbol-Oblique w}{/Times-Italic _0}"
set ylabel "|{/Symbol-Oblique Y}_{/Times-Italic 0}^* {/Symbol-Oblique Y}_{/Times-Italic m}|"
set key bottom center

set yrange [0.93:1.005]
set ylabel "|{/Symbol-Oblique Y}_{/Times-Italic 0}^* {/Symbol-Oblique Y}_{/Times-Italic m}|"
set y2range [0.92:1.0]
set y2tics 0.92,0.03,1.0
set y2label "{/Times-Italic P_{A|B}(}{/Symbol-Oblique r}{/Times-Italic )}"
plot "omegast_overlap.txt" u 1:(abs($4)) axes x1y1 pt 3 lc rgb "#f4a58s" ps 3 notitle, "omegast_overlap.txt" u 1:2 axes x1y2 pt 7 lc rgb "#92c5de" ps 3 title "{/Symbol-Oblique r}={/Symbol-Oblique r}_m", "omegast_overlap.txt" u 1:3 axes x1y2 pt 11 lc rgb "#f46d43" ps 3 title "{/Symbol-Oblique r}={/Symbol-Oblique r}_0"
