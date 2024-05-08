set terminal postscript enhanced landscape colour font "Times-Roman,40"

Q0 = 1.0

kappa1 = -0.006836 
kappa2 = 0.006836 
omegast = 0.002279
omegasc = 0.004116
e1 = -0.00139 
e2 = 0.00139 
lambda = 0.00091153 

a(x,y) = 0.5*omegasc*y**2 + e1 + kappa1*x + 0.5*omegast*x**2
b(x,y)= e2 + kappa2*x + 0.5*omegasc*y**2 + 0.5*omegast*x**2

f(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 - sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))

 g(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 + sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))

unset colorbox
set samples 200
set isosamples 200
set cbrange [0.0:0.2]
set xrange [-6:6]
set yrange [-3:3]
set surface
set ylabel "{/Times-Italic Q_c}/au"
set xlabel "{/Times-Italic Q_t}/au"
set palette defined ( 10 '#f7f7f7', 20 '#d1e5f0', 30 '#92c5de', 40 '#4393c3', 50 '#2166ac', 60 '#053061')
set ytics -3,3,3
set xtics -6,3,6

set output "plot.ps"

set multiplot
set surface
set view map
unset contour
splot "pf_2d_20000_au.txt" u ($2/Q0):($1/Q0):3 with pm3d at b notitle

unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.005,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot

set multiplot
set surface
set view map
unset contour
splot "pf_2d_20000_au.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle

unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.005,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot



set multiplot
set surface
set view map
unset contour
splot "pf_2d_50000_au.txt" u ($2/Q0):($1/Q0):3 with pm3d at b notitle

unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.005,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot

set multiplot
set surface
set view map
unset contour
splot "pf_2d_50000_au.txt" u ($2/Q0):($1/Q0):4 with pm3d at b notitle

unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.005,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot

