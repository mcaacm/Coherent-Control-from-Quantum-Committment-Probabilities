set terminal postscript enhanced colour font "Times-Roman,40"
set output "plot.ps"

Q0 = 1.0

kappa1 = -0.006836 #-0.0050 #-0.003858675*2.0
kappa2 = 0.006836 #0.0050 #0.0054756457*1.6
omegast = 0.002279*0.95 #0.003000 #0.0043364174
omegasc = 0.004116 #0.006400 #0.0027194482
e1 = -0.00139 #-0.00139 #0.0 #0.144792242
e2 = 0.00139 #0.00139 #0.00060 #0.177866612*0.87
lambda = 0.00091153 #0.090*0.096283166 #*0.83

a(x,y) = 0.5*omegasc*y**2 + e1 + kappa1*x + 0.5*omegast*x**2
b(x,y)= e2 + kappa2*x + 0.5*omegasc*y**2 + 0.5*omegast*x**2

f(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 - sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))

 g(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 + sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))


set samples 300
set isosamples 300
set ztics -0.02,0.02,0.02
unset colorbox
set xrange [-6:6]
set yrange [-3:3]
set zrange [-0.02:0.02]
#set cbrange [-0.02:0.02]

set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
set ylabel "{/Times-Italic q_c}/au"
set xlabel "{/Times-Italic q_t}/au"
set xyplane -0.02
set ytics -0.02,0.02,0.02

set terminal postscript enhanced landscape colour font "Times-Roman,40"
reset
unset colorbox

set samples 200
set isosamples 200
set cbrange [-0.4:0.4]
set xrange [-6:6]
set yrange [-3:3]
set surface
set ylabel "{/Times-Italic q_c}/au"
set xlabel "{/Times-Italic q_t}/au"

set multiplot
set ytics -3,3,3
set xtics -6,3,6
set surface
set view map
unset contour
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "wf_2d_0.txt" u ($2/Q0):($1/Q0):5 with pm3d at b notitle

unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.005,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot


set ylabel "{/Times-Italic q_c}/au"
set xlabel "{/Times-Italic q_t}/au"

set multiplot
set surface
set view map
unset contour
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
splot "wf_2d_0.txt" u ($2/Q0):($1/Q0):($7*20) with pm3d at b notitle

unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.005,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot

