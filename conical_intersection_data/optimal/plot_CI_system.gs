# This file is used to plot wavefunctions on top of contour maps.
# These are to be the wavefunctions and are read out
# of files 'wf_2d_1.txt' and 'wf_2d_0.txt.'
# The two electronic states are plotted separately. The real
# part of the wave function is displayed. 
# It reads in the eigenvector energies from several different
# files in order to give them different colors to plot on the
#parabolas.

set terminal postscript enhanced colour font "Times-Roman,30"
set output "parabolas.ps"

# Set up the conical intersection
Q0 = 1.0  # Normalization factor no longer used
kappa1 = -0.006836 
kappa2 = 0.006836 
omegast = 0.002279 
omegasc = 0.004116 
e1 = -0.00139 
e2 = 0.00139 
lambda = 0.00091153 

set xrange [-5.0:5.0]
set yrange [-0.35:0.05] #[3.6:4.66]
set samples 10000
set ylabel "{/Times-Italic E}/eV"# font "Times-Roman,18"
set xlabel "{/Times-Italic q_{t}}/au" # at {/Times-Italic q_{c}}=0" #font "Times-Roman,18"
# Diabatic surfaces
a(x,y) = 0.5*omegasc*y**2 + e1 + kappa1*x + 0.5*omegast*x**2
b(x,y)= e2 + kappa2*x + 0.5*omegasc*y**2 + 0.5*omegast*x**2

# Adiabatic surfaces
f(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 - sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))
 g(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 + sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))


# Plot the diabatic surfaces with the eigenvector energies.
set ytics -0.3,0.1,0
set size 0.7,1.0
plot 27.211*a(x,0) lw 6 lc "black" notitle, 27.211*b(x,0) lw 6 lc "black" notitle, "plot_ev_1.txt" u 2:($3*27.211) with lines lw 6 lc rgb "#3288bd" notitle, "plot_ev_2.txt" u 2:($3*27.211) with lines lw 6 lc rgb "#1a9850" notitle, "plot_ev_3.txt" u 2:($3*27.211) with lines lw 6 lc rgb "#a6d96a" notitle, "plot_ev_4.txt" u 2:($3*27.211) with lines lw 6 lc rgb "#fee08b" notitle, "plot_ev_5.txt" u 2:($3*27.211) with lines lw 6 lc rgb "#fc8d59" notitle, "plot_ev_6.txt" u 2:($3*27.211) with lines lw 6 lc rgb "#d53e4f" notitle, "plot_ev_13.txt" u 2:($3*27.211) with lines dt "-" lw 6 lc rgb "#d53e4f" notitle, "plot_ev_14.txt" u 2:($3*27.211) with lines dt "." lw 6 lc rgb "#d53e4f" notitle


# Plotting the wavefunctions is a whole big thing.
au2ev = 27.21
set samples 500
set isosamples 500
set ztics -0.02,0.02,0.02
unset colorbox
set xrange [-5.5:5.5]
set yrange [-1.7:2.3]
set zrange [-0.35:0.35] #-0.02*au2ev:0.02*au2ev]

set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
set zlabel "{/Times-Italic E}/au"
set ylabel "{/Times-Italic q_c}/au"
set xlabel "{/Times-Italic q_t}/au"
set xyplane -0.02
set ytics -1.5,1.5,1.5
set ztics -0.2,0.2,0.2
set size 1.0,1.0
splot f(x,y)*au2ev with pm3d notitle, g(x,y)*au2ev with pm3d notitle

set terminal postscript enhanced landscape colour font "Times-Roman,46"
reset
set xtics -6,3,6
set ytics -3,3,3
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
splot "wf_2d_0.txt" u ($2/Q0):($1/Q0):7 with pm3d at b notitle

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
splot "wf_2d_1.txt" u ($2/Q0):($1/Q0):5 with pm3d at b notitle

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
splot "wf_2d_1.txt" u ($2/Q0):($1/Q0):7 with pm3d at b notitle

unset surface
set cntrlabel start 5 interval 100 onecolor
set contour
set cntrparam levels incremental -1,0.005,50.0
splot f(x,y) with lines lc "black" notitle

unset multiplot
