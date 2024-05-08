# This file is used to plot the adiabatic PES for a conical 
# intersection, as well as the diabatic surfaces if necessary.

set terminal postscript enhanced colour font "Times-Roman,30"
set output "intersection.ps"

Q0 = 1.0

kappa1 = -0.006836 
kappa2 = 0.006836 
omegast = 0.002279
omegasc = 0.004116 
e1 = -0.00139 
e2 = 0.00139 
lambda = 0.00091153 

set xrange [-5.0:5.0]
set yrange [-0.35:0.05] 
set samples 10000
set ylabel "{/Times-Italic E}/eV"# font "Times-Roman,18"
set xlabel "{/Times-Italic Q_{t}}/au at {/Times-Italic Q_{c}}=0" #font "Times-Roman,18"
# Diabatic surfaces
a(x,y) = 0.5*omegasc*y**2 + e1 + kappa1*x + 0.5*omegast*x**2
b(x,y)= e2 + kappa2*x + 0.5*omegasc*y**2 + 0.5*omegast*x**2

# Adiabatic surfaces
f(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 - sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))
 g(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 + sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))

set samples 600
set isosamples 600
set ztics -0.015,0.015,0.015
unset colorbox
set xrange [-4:4]
set yrange [-2.5:2.5]
set zrange [-0.015:0.015]

set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
unset xlabel
unset ylabel
set xyplane -0.02
splot f(x,y) with pm3d notitle, g(x,y) with pm3d notitle

