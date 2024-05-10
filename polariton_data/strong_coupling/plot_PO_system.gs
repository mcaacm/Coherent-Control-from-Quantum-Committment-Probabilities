# This script is for plotting potential energy
# surface functions for a polariton system. It
# plots energies of various different energy
# eigenstates alongside the PES.

set terminal postscript enhanced colour font "Times-Roman,44"
set output "energy_profile.ps"

# Set up polariton parameters
Q0 = 1.0
coeff_v = -1.7
coeff_y = 3.0
coeff_z = 0.6
cob = 0.8
ceb = 0.05
coeff_c = 4e-3
etac = 0.20
omegasc = 0.0264*0.96

au2ev = 27.21
aon = 1.4

set xlabel "{/Times-Italic R}/au"
set ylabel "{/Times-Italic E}/eV"

set xrange [-1:1]
set yrange [-0.052*au2ev + aon:0.01*au2ev + aon]
set ytics 0.1,0.5,1.6

# Dipole
mu(x) = coeff_v*tanh(coeff_y*x) + coeff_z*x
set xlabel "{/Times-Italic R}/au"
# Potential surface
E(x) = ((cob*x)**4.0)/(16*ceb) - 0.5*(cob*x)**2.0 - coeff_c*x**3.0 + (omegasc**2.0/2.0)*(2/omegasc)*(etac*mu(x))**2

plot E(x)*au2ev + aon with lines lw 6 lc "black" notitle, "plot_ev_b1.txt" u 2:($3*au2ev + aon) with lines lw 6 lc "web-blue" notitle, "plot_ev_b2.txt" u 2:($3*au2ev + aon) with lines lw 6 lc "dark-blue" notitle, "plot_ev_b3.txt" u 2:($3*au2ev + aon) with lines lw 6 lc "spring-green" notitle, "plot_ev_b4.txt" u 2:($3*au2ev + aon) with lines lw 6 lc "dark-green" notitle
reset
set xrange [-1:1]
set ylabel "{/Symbol-Oblique m}/au"
set xlabel "{/Times-Italic R}/au"
plot mu(x) with lines lw 6 lc "black" notitle


f(x,y) = ((cob*x)**4.0)/(16*ceb) - 0.5*(cob*x)**2.0 - coeff_c*x**3.0 + 0.5*(omegasc**2.0)*(y + sqrt(2.0/omegasc)*mu(x)*etac)**2.0 # 

unset colorbox
set terminal postscript enhanced colour font "Times-Roman,34"
#set size 80,80
set samples 1000
set isosamples 1000
set xrange [-1:1]
set yrange [-7:20]
set ytics -5,10,20
set zrange [-0.05*au2ev + aon: 1.6] #0.04*au2ev + aon + 0.1]
set xyplane at -0.06
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')

set ztics 0.1,0.5,1.6 #-0.05*au2ev,0.02*au2ev,0.04*au2ev
set cbtics 0.1,0.5,1.6 #-0.05*au2ev,0.02*au2ev,0.04*au2ev

set xlabel "{/Times-Italic R}/au"
set ylabel "{/Times-Italic q_c}/au"
set zlabel "{/Times-Italic E}/au"

unset xlabel
unset ylabel
unset zlabel

set xlabel "{/Times-Italic R}/au"
set ylabel "{/Times-Italic q_c}/au"
set zlabel "{/Times-Italic E}/au"

splot f(x,y)*au2ev + aon with pm3d notitle
