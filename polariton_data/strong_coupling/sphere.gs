# I learned to set up the sphere plotting code via table assembly using
# the demo on https://stackoverflow.com/questions/61704503/plotting-points-on-a-sphere-with-gnuplot
set terminal postscript enhanced colour font "Times-Roman,24"
set output "sphere.ps"

reset session
set view equal xyz
#set view 60,180, 1, 1 #equal xyz
#set view equal xyz
set xyplane relative 0
R = 1
#unset tics
#unset border
unset colorbox

# sphere prototype data
set xlabel "x"
set ylabel "y"
set zlabel "z"
set parametric
set isosamples 35
set samples 50
set urange [0:2*pi]
set vrange [-pi/2:pi/2]
# Set up points on a sphere.
set table $Sphere
    splot R*cos(u)*cos(v), R*sin(u)*cos(v), R*sin(v)
unset table
unset parametric

# Defining the committor plane for the strong coupling case.
# 0 = ax + by + cz + d
# ax
a = -0.2172  
#by
b = -0.0422 
#cz
c = 0.4295  
#d
d = 0.50 - 0.5  


# Set the color palette and then do the fancy footwork required to make the
# plot look decent.
set palette defined (10 '#67001f', 20 '#b2182b',30 '#d6604d',40 '#f4a582', 50 '#fddbc7', 60 '#f7f7f7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac', 110 '#053061')
set xrange [-1.01:1.01]
set yrange [-1.01:1.01]
set zrange [-1.01:1.01]
set cbrange [-1.0:1.0]
plane(x,y) = (-a*x - b*y - d)/c
plane2(x,y) = (plane(x,y) < -0.31) ? -0.31 : 1/0
plane3(x,y) = -0.31 
set multiplot
set hidden3d
splot plane3(x,y) w pm3d notitle, plane(x,y) w pm3d notitle, plane2(x,y) w pm3d notitle
splot $Sphere u 1:2:3 w l lc "black" notitle, plane3(x,y) w l lc "black" notitle, plane(x,y) notitle linecolor "black", plane2(x,y) notitle linecolor "black"
unset multiplot


# Plot a sphere without a plane and with tics
set terminal postscript enhanced colour font "Times-Roman,34"
set xtics -1,1,1
set ytics -1,1,1
set ztics -1,1,1
unset ztics
splot $Sphere u 1:2:3 w l lc "black" notitle


# Plot a sphere without a plane and without tics
unset xtics
unset ytics
unset ztics
unset border
splot $Sphere u 1:2:3 w l lc "black" notitle

