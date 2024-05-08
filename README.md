The following contains code and data to accompany 'Coherent control from 
quantum committment probabilities'
Michelle C Anderson, Amro Dodin, Thomas P. Fay, David T. Limmer
Copyright 2024
This file is part of QI_COMM
QI_COMM is free software: you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation, either version 3 of the License, 
or (at your option) any later version.
QI_COMM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with QI_COMM
If not, see <https://www.gnu.org/licenses/>.


PAPER DATA

The data for the conical intersection model is found in 
'conical_intersection_data.' It contains a separate
'README.md' file with more detailed descriptions of where
to find paticular pieces of supporting data.


The data for the polariton model is found in 'polariton_data.' 
This folder contains a separate 'README.md' file with more
detailed descriptions of the contents.

--------------------------------------------------------------------------------------------

PROGRAM SUMMARY

QI_COMM contains code to run the partial secular Redfield equations
in the GKSL master equation formulation. Code to calculate committor
values is found in 'nsq.f90' whereas 'nsq_main.f90' is a driver program.
Integration routines for finding the Redfield rates are in 'rkns.f90'
and 'prequel.f90' alongside many utilities for setting up and plotting
systems. System parameters are found in 'params.f90.' The code to set
up the system Hamiltonian and additional operators are found in 
'setup_h.f90.'

NONSECULAR BLOCK SETUP

The nonsecular block binning process is manual in this software and
limited in scope. The user must designate a parameter in 'nsq_main.f90'
which is 'sl' as equal to the minimum energy gap required to separate
nonsecular blocks. This parameter will be used to select nonsecular blocks.
An additional paremeter in 'nsq.f90', 'hard_separate_bottom' will
force the bottom two eigenstates in the system to be members of different
nonsecular blocks if it is set to .TRUE. These are the only two options
available for adjusting blocking in this software. 

The user must inspect the nonsecular blocks which are produced and determine 
whether they are reasonable.  Generally, if the rate to leave any eigenstate under 
the full secular approximation is equal or larger than the gap between that eigenstate 
and any other state, those two states should be in the same nonsecular
block. If the rate is not equal to or larger than the gap, they should not
be in the same nonsecular block. 

In the case that a block structure cannot be assembled  because these
two requirements cannot be satisfied at the same time, it is not
recommended to use the partial secular approximation without some 
means of verifying the veracity of results against a reliable method.


PROGRAMS

The programs can be compiled from the given Makefile with

'make'

Most programs require Lapack and BLAS. Quadpack is also required
but the necessary file is included and compiled by make. Quadpack
is opensource software licensed as specified in 'quadpack.f90'.

--------------------------------------------------------------------------------------------

Programs were constructed and tested on Mac OS Sonoma and compiled with gfortran 
8.2.0. The makefile settings assume that gfortran, Lapack and BLAS are available.
To change the compiler or library linking settings, edit 'Makefile'.

To compile and run a calculation, select the proper 'run_type' in 
'params.f90,' modify 'nsq_main.f90' to carry out the desired tasks 
and run:

make  
./nsq_main

--------------------------------------------------------------------------------------------

Programs in General:

Programs were constructed and tested on Mac OS Mojave and compiled with gfortran 8.2.0. The
makefile settings assume that gfortran, Lapack and BLAS are available.
To change the compiler or library linking settings, edit 'Makefile'.

Note, all programs that call a setup function will write a file 'evic.txt'
which includes information about the average characteristics of energy eigenvectors of in the format
'(coupling coordinate) (tuning coordinate) (energy) (diabatic population 1) (diabatic population 2)'
for the conical intersection case and '(average matter coordinate)', '(average light coordinate)',
'(energy)'. Energies are always in atomic units. The conical intersection positions are
in dimensionless coordinates. The polariton system has positions in atomic units again.

All programs will print a file 'restart.dat' which includes eigenvalues and eigenvectors
which could be used to run another calculation with the same parameters for the 
Hamiltonian without going through the diagonalization process, but 'restart.dat' is
very big.

All programs also produce a file, 'plot_ev.txt' which merely contains 
'(eigenvecor number) (-100 or 100) (eigenvector energy)' which can be used to plot
the energies of the eigenvectors by the gnuplot scripts.

It is recommended to use atomic units for all calculations. Although the value of reduced
Planck's constant can be changed in the parameters file, this option has not been
tested in several years and should be considered deprecated.

The size of the primitive basis is 'm' whereas the trimmed basis, number of eigenstates
involved in the calculation, is 'n.' This and all other simulation parameters
are explained and specified in 'params.f90.'

--------------------------------------------------------------------------------------------

nsq_main

This program calls setup frunctions then organizes into nonsecular blocks based on the size
of 'sl' by calling 'find_sb_jops.' It returns operators and the number of operators, 't_len,'
as well as an object with information about the block structure of rho, 'brho.'

It will then print information about the jump operators including their source and destination
blocks, omegas and which bath produces them as well as the associated one sided Fourier
transform, thpl. It will then print out information about overall jump rates from
one eigenstate to other eigenstates, which may be necessary to determine if blocking
has been done properly.

After printing this information about the nonsecular blocks and rates, actual 
calculations are performed as detailed in the body of the code.

--------------------------------------------------------------------------------------------

wfn_plot

This program will plot wavefunctions and wavefunction densities from a density matrix.

It will begin by opening the file 'in_wfns.txt' and it will attempt to read lines of
length 2*n where n is the trimmed basis size and the real and imaginary parts of each
basis state are given separately. It will put plots of these wavefunctions in
'wfn_plots/wf_2d_#i.txt' where '#i' is the number of wavefunctions read and plotted so
far. In the case that the run is an adiabatic polariton, the output will have three columns,
the wavefunction density then the real part then the imaginary part. In the case that
the run is the conical intersection with two electronic states, there are six columns. The
first two are the electronic state densities, then the real and imaginary part of the
first electronic state, then the real and imaginary part of the second electronic state. 

The program will do the same for eigenstates between bound1 and bound2 and bound3 and
bound4, found in 'params.f90' and write the output to 'wfn_plots/ef_2d_#i.txt'

In the case that a conical intersection is selected, a density matrix can be plotted.
This code is not configured to do polariton density matrices, so do not use in that case.
The program will attempt to open 'density_matrices.txt' and read in an nxn density 
matrix. It will then diagonalize the density matrix and plot a weighted sum
of the resulting eigenvectors with the weights being the eigenvalues, excluding any
with weights less than 0.001. It will print the matrices to 'wfn_plots/pf_2d_#i.txt'
and it will read as many density matrices as it can find. 
