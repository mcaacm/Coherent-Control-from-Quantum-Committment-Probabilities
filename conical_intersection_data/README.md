This folder contains the data for the conical intersection
part of the paper. Most of the data is in 'optimals' where
most data for figure 3 and figure 2 can be found. The folder
'stability_tests' contains data for figures 4 and 5 which
determine how stable the optimization is towards slight
changes to the system.

Optimals contains the 'nsq_main.f90' file which contains
most instructions necessary for regenerating the data in
this folder. The files 'nonsecular_relax_optimal_to_es#.txt'
contain relaxation dynamics for the two optimized wavefunctions
bound for the designated eigenstates (\Psi_A and \Psi_B in
the paper). The files 'optimal_to_#_primitive_basis.txt' 
contain the expressions for the optimized wavefunctions
in the primitive harmonic oscillator basis. This is 
necessary in order to transfer them into calculations
for modified Hamiltonians (which will naturally have
a different eigenbasis). The plotting file 'plot_figure.gs'
is used to plot the population graphs for the 
optimized initial conditions (figure 3 c and d). This
will output in the file 'populations.ps.' The optimal
wavefunction heat map plotting data is found in 
'wf_2d_#.txt' where '#' is either 1 or 0 and the plotting
is handled by 'plot_CI_system.gs' with output going
to 'parabolas.ps.'

The file 'in_wfns.txt' holds the optimal initial 
wavefunctions in the eigenbasis. These can be read in
by 'wfn_plot' and translated into heat map plot files
which will be sent to the subfolder 'wfn_plots.'

The subfolder 'wfn_plots' contains various eigenvalues
in 'plot_ev_#.txt' files to make the PES with overlaid
eigenstate energies in figure 3 a. It also contains the
'wf_2d_#.txt' files which contain information for use
by 'plot_CI_system.gs' to make the heat maps for the
optimized initial wavefunctions in figure 3 b. This will
output in 'parabolas.ps.'

The folder 'stability_tests' contains folders labeled
'r#' where the '#' indicates the factor which multiplies
\omega_t in the perturbed system. Each of these folders
contains the 'nsq_main.f90' file necessary to calculate
overlap factors between the perturbed optimal 
wavefunctions and the unperturbed optimal wavefunctions. 
The 'optimal_evec.txt' contains the unperturbed optimal
wavefunction in the primitive basis which are read in
and then 'optimal_evec_out.txt' contains the optimal
wavefunction for the perturbed system. 

In the cases of 'r0.95' and 'r1.05,' plotting is done
and hence these are the two special examples shown.
These two cases contain subfolders 'wfn_plots' and
these contain 'wf_2d_0.txt' which is the optimal
wavefunction for the perturbed system, which was 
read in by 'wfn_plot' from the file 'in_wfns.txt' and
has been translated into a heat plot form for gnuplot.
The 'plot_CI_system.gs' gnuplot script will produce
a 'plot.ps' file which shows the real part of the 
wavefunction overlaid on energy contours.

The folder 'animation' contains code necessary to 
produce density matrices at 500 au increments which
are output to the file 'dmats_to_es_#.txt" where # 
is either eigenstate 1 or 3 depending on the final
relaxation of the optimized relaxation. These files
can be moved to 'density_matrices.txt' wich will 
be read in by 'wfn_plot' and output as heat map
plotting files in 'wfn_plots/pf_2d_#.txt'. Note
that only the density part is going to make sense
for attempting to plot these. Attempting to plot
the phase will not work, because the signs will
not stay consistent as the system evolves (this is
just a LAPACK thing that could, probably, be 
worked around with some method to control the sign).
This is not desirable and will result in
gibberish when trying to compare results at 
different times.

The folder 'relaxation_plots' contains 'relaxing_to_es_#'
subfolders which contain the generated density matrices
in a file 'density_matrice.txt' which contains the
two density matrices from 20000 au and 50000 au 
(figure 3 b and d) one after another. These are not
labeled because 'wfn_plot' will not read them in 
correctly if they are labeled as it will treat the
label as part of the wavefunction. The heat mapping
files produced from these are found in the subfolder
'wfn_plots' as detailed previously. The gnuplot 
script 'plot.gs' will convert them into the desired
full heat maps.  