This folder contains the data for the polariton calculations. 
The strongest coupling case (figure 1 d in the paper) is in
'strong_coupling.' The moderate case(figure 1 c) is in
'medium_coupling' and the weakest case (figure 1 a) is in
'very_weak_coupling.'

Each folder contains a plotting file 'sphere.gs' which
defines the committor plance for both the partial secular
and fully secular cases at that coupling strength and
plots them together on a single sphere. Each also includes
'nsq_main.f90' and 'params.f90' which are setup and dynamics
propagation files respectively, with 'nsq_main.f90' also
solving for the committor in eigenbasis for each system
and printing out its conclusions. Other debugging information
is printed by these programs as well.

The 'strong_coupling' case contains a number of files titled
'wf_#x_#y_#z.txt' which indicate that these are population
dynamics initialized at (#x,#y,#z) on the Bloch sphere. The
suffix 'coherence_killed' indicates dynamics where coherences
were manually set to zero at each time step. The suffix
'secular' indicates that this is the secular limit of
dynamics. The strong coupling folder also includes the
'plot_PO_system.gs' script for plotting the potential
energy surafce of the system (figure 1 a) and the scripts
to produce population plots (figure 1 e and 1 f) which 
is 'plot_figure.gs.' There are also 'plot_ev#.txt' files
which just have energy eigenvalues in a format convenient
for plotting. 

All three folders contain subfolders called 'confirming_plane'
which include dynamics beginning from three different locations
on the committor planes (secular and nonsecular). The populations
of eigenstates 1 and 2 (g_1 and g_2 in the paper) are equal to
within tolerance after density has departed from the interfering
states 3 and 4 (e_1 and e_2 in the paper) to confirm that the
committor planes as defined in 'sphere.gs' are correct. The 
files 'nsq_main.f90' and parameter files 'params.f90' in these
folders are for running these dynamics.

Files which show population dynamics print out the eigenstate
populations in order g_1, g_2, e_1, e_2 and the first column
is the time in au.
