! "Coherent control from quantum committment probabilities"
! Michelle C Anderson, Amro Dodin, Thomas P. Fay, David T. Limmer
! Copyright 2024

! This file is part of QI_COMM
! QI_COMM is free software: you can redistribute it and/or modify it under the terms of the GNU General 
! Public License as published by the Free Software Foundation, either version 3 of the License, 
! or (at your option) any later version.
! QI_COMM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with QI_COMM
! If not, see <https://www.gnu.org/licenses/>.



! All parameters should be speicifed in atomic units to avoid any
! complications.
! Parameters to specify how the run proceeds
MODULE parameters
IMPLICIT NONE

! In all cases save where otherwise noted, 1 is a 'true' or 'yes' on a binary
! switch and '0' is 'no' or 'false' but usually any number except 1 is ignored.

! Option to print extra information to deal with problems
LOGICAL, PARAMETER :: DEBUG = .TRUE.
LOGICAL, PARAMETER :: RESTART = .FALSE.

! Basic system setup parameters
INTEGER, PARAMETER :: CI = 0  ! Conical intersection
! For the case of a shin_metiu model
INTEGER, PARAMETER :: SM = 1


! What kind of run to perform (CI is the only setup included in this code, so leave unchanged)
INTEGER, PARAMETER :: run_type = SM
! Should the CI off-diagonal coupling be \lambda*|Q_c| rather than \lambda*Q_c 
INTEGER, PARAMETER :: abs_coupling = 0

! Time step for rk4 routine
REAL*8, PARAMETER :: dt = 0.100d0 !duration / (ntimes - 1)

! Constants and conversions
REAL*8, PARAMETER :: pi = 3.1415926d0
REAL*8, PARAMETER :: hbar = 1.0d0 

! Bath parameters
! Number of baths which will be coupled to some coordinate
INTEGER, PARAMETER :: nbath = 2
! Cutoff parameter in spectral density for each bath
REAL*8, PARAMETER, DIMENSION(nbath) :: omegac = (/0.050d0,0.050d0/)
! Strength parameter in spectral density
REAL*8, PARAMETER, DIMENSION(nbath) :: eta = (/0.20d0,0.20d0/) !,0.025d0/) !0.001822766d0/) 
! Debye bath is 1, Ohmic is 0
! ohmic is the default if no valid type is given
INTEGER, PARAMETER, DIMENSION(nbath) :: bath_type = (/0,0/) !,0/) !,1/)
! Temperature 1/(kb*T). To indicate zero temperature set beta = -1.0d0;
! any negative value will do but please use -1.0d0
REAL*8, PARAMETER, DIMENSION(nbath) :: beta = 1052.584412992859d0*(/1.0d0,1.0d0/)

! SM parameters; Initialization WILL ALWAYS USE
! A SINGLE EIGENSTATE; IT WON'T EVEN CHECK THE SWITCH
INTEGER, PARAMETER :: dim_qc = 60 !20 !50  ! Harmonic coordinate dimension; must be even
INTEGER, PARAMETER :: dim_R = 81 ! DVR coordiante dimension; must be odd
REAL*8, PARAMETER :: dx_R = 0.030d0  ! DVR delta x placement information
REAL*8, PARAMETER :: dx_Qc = 0.10d0  ! DVR delta x placement information
REAL*8, PARAMETER :: omegasc = 0.0264d0*0.96d0  ! Harmonic coordinate spring constant
REAL*8, PARAMETER :: mass_R = 1836.0d0  ! DVR coordinate mass
! Parameters for the potential energy surface and H
! for the adiabatic QED Hamiltonian:
! H = P^2/2M + E(R) + H_vib + pc^2/2 + wc^2/2 *
! (qc + sqrt(2/wc) *etac * mu(R))^2 
! mu(R) = coeff_v*tanh(y*R) + coeff_z*R
! E(R) = cob**4/(16 ceb)*R^4 - 0.5*cob*R**2 - c*R**3
REAL*8, PARAMETER :: cob = 0.8d0
REAL*8, PARAMETER :: ceb = 0.05d0
REAL*8, PARAMETER :: coeff_c = 4d-3
REAL*8, PARAMETER :: etac = 0.2000d0 !2193135d0 !0.005d0
REAL*8, PARAMETER :: coeff_v = -1.7d0
REAL*8, PARAMETER :: coeff_y = 3.0d0
REAL*8, PARAMETER :: coeff_z = 0.60d0 

! Conical Intersection Parameters
INTEGER, PARAMETER :: dim_nc = 20 !40   ! Basis set size, coupling coodinate
INTEGER, PARAMETER :: dim_nt = 60 !90   ! Basis set size, tuning coordiante
REAL*8, PARAMETER :: omegas_t = 0.002279d0 !0.0043364174d0 ! Coupling coordinate frequency
REAL*8, PARAMETER :: omegas_c = 0.004116d0 !0.0027194482d0 ! Tuning coordinate frequency
REAL*8, PARAMETER, DIMENSION(2) :: kappa_t = (/-0.006836d0, 0.006836d0/)!-0.105d0*ev2au,0.149d0*ev2au/)
REAL*8, PARAMETER, DIMENSION(2) :: exe = (/-0.00139d0,0.00139d0/) !0.144792242d0, 0.177866612d0/) ! Base energies of both wells
REAL*8, PARAMETER :: lambda_s = 0.00091153d0

! Size parameters
! Dimension of density matrix; must be an even number
INTEGER, PARAMETER :: n = 4 !dim_nc*dim_nt*2 !400 ! Truncated basis size, must be less than or equal to m
INTEGER, PARAMETER :: m = dim_qc*dim_R !dim_nc*dim_nt*2 !dim_qc*dim_R ! Basis size for initial diagonalization (dim_nc*dim_nt*2 for a 2-surface CI)

! Usually the initial state is changed from these default
! initializations before the run occurs and they are 
! irrelevant. Double check to make sure your initial 
! density matrix is what you want. 
! The default initialization is vertical excitation
! Set init_stat_therm to 1 to initialize in an eigenstate
! Populate a single eigenstate
INTEGER, PARAMETER :: init_stat_therm = 0
! Eigenstate to initialize all density in if init_stat_therm = 1
INTEGER, PARAMETER :: init_es = 4
! Initialize in either well 1 of the CI or well 2;
! Well 2 is the one used for vertical excitation usually
INTEGER, PARAMETER :: init_well = 2

! Plotting bounds; the eigenstates between 1&2 and 3&4 are plotted
INTEGER, PARAMETER :: bound1 = 1
INTEGER, PARAMETER :: bound2 = 10
INTEGER, PARAMETER :: bound3 = 11
INTEGER, PARAMETER :: bound4 = 15

END MODULE parameters
