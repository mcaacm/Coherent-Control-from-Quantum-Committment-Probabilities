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
! Spin-Boson for testing cases
INTEGER, PARAMETER :: SB = 2
! 1d SM
INTEGER, PARAMETER :: S1 = 3
! For the Rodopsin style CI initialization
INTEGER, PARAMETER :: RC = 4
! For the Four Level system
INTEGER, PARAMETER :: FL = 5
! For the Thoss HO
INTEGER, PARAMETER :: HO = 6  ! Harmonic oscillator


! What kind of run to perform (CI is the only setup included in this code, so leave unchanged)
INTEGER, PARAMETER :: run_type = CI
! Should the CI off-diagonal coupling be \lambda*|Q_c| rather than \lambda*Q_c 
INTEGER, PARAMETER :: abs_coupling = 0

! When creating a markov state model, how many timesteps to run (tau = dt*msm_steps)
INTEGER, PARAMETER :: msm_steps = 1000

! Only set one of these options to 1 at any time
! Run the secular Redfield approximation (1)
INTEGER, PARAMETER :: secular = 1
! Run lindblad jump dynamics  (1)
INTEGER, PARAMETER :: lb_switch = 0

! Which eigenstates LB propogation should end with
INTEGER, PARAMETER :: flag_estate_1 = 1  ! First eigenstate to stop propogation (should be 1)
INTEGER, PARAMETER :: flag_estate_2 = 4  ! Second eigenstate to stop propogation (should not be 1)

! Number of trajectories for a lindblad run to carry out
INTEGER, PARAMETER :: ntraj = 10000
! Duration for the Redfield propagation (au)
REAL*8, PARAMETER :: duration = 500000.0d0
! Time step for rk4 routine
REAL*8, PARAMETER :: dt = 20.000d0 !duration / (ntimes - 1)
! How many iterations of redfield to do between printings of position/diabatic state
INTEGER, PARAMETER :: print_num = 1

! Constants and conversions
REAL*8, PARAMETER :: pi = 3.1415926d0
REAL*8, PARAMETER :: hbar = 1.0d0 

! Bath parameters
! Number of baths which will be coupled to some coordinate
INTEGER, PARAMETER :: nbath = 2
! Cutoff parameter in spectral density for each bath
REAL*8, PARAMETER, DIMENSION(nbath) :: omegac = (/0.010d0,0.010d0/)
! Strength parameter in spectral density
REAL*8, PARAMETER, DIMENSION(nbath) :: eta = (/0.10d0,0.10d0/) !,0.025d0/) !0.001822766d0/) 
! Debye bath is 1, Ohmic is 0
! ohmic is the default if no valid type is given
INTEGER, PARAMETER, DIMENSION(nbath) :: bath_type = (/0,0/) !,0/) !,1/)
! Temperature 1/(kb*T). To indicate zero temperature set beta = -1.0d0;
! any negative value will do but please use -1.0d0
REAL*8, PARAMETER, DIMENSION(nbath) :: beta = 1052.584412992859d0*(/1.0d0,1.0d0/)

REAL*8, PARAMETER :: factor = 1.0d0
! SM parameters; Initialization WILL ALWAYS USE
! A SINGLE EIGENSTATE; IT WON'T EVEN CHECK THE SWITCH
INTEGER, PARAMETER :: dim_qc = 40 !20 !50  ! Harmonic coordinate dimension; must be even
INTEGER, PARAMETER :: dim_R = 51 ! DVR coordiante dimension; must be odd
REAL*8, PARAMETER :: dx_R = 0.040d0  ! DVR delta x placement information
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
REAL*8, PARAMETER :: coeff_c = 5d-3
REAL*8, PARAMETER :: etac = 0.0000d0 !2193135d0 !0.005d0
REAL*8, PARAMETER :: coeff_v = -1.90249d0
REAL*8, PARAMETER :: coeff_y = 1.26426d0
REAL*8, PARAMETER :: coeff_z = 0.37044d0 


! Conical Intersection Parameters
INTEGER, PARAMETER :: dim_nc = 40   ! Basis set size, coupling coodinate
INTEGER, PARAMETER :: dim_nt = 90   ! Basis set size, tuning coordiante
REAL*8, PARAMETER :: omegas_t = 0.002279d0*0.95d0 !0.0043364174d0 ! Coupling coordinate frequency
REAL*8, PARAMETER :: omegas_c = 0.004116d0 !0.0027194482d0 ! Tuning coordinate frequency
REAL*8, PARAMETER, DIMENSION(2) :: kappa_t = (/-0.006836d0, 0.006836d0/)!-0.105d0*ev2au,0.149d0*ev2au/)
REAL*8, PARAMETER, DIMENSION(2) :: exe = (/-0.00139d0,0.00139d0/) !0.144792242d0, 0.177866612d0/) ! Base energies of both wells
REAL*8, PARAMETER :: lambda_s = 0.00091153d0

! Rhodopsin CI parameters
INTEGER, PARAMETER :: ne = 2
INTEGER, PARAMETER :: nq = 70
INTEGER, PARAMETER :: nphi = 71  ! MUST BE ODD
REAL*8, PARAMETER :: dphi = 2.0d0*pi/(nphi-1.0d0)
REAL*8, PARAMETER :: E1_r = 0.0911383d0
REAL*8, PARAMETER :: W0_r = 0.1322975d0
REAL*8, PARAMETER :: W1_r = 0.0400567d0
REAL*8, PARAMETER :: lambda_r = 0.006982369d0
REAL*8, PARAMETER :: mass_ro = 75.929d0
REAL*8, PARAMETER :: omega_r = 0.00698237d0
REAL*8, PARAMETER :: kappa_r = 0.00367493d0

! HO stuff
REAL*8, PARAMETER :: omegas = 0.05/27.211d0
REAL*8, PARAMETER, DIMENSION(2) :: kappa = (/-omegas,omegas/)


! Size parameters
! Dimension of density matrix; must be an even number
INTEGER, PARAMETER :: n = 15 !dim_nc*dim_nt*2 !400 ! Truncated basis size, must be less than or equal to m
INTEGER, PARAMETER :: m = dim_nc*dim_nt*2 !dim_nc*dim_nt*2 !dim_qc*dim_R ! Basis size for initial diagonalization (dim_nc*dim_nt*2 for a 2-surface CI)

! The default initialization is vertical excitation
! Set init_stat_therm to 1 to initialize in an eigenstate
! Populate a single eigenstate
INTEGER, PARAMETER :: init_stat_therm = 0
! Eigenstate to initialize all density in if init_stat_therm = 1
INTEGER, PARAMETER :: init_es = 4
! Initialize in either well 1 of the CI or well 2;
! Well 2 is the one used for vertical excitation usually
INTEGER, PARAMETER :: init_well = 2

! SB testing
REAL*8, DIMENSION(2) :: energies = (/0.0d0,0.0d0/)
REAL*8 :: coupling = 0.1d0*omegas

! Lindblad jumps tolerance for wfn norm to match a random number
! during a binary search
REAL*8, PARAMETER :: tol = 1.0d-5

! Plotting bounds
INTEGER, PARAMETER :: bound1 = 1
INTEGER, PARAMETER :: bound2 = 10
INTEGER, PARAMETER :: bound3 = 11
INTEGER, PARAMETER :: bound4 = 15

END MODULE parameters
