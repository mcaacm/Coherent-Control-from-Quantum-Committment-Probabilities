PROGRAM nsq


USE parameters
USE prequel
USE rk_ns_utils
USE set_up_H
USE qns

IMPLICIT NONE

! Density matrix, Hamiltonian which will be diagonalized
! right eigenvectors of H, inverse of evec matrix
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma, H, Hls
! Hamiltonian eigenvectors and inverse, matrices for secular Redfield propagation
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: evec, evec_inv
! Starting vectors for the lindblad jumps routine
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: translated_evecs
! Fourier coupling operator
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:,:) ::  couple_op
! Right eigenvalues of Hamiltonian 
COMPLEX*16, DIMENSION(n) :: eval
! Probability to begin jump run from a given vector when more than one is possible
REAL*8, DIMENSION(n) :: estate_probs
! Indexes, number of steps to run in secular redfield
INTEGER :: i, j, k, l
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: G_pls_sig  ! Temporaries for Redfield
REAL*8 :: sl, dsum, rate_sum
TYPE(lb_op), DIMENSION(:), ALLOCATABLE :: operators
REAL*8, DIMENSION(n) :: evals
INTEGER :: t_len, nfound
INTEGER, DIMENSION(3) :: dests

INTEGER :: t_len_secular
TYPE(lb_op), DIMENSION(:), ALLOCATABLE :: operators_secular
TYPE(blocked_rho) :: brho_secular

COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: operator_rhos

COMPLEX*16, DIMENSION(n,1) :: wavefunction_n  ! For outputting optimal wavefunctions
COMPLEX*16, DIMENSION(m,1) :: wavefunction_m
REAL*8, DIMENSION(2*m) :: wavefunction_m_out

REAL*8 :: x, y, z, r, t, p
TYPE(blocked_rho) :: brho

COMPLEX*16, DIMENSION(n,n) :: U, Vt, secular_Hls, nonsecular_Hls
REAL*8, DIMENSION(n) :: sigma_svd
INTEGER :: nsvd

REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: solved_basis
REAL*8, DIMENSION(3) :: committors

INTEGER, DIMENSION(3) :: state_ids
COMPLEX*16, DIMENSION(n,n) :: sigma_original

REAL*8, DIMENSION(:,:), ALLOCATABLE :: interblock_fluxes  ! For flux graph plots

! Allocate memory
ALLOCATE(sigma(n,n))
ALLOCATE(solved_basis(n,n,3))
ALLOCATE(H(n,n))
ALLOCATE(Hls(n,n))
ALLOCATE(evec_inv(n,m))
ALLOCATE(translated_evecs(n,n))
ALLOCATE(evec(m,n))
! Allocate global storage found in prequel
ALLOCATE(couple_op(n,n,nbath))
ALLOCATE(G_min_global(n,n,nbath))
ALLOCATE(G_pls_global(n,n,nbath))
ALLOCATE(G_global(n,n,nbath))
! Allocate global storage in Setup_H
ALLOCATE(Q_x_g(n,n))
ALLOCATE(Q_y_g(n,n))
ALLOCATE(H_g(n,n))
! These are never deallocated but it doesn't matter;
! They're important until the end


IF (DEBUG .EQV. .TRUE.) THEN
  WRITE(*,*) "Running with following parameters: omega_c",  omegac, "eta", eta, "bath type", bath_type, &
           "dim_nc", dim_nc, "dim_nt", dim_nt, "m", m, "n", n, "omegas_c", omegas_c, "omegas_t", omegas_t, &
           "kappa_t", kappa_t, "exe", exe, "lambda_s", lambda_s
END IF


H = (0.0d0, 0.0d0)
evec = (0.0d0, 0.0d0)
sigma = (0.0d0,0.0d0)
eval = (0.0d0,0.0d0)

CALL setup(H,sigma,evec,eval,translated_evecs,estate_probs,couple_op)
CALL normalize_trace(n,sigma) 
CALL set_G(couple_op)

WRITE(*,*) "System eigenvalues are ", eval
! Copy and invert the eigenvector matrix
evec_inv = CONJG(TRANSPOSE(evec))
WRITE(*,*) "Starting sigma is "
CALL print_mat(sigma)
WRITE(*,*) ""

H_g = 0.0d0
DO i = 1, n
  H_g(i,i) = eval(i)
END DO

evals = eval
sl = (evals(11) - evals(10))*0.90d0
CALL find_sb_jops(sl,evals,operators,t_len,couple_op,brho)
ALLOCATE(interblock_fluxes(brho%num_blocks,brho%num_blocks))

WRITE(*,*) "Coupling operator 1 in Eigenbasis "
CALL print_mat(couple_op(1:n,1:n,1))
WRITE(*,*) "Coupling operator 2 in Eigenbasis "
CALL print_mat(couple_op(1:n,1:n,2))

! Printing out blocking and transfer rate information 
DO i = 1, t_len
  DO j = 1, nbath
    WRITE(*,*) "For operator ", i, "bath", j, "omega", operators(i)%omega, "occupation", operators(i)%occupation(j), &
      "Source block", operators(i)%sblock, "Destination block", operators(i)%dblock
    CALL print_mat(operators(i)%A(1:n,1:n,j))
    WRITE(*,*) "got thpl", operators(i)%thpl(j)
  END DO
END DO
DO i = 1, n  ! For each eigenstate source
  rate_sum = 0.0d0
  DO j = 1, t_len  ! For each operator
    DO k = 1, nbath  ! For each bath
      DO l = 1, n  ! For each eigenstate destination 
        IF (ABS(operators(j)%A(l,i,k)) .GT. 1d-10 .AND. &
            ABS(operators(j)%A(l,l,k)) .LT. 1d-10) THEN  ! Occupation confirmed for non-dephasing operator
          WRITE(*,*) "Jump between ", i, "and", l, "with gamma", operators(j)%rate(k), "total rate", &
            ABS(operators(j)%A(l,i,k))*operators(j)%rate(k)
          rate_sum = rate_sum + ABS(operators(j)%A(l,i,k))*operators(j)%rate(k)
        END IF
      END DO
    END DO
  END DO
  WRITE(*,*) "Total departure rate for ", i, "is", rate_sum
END DO



! Assemble Hls
CALL assemble_Hls(t_len,operators,Hls,H)
WRITE(*,*) "Hls assembled"
CALL print_mat(Hls)


! Propagation and flux calculation over 50000 au optimized to ES 3
OPEN(16,FILE="nonsecular_relax_optimal_to_es1.txt")
OPEN(17,FILE="dmats_to_es_1.txt")
sigma = 0.0d0
sigma(13,13) = 0.8658d0
sigma(14,14) = 0.1342d0
sigma(13,14) = 0.34d0
sigma(14,13) = 0.34d0
interblock_fluxes = 0.0d0
DO i = 1, 50000
  CALL rk4_prop(dt,sigma,Hls,t_len,operators,brho,interblock_fluxes)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4))), &
    REAL(REAL(sigma(5,5))), REAL(REAL(sigma(6,6))), REAL(REAL(sigma(7,7))), REAL(REAL(sigma(8,8))), &
    REAL(REAL(sigma(9,9))), REAL(REAL(sigma(10,10))), REAL(REAL(sigma(11,11))), REAL(REAL(sigma(12,12))), &
    REAL(REAL(sigma(13,13))), REAL(REAL(sigma(14,14))), REAL(REAL(sigma(15,15)))
  IF (MOD(i,50) .EQ. 0) THEN  ! Print matrix to file for plotting
    WRITE(17,*) "t=", i*dt
    CALL print_mat(sigma,17)
    WRITE(17,*) ""
  END IF
END DO
CLOSE(16)
CLOSE(17)
WRITE(*,*) "During optimal propogation to eigenstate 1 the following interblock fluxes were recorded: "
DO i = 1, brho%num_blocks
  WRITE(*,*) interblock_fluxes(1:brho%num_blocks,i)
END DO


! Propagation and flux calculation over 50000 au optimized to ES 1
interblock_fluxes = 0.0d0
OPEN(16,FILE="nonsecular_relax_optimal_to_es3.txt")
OPEN(17,FILE="dmats_to_es_3.txt")
sigma = 0.0d0
sigma(14,14) = 0.8658d0
sigma(13,13) = 0.1342d0
sigma(13,14) = -0.34d0
sigma(14,13) = -0.34d0
interblock_fluxes = 0.0d0
DO i = 1, 50000
  CALL rk4_prop(dt,sigma,Hls,t_len,operators,brho,interblock_fluxes)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4))), &
    REAL(REAL(sigma(5,5))), REAL(REAL(sigma(6,6))), REAL(REAL(sigma(7,7))), REAL(REAL(sigma(8,8))), &
    REAL(REAL(sigma(9,9))), REAL(REAL(sigma(10,10))), REAL(REAL(sigma(11,11))), REAL(REAL(sigma(12,12))), &
    REAL(REAL(sigma(13,13))), REAL(REAL(sigma(14,14))), REAL(REAL(sigma(15,15)))
  IF (MOD(i,50) .EQ. 0) THEN  ! Print matrix to file for plotting
    WRITE(17,*) "t=", i*dt
    CALL print_mat(sigma,17)
    WRITE(17,*) ""
  END IF
END DO
CLOSE(16)
CLOSE(17)
WRITE(*,*) "During optimal propogation to eigenstate 3 the following interblock fluxes were recorded: "
DO i = 1, brho%num_blocks
  WRITE(*,*) interblock_fluxes(1:brho%num_blocks,i)
END DO

! Relaxation purely from this eigenstate, 14
OPEN(16,FILE="nonsecular_relax_from_14.txt")
sigma = 0.0d0
sigma(14,14) = 1.0d0
interblock_fluxes = 0.0d0
DO i = 1, 50000
  CALL rk4_prop(dt,sigma,Hls,t_len,operators,brho,interblock_fluxes)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4))), &
    REAL(REAL(sigma(5,5))), REAL(REAL(sigma(6,6))), REAL(REAL(sigma(7,7))), REAL(REAL(sigma(8,8))), &
    REAL(REAL(sigma(9,9))), REAL(REAL(sigma(10,10))), REAL(REAL(sigma(11,11))), REAL(REAL(sigma(12,12))), &
    REAL(REAL(sigma(13,13))), REAL(REAL(sigma(14,14))), REAL(REAL(sigma(15,15)))
END DO
CLOSE(16)


! Relaxation purely from this eigenstate, 13
OPEN(16,FILE="nonsecular_relax_from_13.txt")
sigma = 0.0d0
sigma(13,13) = 1.0d0
interblock_fluxes = 0.0d0
DO i = 1, 50000
  CALL rk4_prop(dt,sigma,Hls,t_len,operators,brho,interblock_fluxes)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4))), &
    REAL(REAL(sigma(5,5))), REAL(REAL(sigma(6,6))), REAL(REAL(sigma(7,7))), REAL(REAL(sigma(8,8))), &
    REAL(REAL(sigma(9,9))), REAL(REAL(sigma(10,10))), REAL(REAL(sigma(11,11))), REAL(REAL(sigma(12,12))), &
    REAL(REAL(sigma(13,13))), REAL(REAL(sigma(14,14))), REAL(REAL(sigma(15,15)))
END DO
CLOSE(16)

sigma = 0.0d0
brho%rho = sigma
ALLOCATE(brho%committor_basis(n,n,3))
state_ids(1) = 1
state_ids(2) = 2
state_ids(3) = 3

dests(1) = 1
dests(2) = 2
dests(3) = 3

! Solve the committor basis.
solved_basis = 0.0d0
state_ids(1) = 1
state_ids(2) = 2
state_ids(3) = 3
CALL solve_basis(dt,n,0.01d0,Hls,t_len,operators,3,state_ids,solved_basis,1,15)
WRITE(*,*) "Basis solved: "
WRITE(*,*) "Basis entry (i,j) is for the imaginary part of coherences (i,j) and (j,i) if above the diagonal,"
WRITE(*,*) "and the real part if below the diagonal."
DO i = 1, n
   WRITE(*,*) ""
   DO j = 1, n
     WRITE(*,*) "i,j", i, j, "committors", solved_basis(i,j,1:3)
   END DO
END DO

! Committors for optimal wavefunction to 1 are verified
sigma = 0.0d0
sigma(13,13) = 0.8658d0
sigma(14,14) = 0.1342d0
sigma(13,14) = 0.34d0
sigma(14,13) = 0.34d0
wavefunction_n = 0.0d0
wavefunction_n(13,1) = SQRT(0.8658d0)
wavefunction_n(14,1) = SQRT(0.1342d0)
wavefunction_m = MATMUL(evec,wavefunction_n)
wavefunction_m_out = org_out(m,wavefunction_m)
OPEN(17,FILE="optimal_to_1_primitive_basis.txt")
WRITE(17,*) wavefunction_m_out
CLOSE(17)
CALL abc_committor(dt,n,sigma,Hls,t_len,operators,3,state_ids,committors,0.01d0)
WRITE(*,*) "Commitors P(1,2,3) via basis approach for coherent sigma: ", committors
CALL calculate_committors_by_basis(n,sigma,solved_basis,3,committors)
WRITE(*,*) "Via brute force propagation: ", committors
WRITE(*,*) "If those do not match within tolerance something is broken."

! Committors for optimal wavefunction to 3 are verified
sigma = 0.0d0
sigma(14,14) = 0.8658d0
sigma(13,13) = 0.1342d0
sigma(13,14) = -0.34d0
sigma(14,13) = -0.34d0
wavefunction_n = 0.0d0
wavefunction_n(14,1) = -SQRT(0.8658d0)
wavefunction_n(13,1) = SQRT(0.1342d0)
wavefunction_m = MATMUL(evec,wavefunction_n)
wavefunction_m_out = org_out(m,wavefunction_m)
OPEN(17,FILE="optimal_to_3_primitive_basis.txt")
WRITE(17,*) wavefunction_m_out
CLOSE(17)
CALL abc_committor(dt,n,sigma,Hls,t_len,operators,3,state_ids,committors,0.01d0)
WRITE(*,*) "Committors P(1,2,3) via basis approach for incoherent sigma: ", committors
CALL calculate_committors_by_basis(n,sigma,solved_basis,3,committors)
WRITE(*,*) "Via brute force propagation: ", committors
WRITE(*,*) "If those do not match within tolerance something is broken."



END PROGRAM
