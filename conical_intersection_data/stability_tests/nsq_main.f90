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

COMPLEX*16, DIMENSION(n,n) :: full_basis_transform, full_basis_transform_inverse, sigma_escape
COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: operator_rhos

REAL*8 :: x, y, z, r, t, p
TYPE(blocked_rho) :: brho

COMPLEX*16, DIMENSION(n,n) :: U, Vt
REAL*8, DIMENSION(n) :: sigma_svd
INTEGER :: nsvd

REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: solved_basis
REAL*8, DIMENSION(3) :: committors
REAL*8, DIMENSION(2*m,1) :: test_vector_raw  ! Read in eigenvector from another system
REAL*8, DIMENSION(2*m) :: temp_vector  ! To write
COMPLEX*16, DIMENSION(m,1) :: test_vector, local_test_vector  ! Read in vector from another system
COMPLEX*16, DIMENSION(n,1) :: test_vectore, local_test_vectore  ! Test vector in eigenbasis
COMPLEX*16, DIMENSION(n,n) :: test_matrix  ! Read and write dmat for basis interconversion
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: full_test_sigma  ! Write pseudo density matrix for coherence dynamics testing
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: test_dmat_full  ! Read in pseudo density matrix for coherence dynamics testing
INTEGER :: reason  ! input output success flag

REAL*8 :: r13, r14, c14, c13, c1314  ! Parameters to solve for the optimal wavefunction of 13/14 to go to ES 3
REAL*8 :: r13sf, r14sf, ref  ! Check for a possible sign flip

INTEGER, DIMENSION(3) :: state_ids


! Allocate memory
ALLOCATE(test_dmat_full(m,m))
ALLOCATE(sigma(n,n))
ALLOCATE(full_test_sigma(m,m))
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

WRITE(*,*) "evals are ", eval
! Copy and invert the eigenvector matrix
evec_inv = CONJG(TRANSPOSE(evec))
WRITE(*,*) "Starting sigma is "
CALL print_mat(sigma)
WRITE(*,*) ""



evals = eval
!CALL sort_ls(evals,operators,t_len,sl,couple_op)

test_vectore = MATMUL(evec_inv,test_vector)
test_matrix = MATMUL(CONJG(test_vectore),TRANSPOSE(test_vectore))

sl = (evals(11) - evals(10))*0.90d0
CALL find_sb_jops(sl,evals,operators,t_len,couple_op,brho)
!CALL find_sb_jops(1.55d-4,evals,operators,t_len,couple_op)

WRITE(*,*) "Coupling operator 1 in Eigenbasis "
CALL print_mat(couple_op(1:n,1:n,1))
WRITE(*,*) "Coupling operator 2 in Eigenbasis "
CALL print_mat(couple_op(1:n,1:n,2))

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



!operators(t_len)%A = 0.0d0  ! Turn on or off global dephasing

! Assemble Hls
CALL assemble_Hls(t_len,operators,Hls,H)
WRITE(*,*) "Hls assembled"
CALL print_mat(Hls)

ALLOCATE(brho%committor_basis(n,n,3))
state_ids(1) = 1
state_ids(2) = 2
state_ids(3) = 3
CALL solve_basis(dt,n,0.01d0,Hls,t_len,operators,3,state_ids,brho%committor_basis(1:15,1:15,1:3),1,15)

! Use solved basis to designate optimal test matrix entry
c14 = brho%committor_basis(14,14,3)
c13 = brho%committor_basis(13,13,3)
c1314 = brho%committor_basis(14,13,3) 
WRITE(*,*) "Found c14 ", c14, "c13", c13, "c1314", c1314
r14 = SQRT((c1314 - 2.0d0*c13)**2.0d0 / ((2.0d0*c14 - c1314)**2.0d0 + (c1314 - 2.0d0*c13)**2.0d0))
r13 = SQRT(1.0d0 - r14**2.0d0)
r13sf = -r13
r14sf = r14
WRITE(*,*) "Found r13 ", r13, "r14 ", r14, "relative signs undetermined"
test_matrix = 0.0d0
test_matrix(14,14) = r14**2.0d0
test_matrix(13,13) = r13**2.0d0
test_matrix(14,13) = r14*r13
test_matrix(13,14) = r14*r13
CALL calculate_committors_by_basis(n,test_matrix,brho%committor_basis,3,committors)
ref = committors(3)
WRITE(*,*) "Committor calculated for optimal wavefunction ", committors
test_matrix = 0.0d0
test_matrix(14,14) = r14sf**2.0d0
test_matrix(13,13) = r13sf**2.0d0
test_matrix(13,14) = r14sf*r13sf
test_matrix(14,13) = r14sf*r13sf
CALL calculate_committors_by_basis(n,test_matrix,brho%committor_basis,3,committors)
WRITE(*,*) "Committor calculated for sign flipped optimal wavefunction", committors
IF (committors(3) .LT. ref) THEN
  WRITE(*,*) "Using sign flip."
  r13 = r13sf 
  r14 = r14sf
END IF

test_vectore = 0.0d0
test_vectore(13,1) = r13
test_vectore(14,1) = r14
local_test_vectore = test_vectore
WRITE(*,*) "Optimal vector to write is ", test_vectore
test_vector = MATMUL(evec,test_vectore)
local_test_vector = test_vector
WRITE(*,*) "Optimal vector in primitive basis is ", test_vector
OPEN(16,FILE="optimal_evec_out.txt")
temp_vector = org_out(m,test_vector)
WRITE(16,*) temp_vector
CLOSE(16)
WRITE(*,*) "Transformed back is ", MATMUL(evec_inv,test_vector)


OPEN(16,FILE="optimal_evec.txt")
READ(16,*,IOSTAT=reason) (test_vector_raw(i,1), i=1, 2*m)
IF (reason .NE. 0) THEN
  WRITE(*,*) "Read error of optimal vector:", reason, "defaulting to 1.0 ...."
  test_vector_raw = 0.0d0
  test_vector_raw(1,1) = 1.0d0
END IF
DO i = 1, m 
  test_vector(i,1) = CMPLX(test_vector_raw(i*2 - 1,1),test_vector_raw(i*2,1))
END DO
CLOSE(16)
WRITE(*,*) "Read in raw test vector ", test_vector
test_vectore = MATMUL(evec_inv,test_vector)
WRITE(*,*) "Optimal vector read in in eigenbasis is ", test_vectore
test_matrix = MATMUL(CONJG(test_vectore),TRANSPOSE(test_vectore))
CALL abc_committor(dt,n,test_matrix,Hls,t_len,operators,3,state_ids,committors,0.01d0)
WRITE(*,*) "Unperturbed optimal wavefunction solved via abc_committor ", committors
CALL calculate_committors_by_basis(n,test_matrix,brho%committor_basis,3,committors)
WRITE(*,*) "Unperturbed optimal wavefunction solved via basis committor ", committors


WRITE(*,*) "Overlap between optimal read in and optimal local is :" , &
  MATMUL(CONJG(TRANSPOSE(local_test_vector)),test_vector)

WRITE(*,*) "Compare with overlap computed in truncated eigenbasis to confirm validity :", &
  MATMUL(CONJG(TRANSPOSE(local_test_vectore)),test_vectore)

END PROGRAM

