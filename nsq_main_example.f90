PROGRAM nsq


USE parameters
USE prequel
USE rk_ns_utils
USE set_up_H
USE qns
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

! This file is a driver file for running the strong coupling polariton
! calculations that involve propagation from particular starting 
! conditions in both partial secular and fully secular cases. Sample
! input files for other kinds of calculations are found in the data
! subfolders alongside outputs.


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
REAL*8 :: sl,  rate_sum
TYPE(lb_op), DIMENSION(:), ALLOCATABLE :: operators
REAL*8, DIMENSION(n) :: evals
INTEGER :: t_len

INTEGER :: t_len_secular
TYPE(lb_op), DIMENSION(:), ALLOCATABLE :: operators_secular
TYPE(blocked_rho) :: brho_secular

TYPE(blocked_rho) :: brho
REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: solved_basis

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

WRITE(*,*) "evals are ", eval
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
!CALL sort_ls(evals,operators,t_len,sl,couple_op)

sl = (evals(4) - evals(3))*1.01d0
CALL find_sb_jops(sl,evals,operators,t_len,couple_op,brho)

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



sigma = 0.0d0
CALL xyz_bloch_to_rho(0.0d0,0.0d0,0.0d0,sigma(3:4,3:4))
WRITE(*,*) "sigma for x,y,z=(0.0,0.0,0.0)"
OPEN(16,FILE="wf_0.0_0.0_0.0.txt")
CALL print_mat(sigma)
DO i = 1, 100000
  CALL rk4_prop(dt,sigma,Hls,t_len,operators)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4)))
END DO
CLOSE(16)

sigma = 0.0d0
CALL xyz_bloch_to_rho(0.376d0,0.100d0,0.2d0,sigma(3:4,3:4))
WRITE(*,*) "sigma for x,y,z=(0.376,0.100,0.2)"
CALL print_mat(sigma)
OPEN(16,FILE="wf_0.376_0.100_0.2.txt")
DO i = 1, 100000
  CALL rk4_prop(dt,sigma,Hls,t_len,operators)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4)))
END DO
CLOSE(16)

sigma = 0.0d0
CALL xyz_bloch_to_rho(0.3954d0,0.0d0,0.2d0,sigma(3:4,3:4))
WRITE(*,*) "sigma for x,y,z=(0.3954,0.0,0.2)"
CALL print_mat(sigma)
OPEN(16,FILE="wf_0.3954_0.0_0.2.txt")
DO i = 1, 100000
  CALL rk4_prop(dt,sigma,Hls,t_len,operators)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4)))
END DO
CLOSE(16)


sl = (evals(4) - evals(3))*0.01d0
CALL find_sb_jops(sl,evals,operators_secular,t_len_secular,couple_op,brho_secular)
Hls = 0.0d0
CALL assemble_Hls(t_len,operators_secular,Hls,H)
WRITE(*,*) "Hls assembled"
CALL print_mat(Hls)
sigma = 0.0d0
CALL xyz_bloch_to_rho(0.2d0,0.2d0,-0.31d0,sigma(3:4,3:4))
OPEN(16,FILE="wf_0.3_0.1_-0.31_secular.txt")
DO i = 1, 100000
  CALL rk4_prop_coherence_kill(dt,sigma,Hls,t_len_secular,operators_secular)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4)))
END DO
CLOSE(16)


sigma = 0.0d0
CALL xyz_bloch_to_rho(0.0d0,0.0d0,-0.31d0,sigma(3:4,3:4))
OPEN(16,FILE="wf_0.0_0.0_-0.31_secular.txt")
DO i = 1, 100000
  CALL rk4_prop_coherence_kill(dt,sigma,Hls,t_len_secular,operators_secular)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4)))
END DO
CLOSE(16)


sigma = 0.0d0
CALL xyz_bloch_to_rho(0.2d0,-0.2d0,-0.31d0,sigma(3:4,3:4))
OPEN(16,FILE="wf_0.2_-0.2_-0.31_secular.txt")
DO i = 1, 100000
  CALL rk4_prop_coherence_kill(dt,sigma,Hls,t_len_secular,operators_secular)
  WRITE(16,*) i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4)))
END DO
CLOSE(16)


END PROGRAM
