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


MODULE qns

USE parameters
USE rk_ns_utils

IMPLICIT NONE

LOGICAL, PARAMETER :: hard_separate_bottom = .TRUE.  ! Hard code bottom two estates in separate, singly occupied secular blocks

! Lindblad operator
TYPE lb_op
  INTEGER, DIMENSION(nbath) :: occupation  ! How many have been added while building this thing
  COMPLEX*16, DIMENSION(n,n,nbath) :: A  
  COMPLEX*16, DIMENSION(n,n,nbath) :: Abasis
  REAL*8 :: omega
  COMPLEX*16, DIMENSION(nbath) :: thpl  ! Theta plus element
  REAL*8, DIMENSION(nbath) :: rate  ! Full jump rate
  INTEGER :: sblock, dblock  ! Source block number, destination block number
END TYPE

! Blocked rho structure.
TYPE blocked_rho
  INTEGER :: num_blocks  ! Number of secular blocks
  INTEGER, DIMENSION(:), ALLOCATABLE :: block_sizes  ! Number of entries in each block (num_blocks)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: block_limits  ! Start and end of each secular block (num_blocks,2)
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: rho  ! Nothing fancy with this yet, but could be done eventually (n,n)
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: committor_basis  !  Committors solved for each basis dmat entry
END TYPE


CONTAINS


! Takes an r, theta, phi location on the Bloch sphere
! phi is in the xy plane and theta is the angle to z
SUBROUTINE rtp_bloch_to_rho(r,t,p,rho)
  REAL*8, INTENT(IN) :: r, t, p  ! Bloch sphere coordinate
  COMPLEX*16, DIMENSION(2,2), INTENT(OUT) :: rho  ! density matrix

  CALL xyz_bloch_to_rho(r*SIN(t)*COS(p),r*SIN(t)*SIN(p),r*COS(t),rho)
END SUBROUTINE


! Takes an x, y, z location on the Bloch sphere
SUBROUTINE xyz_bloch_to_rho(x,y,z,rho)
  REAL*8, INTENT(IN) :: x, y, z  ! Bloch sphere coordinate
  COMPLEX*16, DIMENSION(2,2), INTENT(OUT) :: rho  ! density matrix

  rho = 0.0d0
  rho(1,1) = 1.0d0 + z
  rho(2,2) = 1.0d0 - z
  rho(1,2) = x - CMPLX(0.0d0,y)
  rho(2,1) = CONJG(rho(1,2))
  rho = rho/2.0d0
END SUBROUTINE


! Take a rho location and transform it to x, y, z coordinates
SUBROUTINE bloch_to_xyz(x,y,z,rho)
  REAL*8, INTENT(OUT) :: x, y, z
  COMPLEX*16, DIMENSION(2,2), INTENT(IN) :: rho

  x = 2.0d0*REAL(REAL(rho(1,2)))
  y = 2.0d0*AIMAG(rho(2,1))
  z = rho(1,1) - rho(2,2)
END SUBROUTINE


! Determine which block flux has left and which block it has arrived in due to
! the action of a single operator and add it to the flux tally
! Operator's action is passed in in terms of a single step's derivative
! This will break if more than one block is the ORIGIN of flux
! but can handle more than one destination
SUBROUTINE tally_flux_by_operator(n,brho,deriv_temp,net_fluxes)
  INTEGER :: n  ! deriv_temp size
  TYPE(blocked_rho), INTENT(IN) :: brho ! rho and blocks information
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: deriv_temp  ! Derivative resulting from operator
  REAL*8, DIMENSION(brho%num_blocks,brho%num_blocks), INTENT(INOUT) :: net_fluxes
  INTEGER :: i, j  ! Loop indices/temporaries
  REAL*8, DIMENSION(brho%num_blocks) :: positives  ! Calculate the positive population changes
  REAL*8, DIMENSION(brho%num_blocks) :: negatives  ! Calculate the negative population changes

  positives = 0.0d0
  negatives = 0.0d0
  DO i = 1, n
    IF (ABS(deriv_temp(i,i)) .GT. 1d-12)  THEN  ! Generous zero
      j = identify_block(i,brho)  ! Block id
      !WRITE(*,*) "block is ", j
      IF (REAL(REAL(deriv_temp(i,i))) .GT. 0.0d0) THEN  ! Add to proper tally pos/neg
        positives(j) = positives(j) + REAL(REAL(deriv_temp(i,i)))
      ELSE
        negatives(j) = negatives(j) + REAL(REAL(deriv_temp(i,i)))
      END IF 
    END IF
  END DO

  ! Cancel out interblock transfers to avoid double counts
  DO i = 1, brho%num_blocks
    IF (ABS(negatives(i)) .GT. ABS(positives(i))) THEN
      negatives(i) = negatives(i) + positives(i)
      positives(i) = 0.0d0
    ELSE
      positives(i) = positives(i) + negatives(i)
      negatives(i) = 0.0d0
    END IF
  END DO


  DO i = 1, brho%num_blocks  ! Find the origin of the loss of flux
    IF (ABS(negatives(i)) .GT. 1d-10) THEN
      j = i
      EXIT  ! Done. Skip rest of the loop
    END IF
  END DO 

  DO i = 1, brho%num_blocks  ! Flux went from there to here
    IF (positives(i) .GT. 1d-10) THEN
      net_fluxes(i,j) = positives(i) + net_fluxes(i,j)
      IF (i .EQ. j) THEN
        WRITE(*,*) "Issue encountered: i == j for "
        CALL print_mat(deriv_temp)
        WRITE(*,*) "With positives "
        WRITE(*,*) positives
        WRITE(*,*) "And negatives "
        WRITE(*,*) negatives
      END IF
    END IF
  END DO
END SUBROUTINE
  
! Tell me which block this eigenstate 'ind' index is in
INTEGER FUNCTION identify_block(ind,brho)
  INTEGER, INTENT(IN) :: ind  ! Eigenstate number
  TYPE(blocked_rho), INTENT(IN) :: brho
  INTEGER :: i

  identify_block = brho%num_blocks
  i = 1
  DO WHILE (i .LT. brho%num_blocks)
    !WRITE(*,*) "checking block ", i, "index", ind, "bounds", brho%block_limits(i,2)
    IF (ind .LE. brho%block_limits(i,2)) THEN
      identify_block = i
      EXIT
    END IF
    i = i + 1
  END DO
END FUNCTION


! Take the operators and calculate the lindblad dissipator plus
! Hls propagator to find drho/dt for the system matrix
! Idiot verison; no fancy stuff, just straight up matrix multiplication
! all the way through
SUBROUTINE calculate_derivative(deriv,rho,Hls,t_len,operators,brho,net_fluxes)
  INTEGER, INTENT(IN) :: t_len
  TYPE(lb_op), INTENT(IN), DIMENSION(t_len) :: operators
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: deriv
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: rho, Hls
  COMPLEX*16, DIMENSION(n,n) :: op_temp, op_conj, temp
  TYPE(blocked_rho), INTENT(IN), OPTIONAL :: brho  ! Blocked rho input, optional
  REAL*8, DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: net_fluxes  ! Counting interblock fluxes; size nbxnb

  COMPLEX*16, DIMENSION(n,n) :: deriv_temp
  INTEGER :: i, j

  deriv = 0.0d0
  DO i = 1, t_len  ! Dissipator terms
    DO j = 1, nbath
      op_temp = operators(i)%A(1:n,1:n,j)
      op_conj = CONJG(TRANSPOSE(op_temp))
      temp = MATMUL(op_conj,op_temp)
      deriv_temp =  operators(i)%rate(j)*(MATMUL(op_temp,MATMUL(rho,op_conj)) - &
        0.5d0*(MATMUL(temp,rho) + MATMUL(rho,temp)))

      IF (PRESENT(brho) .AND. PRESENT(net_fluxes)) THEN  ! Calculate the flux between block indexes here
        CALL tally_flux_by_operator(n,brho,deriv_temp,net_fluxes)
      END IF

      deriv = deriv + operators(i)%rate(j)*(MATMUL(op_temp,MATMUL(rho,op_conj)) - &
        0.5d0*(MATMUL(temp,rho) + MATMUL(rho,temp)))
    END DO
  END DO 

  deriv = deriv - CMPLX(0.0d0,1.0d0)*(MATMUL(Hls,rho) - MATMUL(rho,Hls))  ! Lambshift commutator terms
END SUBROUTINE

! Call the calculate derivative routine in order to move
! rho forward in time via the Lindblad equations
! Removes coherences from every rk4 step for use in 
! consistency checking with secular equations of motion
SUBROUTINE rk4_prop_coherence_kill(dt,rho,Hls,t_len,operators,derivative)
  REAL*8, INTENT(IN) :: dt
  INTEGER, INTENT(IN) :: t_len
  TYPE(lb_op), INTENT(IN), DIMENSION(t_len) :: operators
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: rho
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: Hls
  COMPLEX*16, DIMENSION(n,n) :: deriv, temp, k1, k2, k3, k4
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT), OPTIONAL :: derivative
  INTEGER :: i
  COMPLEX*16, DIMENSION(n,n) :: temp_mat, rho_temp

  rho_temp = rho
  rho = 0.0d0
  DO i = 1, n  ! Cut out all initial coherence that may exist
    rho(i,i) = rho_temp(i,i)
  END DO


  CALL calculate_derivative(deriv,rho,Hls,t_len,operators) 
  temp_mat = deriv
  deriv = 0.0d0
  DO i = 1, n
    deriv(i,i) = temp_mat(i,i)
  END DO
  k1 = dt*deriv
  temp = 0.5d0*k1 + rho
  CALL calculate_derivative(deriv,temp,Hls,t_len,operators) 
  temp_mat = deriv
  deriv = 0.0d0
  DO i = 1, n
    deriv(i,i) = temp_mat(i,i)
  END DO
  k2 = dt*deriv
  temp = 0.5d0*k2 + rho
  CALL calculate_derivative(deriv,temp,Hls,t_len,operators) 
  temp_mat = deriv
  deriv = 0.0d0
  DO i = 1, n
    deriv(i,i) = temp_mat(i,i)
  END DO
  k3 = dt*deriv 
  temp = k3 + rho
  CALL calculate_derivative(deriv,temp,Hls,t_len,operators) 
  temp_mat = deriv
  deriv = 0.0d0
  DO i = 1, n
    deriv(i,i) = temp_mat(i,i)
  END DO
  k4 = deriv*dt

  rho = rho_temp
  rho = rho + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4) 

  IF (PRESENT(derivative)) THEN
    derivative = deriv
  END IF
END SUBROUTINE


! Call the calculate derivative routine in order to move
! rho forward in time via the Lindblad equations
SUBROUTINE rk4_prop(dt,rho,Hls,t_len,operators,brho,net_fluxes)
  REAL*8, INTENT(IN) :: dt
  INTEGER, INTENT(IN) :: t_len
  TYPE(lb_op), INTENT(IN), DIMENSION(t_len) :: operators
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: rho
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: Hls
  COMPLEX*16, DIMENSION(n,n) :: deriv, temp, k1, k2, k3, k4
  TYPE(blocked_rho), INTENT(IN), OPTIONAL :: brho  ! Blocked rho input, optional
  REAL*8, DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: net_fluxes  ! Counting interblock fluxes; size nbxnb
  INTEGER :: i, j  ! For checking validity of a density matrix
  REAL*8 :: trace
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: nf1, nf2, nf3, nf4


  IF (PRESENT(net_fluxes) .AND. PRESENT(brho)) THEN
    ALLOCATE(nf1(brho%num_blocks,brho%num_blocks))
    ALLOCATE(nf2(brho%num_blocks,brho%num_blocks))
    ALLOCATE(nf3(brho%num_blocks,brho%num_blocks))
    ALLOCATE(nf4(brho%num_blocks,brho%num_blocks))
    nf1 = 0.0d0
    nf2 = 0.0d0
    nf3 = 0.0d0
    nf4 = 0.0d0
    CALL calculate_derivative(deriv,rho,Hls,t_len,operators,brho,nf1)
    k1 = dt*deriv
    temp = 0.5d0*k1 + rho
    CALL calculate_derivative(deriv,temp,Hls,t_len,operators,brho,nf2)
    k2 = dt*deriv
    temp = 0.5d0*k2 + rho
    CALL calculate_derivative(deriv,temp,Hls,t_len,operators,brho,nf3)
    k3 = dt*deriv 
    temp = k3 + rho
    CALL calculate_derivative(deriv,temp,Hls,t_len,operators,brho,nf4)
    k4 = deriv*dt
   
    net_fluxes = net_fluxes + dt*(1.0d0/6.0d0)*(nf1 + 2.0d0*nf2 + 2.0d0*nf3 + nf4) 
    DEALLOCATE(nf1)
    DEALLOCATE(nf2)
    DEALLOCATE(nf3)
    DEALLOCATE(nf4)
  ELSE
    CALL calculate_derivative(deriv,rho,Hls,t_len,operators) 
    k1 = dt*deriv
    temp = 0.5d0*k1 + rho
    CALL calculate_derivative(deriv,temp,Hls,t_len,operators) 
    k2 = dt*deriv
    temp = 0.5d0*k2 + rho
    CALL calculate_derivative(deriv,temp,Hls,t_len,operators) 
    k3 = dt*deriv 
    temp = k3 + rho
    CALL calculate_derivative(deriv,temp,Hls,t_len,operators) 
    k4 = deriv*dt
  END IF

  rho = rho + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4) 


  trace = 0.0d0
  i = 0  ! I may want to turn this back on
  j = 0

END SUBROUTINE

! Break up the eigenvalue number line into secular blocks;
! evals should already be sorted
SUBROUTINE identify_secular_blocks(lambda,evals,blocknum,blocks)
  INTEGER, INTENT(OUT) :: blocknum  ! How many blocks there be
  INTEGER, DIMENSION(n,2), INTENT(OUT) :: blocks  ! Each entry is start and end of a block
  REAL*8, INTENT(IN), DIMENSION(n) :: evals
  REAL*8, INTENT(IN) ::  lambda  ! Energy difference to define a secular block
  INTEGER :: i

  IF (n .EQ. 2) THEN  ! SB test system may turn up
    IF (evals(2) - evals(1) .GT. lambda) THEN
      blocknum = 2
      blocks(1,1) = 1
      blocks(1,2) = 1
      blocks(2,1) = 2
      blocks(2,2) = 2
    ELSE
      blocknum = 1
      blocks(1,1) = 1
      blocks(1,2) = 2
    END IF
  ELSE

    IF (hard_separate_bottom .EQV. .FALSE.) THEN
      i = 1
      blocknum = 1
      blocks(1,1) = 1
      blocks(1,2) = 1  ! May be overwritten
      DO i = 2, n
        IF (evals(i) - evals(i-1) .GT. lambda) THEN  ! Moved outside the secular block
          blocks(blocknum,2) = i - 1   ! Put in end value
          blocknum = blocknum + 1
          blocks(blocknum,1) = i  ! Put in start value
        END IF 
      END DO
      blocks(blocknum,2) = n
    ELSE

      blocknum = 3  ! Hard code in the bottom two states as firmly separate
      blocks(1,1) = 1
      blocks(1,2) = 1
      blocks(2,1) = 2
      blocks(2,2) = 2
      blocks(3,1) = 3
      blocks(3,2) = 3  ! May be overwritten
      DO i = 4, n
        IF (evals(i) - evals(i-1) .GT. lambda) THEN  ! Moved outside the secular block
          blocks(blocknum,2) = i - 1   ! Put in end value
          blocknum = blocknum + 1
          blocks(blocknum,1) = i  ! Put in start value
        END IF 
      END DO
      blocks(blocknum,2) = n

    END IF
  END IF
END SUBROUTINE


! Given the secular blocks have been defined, find the jump operators between them
SUBROUTINE find_sb_jops(lambda,evals,operators,t_len,couple_op,brho)
  REAL*8, INTENT(IN) :: lambda  ! splitting for secular blocks
  REAL*8, INTENT(IN), DIMENSION(n) :: evals
  INTEGER, INTENT(OUT) :: t_len
  TYPE(blocked_rho), INTENT(INOUT) :: brho
  TYPE(lb_op), INTENT(OUT), DIMENSION(:), ALLOCATABLE :: operators
  INTEGER :: blocknum
  INTEGER, DIMENSION(n,2) :: blocks
  INTEGER :: i, j, k, l, bn, oind
  REAL*8, DIMENSION(:), ALLOCATABLE :: avg_e
  COMPLEX*16, INTENT(IN), DIMENSION(n,n,nbath) :: couple_op

  CALL identify_secular_blocks(lambda,evals,blocknum,blocks) 
  ALLOCATE(avg_e(blocknum))
  WRITE(*,*) "Nonecular blocks have the following bounds:"
  DO i = 1, blocknum
    WRITE(*,*) blocks(i,1), blocks(i,2)
  END DO
 

  DO i = 1, blocknum  ! Calculate the average energy of each block
    j = blocks(i,1)  ! Start entry
    avg_e(i) = 0.0d0
    DO WHILE (j .LE. blocks(i,2))  ! Until end of the block
      avg_e(i) = avg_e(i) + evals(j)   ! Add for every block entry
      j = j + 1
    END DO
    avg_e(i) = avg_e(i) / (blocks(i,2) - blocks(i,1) + 1)  ! Average over entries
  END DO

  CALL setup_rho_blocks(brho,blocknum,blocks)

  t_len = blocknum**2 - blocknum + 1  ! Full number of jump operators between blocks plus dephasing
  ALLOCATE(operators(t_len))
  oind = 1

  ! For each block defined
  DO i = 1, blocknum
    ! Define an operator jumping into that block from every other block
    DO j = 1, blocknum

      IF (i .EQ. j) THEN  ! No interblock jump operators; deal with dephasing elsewhere
        CYCLE
      END IF
 
      operators(oind)%occupation = 0
      operators(oind)%A = 0.0d0
      operators(oind)%Abasis = 0.0d0
      operators(oind)%omega = avg_e(i) - avg_e(j)  ! Is not backwards.
      operators(oind)%sblock = j
      operators(oind)%dblock = i
      DO bn = 1, nbath
        operators(oind)%thpl(bn) = get_thpl_element(operators(oind)%omega,bn,0.0d0,0.0d0)
        operators(oind)%rate(bn) = 2.0d0*REAL(REAL(operators(oind)%thpl(bn)))  ! Does need a 2.0 
      END DO

     
      DO k = blocks(i,1), blocks(i,2)  ! For all entries in the destination block
        DO l = blocks(j,1), blocks(j,2)  ! For all entries in the source block
          DO bn = 1, nbath
            IF (ABS(couple_op(k,l,bn)) .GT. 1d-10) THEN  ! Legitimate entry
              operators(oind)%occupation(bn) = operators(oind)%occupation(bn) + 1
              operators(oind)%A(k,l,bn) = couple_op(k,l,bn)  ! Select the coupling operator entry for k,l,bn 
            END IF
          END DO 
        END DO
      END DO

      oind = oind + 1  ! Next operator index 
    END DO
  END DO 

  IF (oind .NE. t_len) THEN
    WRITE(*,*) "WARNING: You've messed up the indexing of your operators"
  END IF

  ! Dephasing operator 
  operators(oind)%A = 0.0d0
  operators(oind)%omega = 0.0d0
  operators(oind)%occupation = 0
  operators(oind)%sblock = 0  ! All blocks
  operators(oind)%dblock = 0  ! All blocks
  DO bn = 1, nbath
    operators(oind)%thpl(bn) = get_thpl_element(operators(oind)%omega,bn,0.0d0,0.0d0)
    operators(oind)%rate(bn) = 2.0d0*REAL(REAL(operators(oind)%thpl(bn)))  ! Yes. Needs a 2
  END DO
  DO i = 1, blocknum
    DO k = blocks(i,1), blocks(i,2)
      DO j = blocks(i,1), blocks(i,2)
        DO bn = 1, nbath
          IF (ABS(couple_op(k,j,bn)) .GT. 1d-10) THEN  ! Legitimate entry
            operators(oind)%occupation(bn) = operators(oind)%occupation(bn) + 1
            WRITE(*,*) "occupation is now ", operators(oind)%occupation(bn), "from", couple_op(k,j,bn)
            operators(oind)%A(k,j,bn) = couple_op(k,j,bn) ! Is not necessarily diagonal
          END IF
        END DO
      END DO
    END DO
  END DO

  DEALLOCATE(avg_e)
END SUBROUTINE

! Find the lambshift Hamiltonian
SUBROUTINE assemble_Hls(t_len,operators,Hls,H)
  INTEGER, INTENT(IN) :: t_len
  TYPE(lb_op), DIMENSION(t_len), INTENT(IN) :: operators
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: Hls
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: H
  INTEGER :: i, j
  COMPLEX*16, DIMENSION(n,n) :: m_temp

  Hls = H
  DO i = 1, t_len
    DO j = 1, nbath
      m_temp = operators(i)%A(1:n,1:n,j)
      Hls = Hls + AIMAG(operators(i)%thpl(j))*MATMUL(CONJG(TRANSPOSE(m_temp)),m_temp)
    END DO
  END DO

END SUBROUTINE

! Set up a blocked rho structure; must not 
! be allocated already in any way
SUBROUTINE setup_rho_blocks(brho,nblocks,blocks)
  TYPE(blocked_rho), INTENT(INOUT) :: brho
  INTEGER, INTENT(IN) :: nblocks
  INTEGER, INTENT(IN), DIMENSION(n,2) :: blocks
  INTEGER :: i
 
  brho%num_blocks = nblocks
  ALLOCATE(brho%rho(n,n))
  ALLOCATE(brho%block_sizes(nblocks))
  ALLOCATE(brho%block_limits(nblocks,2))

  DO i = 1, nblocks  ! Fill in block limites and sizes
    brho%block_limits(i,1:2) = blocks(i,1:2) 
    brho%block_sizes(i) = blocks(i,2) - blocks(i,1) + 1
    WRITE(*,*) "brho ", i, "block_sizes", brho%block_sizes(i), "limits", brho%block_limits(i,1:2)
  END DO
END SUBROUTINE


! Helper function print out the blocks of a blocked rho matrix
! must be allocated
SUBROUTINE print_blocks(brho,b)
  TYPE(blocked_rho), INTENT(IN) :: brho
  INTEGER, INTENT(IN), OPTIONAL :: b  ! Optionally only print one block
  INTEGER :: i, l1, l2

  IF (PRESENT(b)) THEN
    WRITE(*,*) "printing block ", b, " of ", brho%num_blocks
    l1 = brho%block_limits(b,1)
    l2 = brho%block_limits(b,2)
    CALL print_mat(brho%rho(l1:l2,l1:l2))
  ELSE
    DO i = 1, brho%num_blocks
      l1 = brho%block_limits(i,1)
      l2 = brho%block_limits(i,2)
      WRITE(*,*) "printing block ", i, " of ", brho%num_blocks
      CALL print_mat(brho%rho(l1:l2,l1:l2))
    END DO
  END IF
END SUBROUTINE

! Divide to set the trace to 1; assumed to have real diagonal entries
SUBROUTINE normalize_trace(d,rho)
  INTEGER, INTENT(IN) :: d
  COMPLEX*16, DIMENSION(d,d), INTENT(INOUT) :: rho
  REAL*8 :: dsum
  INTEGER :: i

  dsum = 0.0d0
  DO i = 1, d
    dsum = dsum + REAL(REAL(rho(i,i)))
  END DO
  WRITE(*,*) "DSUM", dsum
  rho = rho/dsum
END SUBROUTINE


! Run with absorbing boundary conditions to specified states
! until a certain amount of density is localized in the requested
! states. Return the committor vector; zeroes out coherences on every 
! step
! The destination states really should be in separate secular blocks
SUBROUTINE abc_committor_noco(dt,d,sigma,Hls,t_len,operators,numstates,state_ids,states,threshold)
  REAL*8, INTENT(IN) :: dt, threshold  ! Time step, when to conclude propogation complete
  INTEGER, INTENT(IN) :: d, t_len  ! Dimension of the matrix, number of operators
  INTEGER, INTENT(IN) :: numstates ! Number of final, concluding states
  TYPE(lb_op), DIMENSION(t_len), INTENT(IN) :: operators  ! LB jump operators
  REAL*8, DIMENSION(numstates), INTENT(OUT) :: states   ! Density that ends in state 1 and 2 respectively
  INTEGER, DIMENSION(numstates), INTENT(IN) :: state_ids  ! Which states to use
  COMPLEX*16, INTENT(IN), DIMENSION(d,d) :: sigma  ! Input density matrix
  COMPLEX*16, INTENT(IN), DIMENSION(d,d) :: Hls   ! Lamb shif hamiltonian
  COMPLEX*16, DIMENSION(d,d) :: sigma_survive
  REAL*8 :: survived
  INTEGER :: i, j, k

  sigma_survive = sigma
  states = 0.0d0
  DO  ! While density remains to evolve

    DO j = 1, d  ! Zero out all coherences before evolving a timestep
      DO k = 1, d
        IF (j .NE. k) THEN
          sigma_survive(j,k) = 0.0d0
        END IF
      END DO
    END DO
 
    CALL rk4_prop(dt,sigma_survive,Hls,t_len,operators)

    DO i = 1, numstates  ! Copy in the new committor contribution; zero out in evolving matrix
      states(i) = states(i) + sigma_survive(state_ids(i),state_ids(i))
    
      DO j = 1, d
        sigma_survive(state_ids(i),j) = 0.0d0
        sigma_survive(j,state_ids(i)) = 0.0d0
      END DO
   END DO

    survived = 0.0d0
    DO i = 1, n  ! Count up surviving probability
      DO j = 1, n
        survived = survived + ABS(sigma_survive(i,j))
      END DO
    END DO

    IF (survived .LT. threshold) THEN
      EXIT  ! Decayed to the end
    END IF
  END DO

END SUBROUTINE



! Run with absorbing boundary conditions to specified states
! until a certain amount of density is localized in the requested
! states. Return the committor vector
! The destination states really should be in separate secular blocks
SUBROUTINE abc_committor(dt,d,sigma,Hls,t_len,operators,numstates,state_ids,states,threshold)
  REAL*8, INTENT(IN) :: dt, threshold  ! Time step, when to conclude propogation complete
  INTEGER, INTENT(IN) :: d, t_len  ! Dimension of the matrix, number of operators
  INTEGER, INTENT(IN) :: numstates ! Number of final, concluding states
  TYPE(lb_op), DIMENSION(t_len), INTENT(IN) :: operators  ! LB jump operators
  REAL*8, DIMENSION(numstates), INTENT(OUT) :: states   ! Density that ends in state 1 and 2 respectively
  INTEGER, DIMENSION(numstates), INTENT(IN) :: state_ids  ! Which states to use
  COMPLEX*16, INTENT(IN), DIMENSION(d,d) :: sigma  ! Input density matrix
  COMPLEX*16, INTENT(IN), DIMENSION(d,d) :: Hls   ! Lamb shif hamiltonian
  COMPLEX*16, DIMENSION(d,d) :: sigma_survive
  REAL*8 :: survived
  INTEGER :: i, j

  sigma_survive = sigma
  states = 0.0d0
  DO
    CALL rk4_prop(dt,sigma_survive,Hls,t_len,operators)

    DO i = 1, numstates  ! Copy in the new committor contribution; zero out in evolving matrix
      states(i) = states(i) + sigma_survive(state_ids(i),state_ids(i))
    
      DO j = 1, d
        sigma_survive(state_ids(i),j) = 0.0d0
        sigma_survive(j,state_ids(i)) = 0.0d0
      END DO
    END DO

    survived = 0.0d0
    DO i = 1, n  ! Count up surviving probability
      DO j = 1, n
        survived = survived + ABS(sigma_survive(i,j))
      END DO
    END DO

    IF (survived .LT. threshold) THEN
      EXIT  ! Decayed to the end
    END IF
  END DO

END SUBROUTINE

! Run abc_committor to solve for all the 'basis states' necessary to describe a
! density matrix
! The solved basis is a matrix of committors; the lower triangle is real parts of the
! coherence and the upper triangle is the imaginary parts of the coherence with the
! positive of the complex conjugate pair on the lower triangle
! You can select out a subset of the basis to assemble using bound1/bound2 as numbers other than
! 1 and n
SUBROUTINE solve_basis(dt,d,threshold,Hls,t_len,operators,numstates,state_ids,solved_basis,bound1,bound2)
  REAL*8, INTENT(IN) :: dt, threshold  ! Time step, when to conclude propogation complete
  INTEGER, INTENT(IN) :: d, t_len  ! Dimension of the matrix, number of operators
  INTEGER, INTENT(IN) :: numstates ! Number of final, concluding states
  INTEGER, INTENT(IN) :: bound1, bound2  ! Subset of matrix basis to run
  TYPE(lb_op), DIMENSION(t_len), INTENT(IN) :: operators  ! LB jump operators
  INTEGER, DIMENSION(numstates), INTENT(IN) :: state_ids  ! Which states to use
  REAL*8, INTENT(OUT), DIMENSION(d,d,numstates) :: solved_basis 
  COMPLEX*16, INTENT(IN), DIMENSION(d,d) :: Hls   ! Lamb shif hamiltonian
  COMPLEX*16, DIMENSION(d,d) :: sigma_temp
  REAL*8, DIMENSION(numstates) :: states

  LOGICAL :: is_end
  INTEGER :: i, j, k
  
  solved_basis = 0.0d0

  DO i = 1, numstates  ! Committor of 100% for the destination state we are currently in
    solved_basis(state_ids(i),state_ids(i),i) = 1.0d0
  END DO

  DO i = bound1, bound2  ! row
    is_end = .FALSE.
    DO k = 1, numstates
      IF (i .EQ. state_ids(k)) THEN  ! Don't run on an end state
        is_end = .TRUE.
      END IF
    END DO
    IF (is_end .EQV. .TRUE.) THEN  ! Skip this part of the loop
      CYCLE
    END IF    


    DO j = bound1, bound2  ! column
      WRITE(*,*) "Solving ", i, j
      is_end = .FALSE.
      DO k = 1, numstates
        IF (j .EQ. state_ids(k)) THEN  ! Don't run on an end state
          is_end = .TRUE.
        END IF
      END DO
      IF (is_end .EQV. .TRUE.) THEN  ! Skip this part of the loop
        CYCLE
      END IF    


      IF (i .EQ. j) THEN
        states = 0.0d0
        sigma_temp = 0.0d0
        sigma_temp(i,i) = 1.0d0 
        CALL abc_committor(dt,d,sigma_temp,Hls,t_len,operators,numstates,state_ids,states,threshold)
        WRITE(*,*) "solved basis ", i, j, "is", states
        solved_basis(i,i,1:numstates) = states
      ELSE IF (i .LT. j) THEN  ! Upper triangle
        states = 0.0d0
        sigma_temp = 0.0d0
        sigma_temp(i,j) = CMPLX(0.0d0,-1.0d0)
        sigma_temp(j,i) = CMPLX(0.0d0,1.0d0)
        CALL abc_committor(dt,d,sigma_temp,Hls,t_len,operators,numstates,state_ids,states,threshold)
        WRITE(*,*) "solved basis ", i, j, "is", states
        solved_basis(i,j,1:numstates) = states
      ELSE  ! Lower triangle
        states = 0.0d0
        sigma_temp = 0.0d0
        sigma_temp(i,j) = 1.0d0
        sigma_temp(j,i) = 1.0d0
        CALL abc_committor(dt,d,sigma_temp,Hls,t_len,operators,numstates,state_ids,states,threshold)
        WRITE(*,*) "solved basis ", i, j, "is", states
        solved_basis(i,j,1:numstates) = states
      END IF
    END DO
  END DO 

END SUBROUTINE

! Use the calculated n^2 size basis for the density matrix to
! determine the committor value
SUBROUTINE calculate_committors_by_basis(d,sigma,solved_basis,numstates,comms)
  INTEGER, INTENT(IN) :: d, numstates  ! Dimension, number of final exit states
  COMPLEX*16, DIMENSION(d,d), INTENT(IN) :: sigma  ! The density matrix to calculate committors on
  REAL*8, DIMENSION(d,d,numstates), INTENT(IN) :: solved_basis  ! committors by basis entry 
  REAL*8, DIMENSION(numstates), INTENT(OUT) :: comms
  INTEGER :: i, j

  comms = 0.0d0
  DO i = 1, d  ! Scan over the upper triangle
    DO j = i, d
      IF (i .EQ. j) THEN  ! Put together contributions from all populations and coherences
        comms = comms + REAL(REAL(sigma(i,j)))*solved_basis(i,j,1:numstates)
      ELSE
        comms = comms + REAL(REAL(sigma(i,j)))*solved_basis(j,i,1:numstates)
        comms = comms - AIMAG(sigma(i,j))*solved_basis(i,j,1:numstates)
      END IF
    END DO
  END DO
END SUBROUTINE
  


! Run abc_committor to solve for all the 'basis states' necessary to describe a
! density matrix; zero out coherences on every step
! The solved basis is a matrix of committors; the lower triangle is real parts of the
! coherence and the upper triangle is the imaginary parts of the coherence with the
! positive of the complex conjugate pair on the lower triangle
! You can select out a subset of the basis to assemble using bound1/bound2 as numbers other than
! 1 and n
SUBROUTINE solve_basis_noco(dt,d,threshold,Hls,t_len,operators,numstates,state_ids,solved_basis,bound1,bound2)
  REAL*8, INTENT(IN) :: dt, threshold  ! Time step, when to conclude propogation complete
  INTEGER, INTENT(IN) :: d, t_len  ! Dimension of the matrix, number of operators
  INTEGER, INTENT(IN) :: numstates ! Number of final, concluding states
  INTEGER, INTENT(IN) :: bound1, bound2  ! Subset of matrix basis to run
  TYPE(lb_op), DIMENSION(t_len), INTENT(IN) :: operators  ! LB jump operators
  INTEGER, DIMENSION(numstates), INTENT(IN) :: state_ids  ! Which states to use
  REAL*8, INTENT(OUT), DIMENSION(d,d,numstates) :: solved_basis 
  COMPLEX*16, INTENT(IN), DIMENSION(d,d) :: Hls   ! Lamb shif hamiltonian
  COMPLEX*16, DIMENSION(d,d) :: sigma_temp
  REAL*8, DIMENSION(numstates) :: states

  LOGICAL :: is_end
  INTEGER :: i, j, k
  
  solved_basis = 0.0d0

  DO i = 1, numstates  ! Committor of 100% for the destination state we are currently in
    solved_basis(state_ids(i),state_ids(i),i) = 1.0d0
  END DO

  DO i = bound1, bound2  ! row
    is_end = .FALSE.
    DO k = 1, numstates
      IF (i .EQ. state_ids(k)) THEN  ! Don't run on an end state
        is_end = .TRUE.
      END IF
    END DO
    IF (is_end .EQV. .TRUE.) THEN  ! Skip this part of the loop
      CYCLE
    END IF    


    DO j = bound1, bound2  ! column
      WRITE(*,*) "Solving ", i, j
      is_end = .FALSE.
      DO k = 1, numstates
        IF (j .EQ. state_ids(k)) THEN  ! Don't run on an end state
          is_end = .TRUE.
        END IF
      END DO
      IF (is_end .EQV. .TRUE.) THEN  ! Skip this part of the loop
        CYCLE
      END IF    


      IF (i .EQ. j) THEN
        states = 0.0d0
        sigma_temp = 0.0d0
        sigma_temp(i,i) = 1.0d0 
        CALL abc_committor_noco(dt,d,sigma_temp,Hls,t_len,operators,numstates,state_ids,states,threshold)
        WRITE(*,*) "solved basis ", i, j, "is", states
        solved_basis(i,i,1:numstates) = states
      ELSE IF (i .LT. j) THEN  ! Upper triangle
        states = 0.0d0
        sigma_temp = 0.0d0
        sigma_temp(i,j) = CMPLX(0.0d0,-1.0d0)
        sigma_temp(j,i) = CMPLX(0.0d0,1.0d0)
        CALL abc_committor_noco(dt,d,sigma_temp,Hls,t_len,operators,numstates,state_ids,states,threshold)
        WRITE(*,*) "solved basis ", i, j, "is", states
        solved_basis(i,j,1:numstates) = states
      ELSE  ! Lower triangle
        states = 0.0d0
        sigma_temp = 0.0d0
        sigma_temp(i,j) = 1.0d0
        sigma_temp(j,i) = 1.0d0
        CALL abc_committor_noco(dt,d,sigma_temp,Hls,t_len,operators,numstates,state_ids,states,threshold)
        WRITE(*,*) "solved basis ", i, j, "is", states
        solved_basis(i,j,1:numstates) = states
      END IF
    END DO
  END DO 

END SUBROUTINE


END MODULE

