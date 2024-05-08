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


! Integration functions and machinery for Lindblad operators.

! Theta, G and Theta+/Theta- notation in this file follows the derivation in 
! Pollard, W. T. and Frisner, R. A. "Solutions of the Redfield equation for the dissipative quantum dynamics of 
! multilevel systems", J. Chem. Phys. 100, 1994.

MODULE set_up_H
USE parameters
USE prequel
USE rk_ns_utils

IMPLICIT NONE

CONTAINS 


! Setup a DVR sinc basis kinetic energy matrix, 1d
! x_i = i \dx*i
SUBROUTINE sinc_ke(dx,dim_coor,mass,kemat)
  REAL*8, INTENT(IN) :: dx  ! delta x interval
  REAL*8, INTENT(IN) :: mass ! mass
  COMPLEX*16, DIMENSION(dim_coor,dim_coor), INTENT(OUT) :: kemat 
  INTEGER, INTENT(IN) :: dim_coor

  INTEGER :: i, j, lim
  REAL*8 :: temp1, temp2
  temp1 = (hbar**2.0d0)/(2.0d0*mass*(dx**2.0d0))
  temp2 = (pi**2.0d0)/3.0d0

  lim = FLOOR(dim_coor/2.0d0)
  kemat = 0.0d0
  WRITE(*,*) "lim is ", lim
  DO i = -lim, lim
    DO j = -lim, lim
      IF (i .EQ. j) THEN
        kemat(i+lim+1,j+lim+1) = temp1*temp2*((-1.0d0)**(i-j))
      ELSE
        kemat(i+lim+1,j+lim+1) = temp1*((-1.0d0)**(i-j))*2.0d0/((i-j)**2.0d0)
      END IF
    END DO 
  END DO 

END SUBROUTINE



! Helper function for the DVR sinc basis; just returns the value
! of the full potential at the particular point in question, i*dx = R
! This is the part that will be tensor producted with the identity
! cob**4/(16 ceb)*R^4 - 0.5*cob*R**22 - c*R**3 + omegas_c*(\eta_c^2*mu(R))^2
COMPLEX*16 FUNCTION sinc_V(i,dx)
  INTEGER, INTENT(IN) :: i
  REAL*8, INTENT(IN) :: dx
  REAL*8 :: R, R2, R3, R4, tanhR  ! position and position squared

  R = dx*i
  R2 = R**2.0d0
  R3 = R*R2
  R4 = R3*R
  tanhR = coeff_v*TANH(coeff_y*R) + coeff_z*R
  sinc_V = (cob**4.0d0)*R4/(16.0d0*ceb) - 0.5d0*(cob**2.0d0)*R2 - coeff_c*R3 + omegasc*(etac*tanhR)**2.0d0
END FUNCTION

! Helper function for the DVR sinc basis; returns the value
! of the dipole at the particular point in question
! This is the part that will be tensor producted with q_c
! Includes prefacotrs!
! (omegasc^2)*SQRT(2.0d0/omegasc)*eta_c*mu
COMPLEX*16 FUNCTION sinc_V_cross(i,dx)
  INTEGER, INTENT(IN) :: i
  REAL*8, INTENT(IN) :: dx
  REAL*8 :: loc

  loc = dx*i
  sinc_V_cross = coeff_z*loc + coeff_v*TANH(coeff_y*loc)
  sinc_V_cross = sinc_V_cross*SQRT(2.0d0/omegasc)*etac*(omegasc**2.0d0)
END FUNCTION

! Output information about a Shin Menitu model
! gives position and energy because diabatic states are not a thing
SUBROUTINE output_info_SM_dmat(time,fd,dmat,H)
  INTEGER, INTENT(IN) :: fd
  COMPLEX*16, INTENT(IN), DIMENSION(n,n) :: dmat  ! Density matrix
  COMPLEX*16, INTENT(IN), DIMENSION(n,n) :: H
  REAL*8, INTENT(IN) :: time
  REAL*8 :: Q_x, Q_y, H_temp
  REAL*8, DIMENSION(2*n) :: output_wfn
  INTEGER :: i

  ! Usual wavefunction format for plotting via wfn_plot, even
  ! though it's a mixed state it should be fine
  DO i = 1, n
    output_wfn(2*i - 1) = dmat(i,i)
    output_wfn(2*i) = 0.0d0
  END DO

  Q_x = trace_mat(n,MATMUL(Q_x_g,dmat),1,n)
  Q_y = trace_mat(n,MATMUL(Q_y_g,dmat),1,n)
  H_temp = trace_mat(n,MATMUL(H,dmat),1,n)
  WRITE(fd,*) time, Q_x, Q_y, H_temp, output_wfn
END SUBROUTINE

! Output information about a conical intersection
! density matrix; gives the time, the diabatic popluations
! and the coupling/tuning positions
SUBROUTINE output_info_CI_dmat(time,fd,dmat,evec,evec_inv)
  INTEGER, INTENT(IN) :: fd
  COMPLEX*16, INTENT(IN), DIMENSION(n,n) :: dmat  ! Density matrix
  COMPLEX*16, INTENT(IN), DIMENSION(m,n) :: evec  ! Energy eigenvectors
  COMPLEX*16, INTENT(IN), DIMENSION(n,m) :: evec_inv
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: dmat_temp
  REAL*8, INTENT(IN) :: time
  REAL*8 :: Q_x, Q_y, temp, tm1, tm2

  ALLOCATE(dmat_temp(m,m))
  Q_x = trace_mat(n,MATMUL(Q_x_g,dmat),1,n)
  Q_y = trace_mat(n,MATMUL(Q_y_g,dmat),1,n)
  CALL from_basis(dmat,dmat_temp,evec,evec_inv)
  tm1 = trace_mat(m,dmat_temp,1,m/2)
  tm2 = trace_mat(m,dmat_temp,m/2+1,m)
  temp = 1.0d0/(tm1 + tm2)
  WRITE(fd,*) time, trace_mat(m,dmat_temp,1,m/2), trace_mat(m,dmat_temp,m/2+1,m), Q_x, Q_y 
  DEALLOCATE(dmat_temp)
END SUBROUTINE


! Call the setup routines for specific types; when multiple different models were supported,
! this had a much bigger role
SUBROUTINE setup(H, sigma, evec, eval, H_trans, estate_probs,couple_op,Qt_full)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: H, sigma  ! Hamiltonian, starting denisty matrix
  COMPLEX*16, DIMENSION(m,n), INTENT(OUT) :: evec  ! Energy eigenvectors
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: couple_op  ! System Bath coupling operator
  COMPLEX*16, DIMENSION(m,m), INTENT(OUT), OPTIONAL :: Qt_full  ! Full basis position operator
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval  ! Energy eigenvalues
  ! Boltzman likelihood to start in oscillator eigenstate; used by thermal init and DEPRECATED
  REAL*8, INTENT(OUT), DIMENSION(n) :: estate_probs
  ! Whose translation to the unshifted oscillator basis
  ! looks like this; each vector is its own column; used by thermal init and DEPRECATED
  COMPLEX*16, INTENT(OUT), DIMENSION(n,n) :: H_trans

  IF (PRESENT(Qt_full) .AND. run_type .NE. CI) THEN
    WRITE(*,*) "Error: optional argument Qt_full not implemented for this system type (function call to setup)."
    STOP
  END IF 

  IF (run_type .EQ. CI) THEN
    IF (nbath .EQ. 2) THEN
      IF (PRESENT(Qt_full)) THEN
        CALL setup_H_CI_2bath(H,sigma,evec,eval,H_trans,estate_probs,couple_op,Qt_full)
      ELSE
        CALL setup_H_CI_2bath(H,sigma,evec,eval,H_trans,estate_probs,couple_op)
      END IF
    ELSE
      WRITE(*,*) "Error. One bath no longer supported. Use two baths; set second eta to zero."
      STOP
    END IF
  ELSE IF (run_type .EQ. SM) THEN 
    IF (PRESENT(Qt_full)) THEN
      CALL setup_H_shin_metiu(H,sigma,evec,eval,H_trans,estate_probs,couple_op,Qt_full)
    ELSE
      CALL setup_H_shin_metiu(H,sigma,evec,eval,H_trans,estate_probs,couple_op)
    END IF
  ELSE
    WRITE(*,*) "Unrecognized runtype"
    STOP
  END IF
END SUBROUTINE


! Set up a conical intersection Hamiltonian.
! Used to setup two baths, one on the coupling mode and 
! one on the tuning mode.
SUBROUTINE setup_H_CI_2bath(H, sigma, evec, eval, H_trans, estate_probs,couple_op,Qt_full)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: H, sigma  ! Hamiltonian, initial density matrix
  COMPLEX*16, DIMENSION(m,n), INTENT(OUT) :: evec  ! H eigenvectors
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: couple_op  ! Coupling operator for each bath
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval  ! Hamiltonian eigenvalues
  ! Boltzman likelihood to start in oscillator eigenstate (no longer used)
  REAL*8, INTENT(OUT), DIMENSION(n) :: estate_probs  ! DEPRECATED; first is always 1
  COMPLEX*16, OPTIONAL, INTENT(OUT), DIMENSION(m,m) :: Qt_full  ! An optional full Q operator
  ! Whose translation to the unshifted oscillator basis
  ! looks like this; each vector is its own column
  COMPLEX*16, INTENT(OUT), DIMENSION(n,n) :: H_trans  ! DEPRECATED Only the first vector is used now.
  ! Operators used during setup
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: Q_c, Q_t, H_temp, proj1m, proj2m
  ! For diagonalization of Q_c to get |Q_c| if needed
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: Qc_evec, Qc_sqr, Qc_abs, Q_temp
  COMPLEX*16, DIMENSION(m/2) :: Qc_eval
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: proj1n, proj2n
  ! Bounds and loop indexes to fill the matrices
  INTEGER :: k, k_p, nc, nt, nc_p, nt_p, index_1, index_2, i, bound1, bound2
  REAL*8 :: total
  TYPE(wfunc) :: wfn

  ! For building in the larger Hilbert space before reduction
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: H_m, sigma_m, evec_m
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: evec_inv
  COMPLEX*16, DIMENSION(m) :: evals_m
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: H_trans_m, tvecs
  REAL*8, DIMENSION(:), ALLOCATABLE :: tprobs

  REAL*8, DIMENSION(10000) :: bins  ! For density of states plots
  REAL*8 :: bin_highest, bin_lowest, bin_width
  INTEGER :: binind

  ! Set up a one dimensional position operator 
  ALLOCATE(Q_temp(m,m))
  ALLOCATE(tprobs(m/2))
  ALLOCATE(tvecs(m/2,m/2))
  ALLOCATE(H_m(m,m))
  ALLOCATE(sigma_m(m,m))
  ALLOCATE(evec_m(m,m))
  ALLOCATE(evec_inv(n,m))
  ALLOCATE(H_trans_m(m,m))
  ALLOCATE(proj1n(n,n))
  ALLOCATE(proj2n(n,n))
  ALLOCATE(Q_c(m,m))
  ALLOCATE(Q_t(m,m))
  ALLOCATE(H_temp(m,m))
  ALLOCATE(proj1m(m,m))
  ALLOCATE(proj2m(m,m))


  ! Setup a Q matrix
  DO nc = 1, m 
    DO nc_p = 1, m 
      Q_temp(nc,nc_p) = del(nc,nc_p+1)*SQRT((nc_p)/2.0d0) + &
        del(nc,nc_p-1)*SQRT((nc_p - 1.0)/2.0d0)
     END DO
  END DO

  ! Set up the |Q_c| matrix if absolute coupling was requested
  IF (abs_coupling .EQ. 1) THEN
    ALLOCATE(Qc_evec(m/2,m/2))
    ALLOCATE(Qc_abs(m/2,m/2))
    ALLOCATE(Qc_sqr(m/2,m/2))

    Qc_sqr = Q_temp(1:m/2,1:m/2)


    CALL diagonalize_m(m/2,Qc_sqr,Qc_eval,Qc_evec)
    Qc_sqr = 0.0d0

    DO nc = 1, m/2
      Qc_sqr(nc,nc) = ABS(REAL(REAL(Qc_eval(nc))))
    END DO
    ! Reassemble into the absolute value, named Qc_sqr because it
    ! used to be the square root of the square but this is a neater construction
    Qc_sqr = MATMUL(Qc_evec,MATMUL(Qc_sqr,CONJG(TRANSPOSE(Qc_evec))))
  END IF

  H_m = 0.0d0
  H_trans_m = 0.0d0
  H = 0.0d0
  Q_c = 0.0d0
  Q_t = 0.0d0
  ! Set up full H as in Faraday. Discus. 
  DO k = 1, 2
    DO k_p = 1, 2
      DO nc = 1, dim_nc
        DO nc_p = 1, dim_nc
          DO nt = 1, dim_nt
            DO nt_p = 1, dim_nt
              index_1 = (k-1)*(dim_nc*dim_nt) + (nc-1)*dim_nt + nt
              index_2 = (k_p-1)*(dim_nc*dim_nt) + (nc_p-1)*dim_nt + nt_p
              IF (abs_coupling .EQ. 1) THEN   ! Couple to |Q_c|
                H_m(index_1,index_2) = del(k,k_p)*del(nc,nc_p)*del(nt,nt_p)*  &
                     (exe(k) + hbar*omegas_c*(nc - 0.5) + hbar*omegas_t*(nt - 0.5)) & 
                      + del(k,k_p)*del(nc,nc_p)*kappa_t(k)*(del(nt,nt_p+1)*SQRT((nt_p)/2.0d0) + &
                      del(nt,nt_p-1)*SQRT((nt_p - 1.0)/2.0d0)) + (del(k,k_p-1) + del(k,k_p+1))*&
                      del(nt,nt_p)*lambda_s*Qc_sqr(nc,nc_p)  ! Get the absolute coupling matrix entry
              ELSE  ! Couple to Q_c
                H_m(index_1,index_2) = del(k,k_p)*del(nc,nc_p)*del(nt,nt_p)*  &
                     (exe(k) + hbar*omegas_c*(nc - 0.5) + hbar*omegas_t*(nt - 0.5)) & 
                      + del(k,k_p)*del(nc,nc_p)*kappa_t(k)*(del(nt,nt_p+1)*SQRT((nt_p)/2.0d0) + &
                      del(nt,nt_p-1)*SQRT((nt_p - 1.0)/2.0d0)) + (del(k,k_p-1) + del(k,k_p+1))*&
                      del(nt,nt_p)*lambda_s*Q_temp(nc,nc_p)

              END IF

              ! Coupling mode position
              Q_c(index_1,index_2) = del(k,k_p)*del(nt,nt_p)*Q_temp(nc,nc_p) 
              ! Tuning mode position
              Q_t(index_1,index_2) = del(k,k_p)*del(nc,nc_p)*Q_temp(nt,nt_p)
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO


  H_temp = H_m


  IF (PRESENT(Qt_full)) THEN
    Qt_full = Q_t
  END IF

  CALL diagonalize_m(m,H_m,evals_m,evec_m)

  evec = evec_m(1:m,1:n)
  eval = evals_m(1:n)
  H = 0.0d0
  DO k = 1, n
    H(k,k) = eval(k)
  END DO
  evec_inv = CONJG(TRANSPOSE(evec))
  ! Put global Q matrices in the eigenbasis

  CALL to_basis(Q_c,Q_x_g,evec,evec_inv)
  CALL to_basis(Q_t,Q_y_g,evec,evec_inv)

  ! Setup initial density matrix with either thermal or vertical
  ! state in either well 2 or well 1
  IF (init_well .EQ. 1) THEN
    bound1 = 1
    bound2 = m/2
  ELSE ! Initialize in well 2
    bound1 = m/2 + 1
    bound2 = m
  END IF

  sigma_m = 0.0d0
  total = 0.0d0
  estate_probs = 0.0d0
  H_trans = 0.0d0
  ! Stationary thermal initalization initializes lowest energy eigenstate at the moment
  IF (init_stat_therm .EQ. 1) THEN
    H_trans = 0.0d0
    H_trans(init_es,1) = 1.0d0
    estate_probs = 0.0d0
    estate_probs(1) = 1.0d0
    sigma(init_es,init_es) = 1.0d0
  ELSE
    estate_probs = 0.0d0
    estate_probs(1) = 1.0d0
    H_trans_m(bound1,1) = 1.0d0
    H_trans(1:n,1) = MATMUL(evec_inv,H_trans_m(1:m,1))
    CALL wfn_to_rho_2(H_trans(1:n,1),sigma)
    sigma_m(bound1,bound1) = 1.0d0
  END IF
   

  IF (init_stat_therm .EQ. 0) THEN
    WRITE(*,*) "Sigma prior to transform"
    CALL print_mat(sigma_m)
    CALL to_basis(sigma_m,sigma,evec,evec_inv)
    WRITE(*,*) "Sigma after transform"
    CALL print_mat(sigma)
  END IF

  DO k = 1, n
    IF (REAL(REAL(sigma(k,k))) .LT. -1.0d-8) THEN
      WRITE(*,*) "Error: negative diagonal entry in density matrix."
      STOP
    END IF
  END DO

  ! Set two separate baths
  couple_op(1:n,1:n,1) = Q_x_g
  couple_op(1:n,1:n,2) = Q_y_g

  proj1m = 0.0d0
  proj2m = 0.0d0
  DO i = 1, m/2
    proj1m(i,i) = 1.0d0 
  END DO
  DO i = m/2 + 1, m
    proj2m(i,i) = 1.0d0
  END DO
  CALL to_basis(proj1m,proj1n,evec,evec_inv)
  CALL to_basis(proj2m,proj2n,evec,evec_inv)

  ! For use with the non-secular trajectory method for plotting purposes
  ALLOCATE(proj1_g(n,n))
  ALLOCATE(proj2_g(n,n))
  proj1_g = proj1n
  proj2_g = proj2n

  OPEN(62,FILE="eigenvector_info.txt")
  OPEN(63,FILE="plot_ev.txt")
  OPEN(64,FILE="evic.txt")
  ! Write information about the eigenvectors
  CALL init_wfn(wfn)
  DO i = 1, n
    wfn%fe = 0.0d0
    wfn%fe(i,1) = 1.0d0
    WRITE(64,*) op_avg(wfn,couple_op(1:n,1:n,1)), op_avg(wfn,couple_op(1:n,1:n,2)), op_avg(wfn,H), &
                op_avg(wfn,proj1n),  op_avg(wfn,proj2n)

    WRITE(62,*) "Eigenvector ", i, "pos1", op_avg(wfn,couple_op(1:n,1:n,1)), "pos2", op_avg(wfn,couple_op(1:n,1:n,2)),&
               "energy", op_avg(wfn,H), "p in 1", op_avg(wfn,proj1n), "p in 2", op_avg(wfn,proj2n)
    WRITE(63,*) i, 0.0, op_avg(wfn,H)
    WRITE(63,*) i, 100.0, op_avg(wfn,H)
    WRITE(63,*) ""
    WRITE(63,*) ""

  END DO
  CALL destroy_wfn(wfn)
  CLOSE(62)
  CLOSE(63)
  CLOSE(64)


  bin_highest = eval(n)
  bin_lowest = eval(1)
  bin_width = (bin_highest - bin_lowest)/99.0d0
  bins = 0.0d0
  WRITE(*,*) "highest", eval(n), "lowest", eval(1), "bin width", bin_width
  DO i = 1, n
    binind = FLOOR(REAL(REAL((eval(i) - eval(1))/(1.0d0*bin_width)))) + 1 
    bins(binind) = bins(binind) + 1.0d0
  END DO
  bins = bins/(1.0d0*n)

  OPEN(65,FILE="dos.txt")
  DO i = 1, 100
    WRITE(65,*) bin_width*(i-1) + REAL(REAL(eval(1))), bins(i)
    WRITE(65,*) bin_width*(i) + REAL(REAL(eval(1))), bins(i)
  END DO
  


  DEALLOCATE(H_m)
  DEALLOCATE(sigma_m)
  DEALLOCATE(evec_m)
  DEALLOCATE(evec_inv)
  DEALLOCATE(H_trans_m)
  DEALLOCATE(proj1n)
  DEALLOCATE(proj2n)
  DEALLOCATE(Q_c)
  DEALLOCATE(Q_t)
  DEALLOCATE(H_temp)
  DEALLOCATE(proj1m)
  DEALLOCATE(tprobs)
  DEALLOCATE(tvecs)
  DEALLOCATE(proj2m)

  IF (ALLOCATED(Qc_sqr)) THEN
    DEALLOCATE(Qc_sqr)
  END IF
  IF (ALLOCATED(Qc_evec)) THEN
    DEALLOCATE(Qc_evec)
  END IF
  IF (ALLOCATED(Qc_abs)) THEN
    DEALLOCATE(Qc_abs)
  END IF
END SUBROUTINE

! Set a Shin-Metiu Hamiltonian describing an adiabatic setup
! where a proton/electron pair is confined in an optical cavity
! with two charges resulting in a quartic PES with the dipole
! moment serving as the coupling between the proton and photon modes
SUBROUTINE setup_H_shin_metiu(H,sigma,evec,eval,H_trans,estate_probs,couple_op,Qt_full)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: H, sigma  ! Hamiltonian, density matrix
  COMPLEX*16, DIMENSION(m,n), INTENT(OUT) :: evec  ! H eigenvectors
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: couple_op  ! System-bath coupling operator
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval  ! H evals
  ! Boltzman likelihood to start in oscillator eigenstate
  REAL*8, INTENT(OUT), DIMENSION(n) :: estate_probs  ! Initial probabilities for "vertical" excitation
  COMPLEX*16, OPTIONAL, INTENT(OUT), DIMENSION(m,m) :: Qt_full  ! Full H prior to basis trimming
  ! Translation of initial eigenvectors for "vertical excitation to the unshifted oscillator basis
  ! looks like this; each vector is its own column
  COMPLEX*16, INTENT(OUT), DIMENSION(n,n) :: H_trans
  ! Raising/lowering, q, q^2, p, p^2 operators for the photon mode
  ! Potential, matter-light coupling component for matter, kinetic energy and lambshift for the
  ! proton mode
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: raise_op, lower_op, q_op, q2_op, mom_op, &
    mom2_op, R_V, R_Vc, R_ke
  ! Full size hamiltonian, temporaries and identity matrices for performing tensor products 
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H_m, temp1, temp2, ident_qc, ident_R
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: evec_m, evec_inv
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: evals_m
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: evals_m_2d
  INTEGER, DIMENSION(1) :: jtemp
  INTEGER :: i, j, sub
  TYPE(wfunc) :: wfn
  ! Full size density matrix and vertical excitation start vectors
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma_m, H_trans_m
  COMPLEX*16 :: temp
  REAL*8 :: entanglement

  ALLOCATE(raise_op(dim_qc,dim_qc))
  ALLOCATE(lower_op(dim_qc,dim_qc))
  ALLOCATE(q_op(dim_qc,dim_qc))
  ALLOCATE(q2_op(dim_qc,dim_qc))
  ALLOCATE(mom_op(dim_qc,dim_qc))
  ALLOCATE(mom2_op(dim_qc,dim_qc))
  ALLOCATE(R_Vc(dim_R,dim_R))
  ALLOCATE(R_V(dim_R,dim_R))
  ALLOCATE(R_ke(dim_R,dim_R))
  ALLOCATE(H_trans_m(m,m))
  ALLOCATE(H_m(m,m))
  ALLOCATE(sigma_m(m,m))

  CALL fill_raise(dim_qc,raise_op)
  CALL fill_lower(dim_qc,lower_op)

  q_op = SQRT(hbar/(2.0d0*omegasc))*(raise_op + lower_op)
  mom_op = SQRT(hbar*omegasc/2.0d0)*(raise_op - lower_op)*CMPLX(0.0d0,1.0d0)
  q2_op = MATMUL(q_op,q_op)
  mom2_op = MATMUL(mom_op,mom_op)

  ! Set up R coordinate sinc function basis ke and V matrices
  R_V = 0.0d0
  R_Vc = 0.0d0
  CALL sinc_ke(dx_R,dim_R,mass_R,R_ke)
  sub = CEILING(1.0d0*dim_R/2.0d0)
  DO i = 1, dim_R
    R_V(i,i) = sinc_V(i - sub,dx_R)
    R_Vc(i,i) = sinc_V_cross(i - sub,dx_R)
  END DO

  ALLOCATE(evals_m(m))
  ALLOCATE(evals_m_2d(m,1))
  ALLOCATE(evec_m(m,m))
  ALLOCATE(evec_inv(n,m))

  ALLOCATE(temp1(m,m))
  ALLOCATE(temp2(m,m))
  ALLOCATE(ident_qc(1:dim_qc,1:dim_qc))
  ALLOCATE(ident_R(1:dim_R,1:dim_R))  ! m = dim_qc*dim_R


  ident_qc = 0.0d0
  ident_R = 0.0d0
  DO i = 1, dim_qc
    ident_qc(i,i) = 1.0d0
  END DO 
  DO i = 1, dim_R
    ident_R(i,i) = 1.0d0
  END DO

  ! Add in the kinetic and potential term from the R sinc basis
  CALL tensor_product(dim_qc,dim_R,ident_qc,R_V + R_ke,temp1)
  H_m = temp1 
  ! Add in the kinetic and potential term from the qc HO basis
  CALL tensor_product(dim_qc,dim_R,mom2_op/2.0d0 + 0.5d0*omegasc**2.0d0*q2_op,ident_R,temp1)
  H_m = H_m + temp1
  ! Add in the coupling term between R and qc
  CALL tensor_product(dim_qc,dim_R,q_op,R_Vc,temp1)
  H_m = H_m + temp1

  IF (RESTART .EQV. .FALSE.) THEN
    CALL diagonalize_m(m,H_m,evals_m,evec_m)

    DO i = 1, n  ! Hard code in that the largest element of every eigenvector shall have its larger real/imag component positive
      jtemp = MAXLOC(ABS(evec_m(1:m,i:i)),1) 
      WRITE(*,*) "maxloc", jtemp
      j = jtemp(1)
      IF (ABS(AIMAG(evec_m(j,i))) .LT. ABS(REAL(REAL(evec_m(j,i)))) .AND. REAL(REAL(evec_m(j,i))) .LT. 0) THEN
        evec_m(1:m,i) = -1.0d0*evec_m(1:m,i)
      ELSE IF (ABS(AIMAG(evec_m(j,i))) .GT. ABS(REAL(REAL(evec_m(j,i)))) .AND. AIMAG(evec_m(j,i)) .LT. 0) THEN
        evec_m(1:m,i) = -1.0d0*evec_m(1:m,i)
      END IF
    END DO


    evals_m_2d(1:m,1) = evals_m(1:m)
    OPEN(22,FILE="restart.dat")
    CALL print_mat_C16(evals_m_2d,22)
    CALL print_mat_C16(evec_m,22)
    CLOSE(22) 
  ELSE
    OPEN(22,FILE="restart.dat")
    i = read_matrix(evals_m_2d,22)
    i = i + read_matrix(evec_m,22)
    evals_m(1:m) = evals_m_2d(1:m,1)
    IF (i .NE. 2) THEN
      WRITE(*,*) "Failed to read in eigenvalues/vectors for restart"
      STOP
    END IF 
    CLOSE(22) 
  END IF
  
  evec = evec_m(1:m,1:n)
  eval = evals_m(1:n)
  H = 0.0d0
  DO i = 1, n
    H(i,i) = eval(i)
  END DO
  WRITE(*,*) "Eigenvalues ", eval
  evec_inv = CONJG(TRANSPOSE(evec))
  ! The coupling operator is linear in R?
  ! Single bath system I believe
  DO i = 1, dim_R
    R_V(i,i) = dx_R*(i - sub)
  END DO
  CALL tensor_product(dim_qc,dim_R,ident_qc,R_V,temp1)
  IF (PRESENT(Qt_full)) THEN
    Qt_full = temp1
  END IF
  CALL to_basis(temp1,Q_x_g,evec,evec_inv)
  CALL tensor_product(dim_qc,dim_R,q_op,ident_R,temp1)
  CALL to_basis(temp1,Q_y_g,evec,evec_inv)


  ! The specified init_es is presumed to be an eigenstate
  IF (init_stat_therm .EQ. 1) THEN

    sigma = 0.0d0
    H_trans = 0.0d0
    estate_probs = 0.0d0
    sigma(init_es,init_es) = 1.0d0
    estate_probs = 1.0d0
    estate_probs(init_es) = 1.0d0
    H_trans(init_es,1) = 1.0d0

  ELSE  ! Use the init_es as a DVR position index and diagonalize
    estate_probs = 0.0d0
    estate_probs(1) = 1.0d0
    ident_R = 0.0d0  ! Reused as a temporary
    ident_qc = 0.0d0  ! Reused as a temporary
    ident_qc(1,1) = 1.0d0  ! Lowest HO state
    ident_R(init_es,init_es) = 1.0d0  ! Selected DVR center state
    CALL tensor_product(dim_qc,dim_R,ident_qc,ident_R,sigma_m)  ! Full proof but simpler ways exist
    CALL to_basis(sigma_m,sigma,evec,evec_inv)  ! Complete formation of sigma
    H_trans_m = 0.0d0
    H_trans_m(dim_qc*(init_es - 1) + 1, 1) = 1.0d0
    H_trans(1:n,1) = MATMUL(evec_inv,H_trans_m(1:m,1))
    OPEN(32,FILE="init_wfn.txt")
    WRITE(32,*) org_out(n,H_trans(1:n,1))
    CLOSE(32)
    temp = 0.0d0
    DO i = 1, n
      temp = temp + sigma(i,i)
    END DO
    WRITE(*,*) "temp is", temp, "ABS", ABS(temp)
    sigma = sigma/temp  ! Normalize the trace
  END IF


  couple_op(1:n,1:n,1) = Q_x_g
  couple_op(1:n,1:n,2) = Q_y_g 

  ALLOCATE(proj1_g(n,n))  ! These are zero as there is only one electronic state
  ALLOCATE(proj2_g(n,n))
  proj1_g = 0.0d0
  proj2_g = 0.0d0
 

  ! Write out some wavefunction information
  OPEN(63,FILE="plot_ev.txt")
  OPEN(64,FILE="evic.txt")
  ! Write information about the eigenvectors
  CALL init_wfn(wfn)
  DO i = 1, n
    entanglement = 0.0d0
    DO j = 1, m
      entanglement = entanglement - evec(j,i)**2.0d0*LOG(evec(j,i)**2.0d0)
    END DO
    wfn%fe = 0.0d0
    wfn%fe(i,1) = 1.0d0
    WRITE(64,*) i, op_avg(wfn,Q_x_g), op_avg(wfn,Q_y_g), op_avg(wfn,H)
    WRITE(63,*) i, 0.0, op_avg(wfn,H)
    WRITE(63,*) i, 100.0, op_avg(wfn,H)
    WRITE(63,*) ""
    WRITE(63,*) ""

    WRITE(*,*) "Eigenvector ", i, "is", evec(1:m,i)
  END DO
  CALL destroy_wfn(wfn)
  CLOSE(62)
  CLOSE(63)
  CLOSE(64)



  DEALLOCATE(H_m)
  DEALLOCATE(evec_m)
  DEALLOCATE(evals_m)
  DEALLOCATE(evec_inv)
  DEALLOCATE(raise_op)
  DEALLOCATE(lower_op)
  DEALLOCATE(q_op)
  DEALLOCATE(q2_op)
  DEALLOCATE(mom_op)
  DEALLOCATE(mom2_op)
  DEALLOCATE(R_Vc)
  DEALLOCATE(R_V)
  DEALLOCATE(R_ke)
  DEALLOCATE(H_trans_m)
  DEALLOCATE(sigma_m)

END SUBROUTINE



END MODULE



