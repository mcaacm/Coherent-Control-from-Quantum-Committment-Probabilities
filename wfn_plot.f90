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


! Reads specifications in the file in_wfns.txt and eigenstate
! wavefunction specifications in params.f90 and plots the wavefunctions
! in CI_data/wfn_plots/
! Also looks at density matrices in "density_matrices.txt"

PROGRAM plot_wf

USE parameters
USE prequel
USE set_up_H
USE rk_ns_utils
USE qns

IMPLICIT NONE
INTEGER :: i, j, k, reason  ! Loop indexes, read error indicator
CHARACTER(LEN=100) :: fname   ! File name
CHARACTER(LEN=20) :: fstring  ! More file name stuff
REAL*8, DIMENSION(sz,sz) :: hf  ! Hermite polynomial factors
REAL*8, DIMENSION(:,:), ALLOCATABLE :: grid  ! Hermite evaluation on grid points
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma, H, evec  ! Denisty matrix, Hamiltonian, energy eigenvectors
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:,:) :: couple_op  ! Coupling operator
REAL*8, DIMENSION(n) :: estate_probs  ! Unneeded setup that used to be important and could be again
COMPLEX*16, DIMENSION(n) :: eval  ! Energy eigenvalues
COMPLEX*16, DIMENSION(m,1) :: wfn_m  ! Wavefunctions in full vs reduced basis
COMPLEX*16, DIMENSION(n,1) :: wfn
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: translated_evecs  ! Unneeded setup
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: plot  ! For wavefunction plotting
REAL*8 :: low, high, lowDVR, highDVR, integral  ! Plot grid limits and an integration temporary
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: out_grid  
REAL*8, DIMENSION(2*n) :: temp
! For an attempt to plot the density matrix wavefunction
COMPLEX*16, DIMENSION(n,n) :: dmat, dmat_evecs
COMPLEX*16, DIMENSION(n) :: dmat_evals
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: plot_acc  ! Grid to accumulate density for printing out the density matrix wavefunction
INTEGER :: iteration, success ! For plotting density matrices
REAL*8 :: full_density


IF (bound1 .LT. 0 .OR. bound2 .LT. 0 .OR. bound2 .LT. bound1) THEN
  WRITE(*,*) "Invalid plot bounds:", bound1, bound2
  STOP
END IF
IF (bound3 .LT. 0 .OR. bound4 .LT. 0 .OR. bound4 .LT. bound3) THEN
  WRITE(*,*) "Invalid plot bounds:", bound3, bound4
  STOP
END IF

ALLOCATE(plot(res,res,6))
ALLOCATE(grid(sz,res))
ALLOCATE(out_grid(res,2))

ALLOCATE(sigma(n,n))
ALLOCATE(H(n,n))
ALLOCATE(evec(m,n))
ALLOCATE(couple_op(n,n,nbath))
ALLOCATE(translated_evecs(n,n))
! Allocate global storage found in prequel
ALLOCATE(G_min_global(n,n,nbath))
ALLOCATE(G_pls_global(n,n,nbath))
ALLOCATE(G_global(n,n,nbath))
! Allocate global storage in Setup_H
ALLOCATE(Q_x_g(n,n))
ALLOCATE(Q_y_g(n,n))


! Setup system
CALL setup(H,sigma,evec,eval,translated_evecs,estate_probs,couple_op)

! Fill a grid of Hermite polynomials for the basis
low = -6.0d0 !-70.0d0
high = 6.0d0 !70.0d0
lowDVR = -6.0d0 !-pi
highDVR = 6.0d0 !pi
CALL fill_hermite(sz,hf)
CALL fill_wf_grid(sz,res,hf,grid,low,high)
WRITE(*,*) "Filled hermite grid"
WRITE(*,*) hf

! Read specific wavefunction data from "in_wfns.txt" and plot as
! wfn_plots/wfn_2d_i.txt
i = 0
OPEN(11,FILE="in_wfns.txt")
DO
  READ(11,*,IOSTAT=reason) (temp(k), k=1, 2*n)

  IF (reason .NE. 0) THEN
    WRITE(*,*) "Read error:", reason, "exiting loop after", i, "reads"
    EXIT
  END IF

  DO j = 1, n
    wfn(j,1) = CMPLX(temp(2*j - 1),temp(2*j))
  END DO
  WRITE(*,*) "WFN to plot ", wfn
  wfn_m = MATMUL(evec,wfn)
  WRITE(*,*) "Basis conversion ", wfn_m
  WRITE(*,*) "m/2+1", wfn_m(m/2 - 1:m/2 + 1,1)
  !wfn_m = 0.0d0
  !wfn_m(m/2 + 1,1) = 1.0d0
  IF (run_type .EQ. CI) THEN
    CALL fill_wf_plot(sz,res,grid,low,high,wfn_m,plot)
  ELSE
    CALL fill_wf_plot_ho_sinc(sz,res,grid,low,high,lowDVR,highDVR,wfn_m,plot(1:res,1:res,1:3))
  END IF

  IF (i .LT. 10) THEN
    fstring = "(A16,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A16,I2,A4)"
  ELSE 
    fstring = "(A16,I3,A4)"
  END IF 

  IF (run_type .EQ. CI) THEN  ! Conical intersection HOxHO basis
    WRITE(fname,fstring) "wfn_plots/wf_2d_", i, ".txt"
    OPEN(43,FILE=TRIM(fname))
    integral = integrate_grid(res,plot,low,high)
    CALL print_wf_plot(43,plot,low,high,res)
    WRITE(43,*) "Integrated to", integral
    CLOSE(43)
  ELSE  ! Polariton SINCxHO basis
    WRITE(fname,fstring) "wfn_plots/wf_2d_", i, ".txt"
    OPEN(43,FILE=TRIM(fname))
    CALL print_wf_plot_ho_sinc(43,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR,res)
    integral = integrate_grid_ho_sinc(res,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR)
    WRITE(43,*) "Integrated to", integral
    CLOSE(43)
  END IF

  i = i + 1
END DO


! Plot lower energy eigenstates of relevance
DO i = bound1, bound2
  wfn = 0.0d0
  wfn(i,1) = 1.0d0
  wfn_m = MATMUL(evec,wfn)
  ! Fill and integrate a grid for each eigenstate
  IF (run_type .EQ. CI) THEN
    CALL fill_wf_plot(sz,res,grid,low,high,wfn_m,plot)
    integral = integrate_grid(res,plot,low,high)
  ELSE  ! Slightly different routines required for the SINC function
    CALL fill_wf_plot_ho_sinc(sz,res,grid,low,high,lowDVR,highDVR,wfn_m,plot(1:res,1:res,1:3))
    integral = integrate_grid_ho_sinc(res,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR)
  END IF
  WRITE(*,*) "Integrated density eigen fn ", i, "is", integral


  IF (i .LT. 10) THEN
    fstring = "(A16,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A16,I2,A4)"
  ELSE 
    fstring = "(A16,I3,A4)"
  END IF 


  WRITE(fname,fstring) "wfn_plots/ef_2d_", i, ".txt"
  IF (run_type .EQ. CI) THEN 
    OPEN(43,FILE=TRIM(fname))
    CALL print_wf_plot(43,plot,low,high,res)
    CLOSE(43)
  ELSE
    OPEN(43,FILE=TRIM(fname))
    CALL print_wf_plot_ho_sinc(43,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR,res)
    CLOSE(43)
  END IF
END DO


! Check normalization of individual eigenstates
! Plot higher energy eigenstates of relevance
DO i = bound3, bound4
  wfn = 0.0d0
  wfn(i,1) = 1.0d0
  wfn_m = MATMUL(evec,wfn)
  IF (run_type .EQ. CI) THEN
    CALL fill_wf_plot(sz,res,grid,low,high,wfn_m,plot)
    integral = integrate_grid(res,plot,low,high)
  ELSE
    CALL fill_wf_plot_ho_sinc(sz,res,grid,low,high,lowDVR,highDVR,wfn_m,plot(1:res,1:res,1:3))
    integral = integrate_grid_ho_sinc(res,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR)
  END IF
  WRITE(*,*) "Integrated density eigen fn ", i, "is", integral

  IF (i .LT. 10) THEN
    fstring = "(A16,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A16,I2,A4)"
  ELSE 
    fstring = "(A16,I3,A4)"
  END IF 

  
  WRITE(fname,fstring) "wfn_plots/ef_2d_", i, ".txt"
  IF (run_type .EQ. CI) THEN
    OPEN(43,FILE=TRIM(fname))
    CALL print_wf_plot(43,plot,low,high,res)
    CLOSE(43)
  ELSE
    OPEN(43,FILE=TRIM(fname))
    CALL print_wf_plot_ho_sinc(43,plot(1:res,1:res,1:3),low,high,lowDVR,highDVR,res)
    CLOSE(43)
  END IF

END DO


CLOSE(11)



! While there are more density matrices to read
! Read in a density matrix of size nxn
! Diagonalize the density matrix
! Plot each of the wavefunctions
! Put them together into a weighted average
! Print to a file
OPEN(20,FILE="density_matrices.txt")
success = 1
iteration = 0
ALLOCATE(plot_acc(res,res,6))
DO WHILE (success .GT. 0)
  success = read_matrix(dmat,20)
  CALL diagonalize_m(n,dmat,dmat_evals,dmat_evecs)
  plot_acc = 0.0d0
  full_density = 0.0d0
  DO i = 1, n
    IF (ABS(dmat_evals(i)) .GT. 1d-3) THEN
      wfn(1:n,1) = dmat_evecs(1:n,i)
      wfn_m = MATMUL(evec,wfn)
      CALL fill_wf_plot(sz,res,grid,low,high,wfn_m,plot)
      integral = integrate_grid(res,plot,low,high)
      WRITE(*,*) "INTEGRAL", iteration, "is", integral
     plot_acc = dmat_evals(i)*plot + plot_acc  ! Weight it in
    END IF
    full_density = dmat_evals(i)*integral + full_density
  END DO


  IF (i .LT. 10) THEN
    fstring = "(A16,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A16,I2,A4)"
  ELSE 
    fstring = "(A16,I3,A4)"
  END IF 

  WRITE(fname,fstring) "wfn_plots/pf_2d_", iteration, ".txt"
  WRITE(*,*) "Printing ", i, "to", fname, "integral", full_density, "evals", dmat_evals
  OPEN(43,FILE=fname)
  CALL print_wf_plot(43,plot_acc,low,high,res)
  CLOSE(43)
  iteration = iteration + 1
END DO


END PROGRAM

