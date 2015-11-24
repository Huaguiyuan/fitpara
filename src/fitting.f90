MODULE fitting
  !
  USE constants,  only : dp
  !
  IMPLICIT NONE
  !
  integer, parameter ::  npar=18
  real(dp), allocatable :: par(:)
  integer,  allocatable :: opt_switch(:)
  real(dp), allocatable :: parder(:)
  real(dp), allocatable :: parderold(:)
  integer mode
  real(dp) tau, step
  !
 CONTAINS
  !
SUBROUTINE read_param()
  !
  USE constants,  only : fin
  !
  IMPLICIT NONE
  !
  integer ii
  real(dp),dimension(1:10) :: tt
  !
  allocate(par(1:npar), opt_switch(1:npar), parder(1:npar), parderold(1:npar))
  !
  open(unit=fin, file="hopping.dat")
  !
  read(fin, *) step, tau
  !
  do ii=0, 1
    read(fin, *) tt(:)
    par(ii*5+1:ii*5+5)=tt(1:5)
    opt_switch(ii*5+1:ii*5+5)=tt(6:10)
  enddo
  !
  do ii=0, 1
    read(fin, *) tt(1:8)
    par(ii*4+11:ii*4+14)=tt(1:4)
    opt_switch(ii*4+11:ii*4+14)=tt(5:8)
  enddo
  !
  close(unit=fin)
  !
END SUBROUTINE

SUBROUTINE finalize_param()
  !
  IMPLICIT NONE
  !
  deallocate(par, opt_switch, parder, parderold)
  !
END SUBROUTINE

SUBROUTINE write_param()
  !
  USE constants,      only : stdout
  !
  IMPLICIT NONE
  !
  integer ii
  !
  write(stdout, '(2F22.16)') step, tau
  !
  do ii=0, 1
    write(stdout, '(5F16.12,5I3)') par(ii*5+1:ii*5+5), opt_switch(ii*5+1:ii*5+5)
  enddo
  !
  do ii=0, 1
    write(stdout, '(4F16.12,16X,4I3)') par(ii*4+11:ii*4+14), opt_switch(ii*4+11:ii*4+14)
  enddo
  !
END SUBROUTINE
  !
FUNCTION calc_deviation(t)
  !
  USE constants,   only : dp
  USE banddata,    only : kvec, eig, wgt, ef, nbnd, nkpt
  USE model,       only : calc_band_ek
  !
  IMPLICIT NONE
  !
  real(dp) t(:)
  real(dp) calc_deviation
  !
  integer ik, ib
  real(dp), allocatable :: ek_fit(:)
  real(dp), allocatable :: dev(:)
  !
  allocate(ek_fit(1:nbnd))
  allocate(dev(1:nbnd))
  dev(:)=0.d0
  !
  do ik=1, nkpt
    call calc_band_ek(ek_fit, kvec(:, ik), t)
    dev(:)=dev(:)+(ek_fit(:)-eig(:, ik))*(ek_fit(:)-eig(:, ik))*wgt(:, ik)
  enddo
  !
  calc_deviation=SUM(dev(:))
  !
  deallocate(ek_fit,dev)
  !
END FUNCTION

SUBROUTINE calc_partial()
  !
  USE constants,    only : dp, eps8
  !
  IMPLICIT NONE
  !
  integer ip
  real(dp), allocatable :: tt(:)
  real(dp) dev, dev0
  !
  allocate(tt(1:npar))
  !
  dev0=calc_deviation(par)
  !
  do ip=1, npar
    tt(:)=par(:)
    if (opt_switch(ip)>0) then
      tt(ip)=par(ip)+eps8
      dev=calc_deviation(tt)
      parder(ip)=(dev-dev0)/eps8
    else
      parder(ip)=0.d0
    endif
  enddo
  !
  deallocate(tt)
  !
END SUBROUTINE

FUNCTION minimize_dev(max_iter, ireport, eps)
  !
  USE constants,     only : dp, stdout
  !
  IMPLICIT NONE
  !
  integer max_iter
  integer ireport
  real(dp) eps
  real(dp) minimize_dev
  !
  integer ii
  real(dp) norm
  !
  do ii=1, max_iter
    !
    par(:)=par(:)-step*parder(:)
    !
!    CALL fix_chempot()
    !
    CALL calc_partial()
    !
    norm=sqrt(SUM(parder(:)*parder(:)))
    !
    if (mod(ii-1,ireport)==0) then
      write(stdout, '(A,1I8,A,1F14.9,A)') "==Iteration:", ii, "  Norm:", norm, "=="
      write(stdout, *) " Derivatives are: "
      write(stdout, '(5F16.12)') parder(1:10)
      write(stdout, '(4F16.12)') parder(11:18)
      CALL write_param()
    endif
    !
    if (norm<eps) exit
    !
  enddo
  !
  minimize_dev=norm
  !
END FUNCTION
  !
SUBROUTINE fix_chempot()
  !
  IMPLICIT NONE
  !
  par(4)= 0.45582371-3*par(3)
  par(9)=-0.16697780-3*par(8)
  par(1)= 2.19153437-3*par(2)
  par(6)= 0.12980667-3*par(7)
  !
END SUBROUTINE
  !
END MODULE
