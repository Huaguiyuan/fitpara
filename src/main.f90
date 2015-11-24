PROGRAM fit_param
  !
  USE constants,    only : dp, eps6, stdout
  USE fitting
!  USE fitting,      only : minimize_dev, calc_deviation, read_param, calc_partial, finalize_param, write_param, par, tau
  USE banddata,     only : read_banddata, finalize_banddata, calc_weight
  !
  IMPLICIT NONE
  !
  integer niter, iter
  integer method, iflag, irest
  logical finish
  integer,dimension(1:2) :: iprint
  real(dp) dev, norm, tlev
  real(dp), allocatable :: d(:), w(:)
  character(len=20) arg
  !
  iprint(1) = 1
  iprint(2) = 10
  irest =  1
  method = 3
  allocate(d(1:npar), w(1:npar))
  !
  CALL getarg(1, arg)
  read(arg, *) niter
  !
  CALL read_param()
  !
  CALL read_banddata()
  !
  CALL calc_weight(tau)
  !
  dev=calc_deviation(par)
  !
  write(stdout, '(A,1F16.9)') " Initial deviation:", dev
  !
!  norm=minimize_dev(niter, 10, eps6)
  do iter=1, niter
    !
    if (iflag.eq.1 .or. iter.eq.1) then
      dev=calc_deviation(par)
      CALL calc_partial()
    endif
    !
    CALL cgfam(npar, par, dev, parder, d, parderold, iprint, eps6, w, &
              iflag, irest, method, finish)
    !
    norm=SUM(par(:)*par(:))
    !
    if (iflag.le.0) exit
    !
    if (iflag.eq.2) then
      if (norm<eps6) then
        finish=.true.
      endif
    endif
    !
  enddo
  !
  dev=calc_deviation(par)
  !
  if (norm<eps6) then
    write(stdout, *) "  Convergence Achieved!!!"
  else
    write(stdout, *) "  Convergence Not Achieved!!!"
  endif
  !
  write(stdout, '(A,1F16.9)') " Final deviation:", dev
  !
  write(stdout, *) "  Final Parameters:"
  !
  CALL write_param()
  !
  CALL finalize_param()
  CALL finalize_banddata()
  !
END PROGRAM
