PROGRAM test
  !
  USE banddata
  USE model
  USE fitting
  !
  IMPLICIT NONE
  !
  integer ik, ib
  real(dp),dimension(1:3) :: ek
  !
  CALL read_banddata()
  CALL read_param()
  !
  write(*,'(6F16.9)') par(:)
  !
  do ik=1, nkpt
    CALL calc_band_ek(ek, kvec(:, ik), par)
    do ib=1, nbnd
      write(*,'(3F16.9,2X,1F16.9)') kvec(:, ik), ek(ib)
    enddo
  enddo
  !
  CALL finalize_banddata()
  CALL finalize_param()
  !
END PROGRAM
