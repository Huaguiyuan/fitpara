MODULE banddata
  !
  USE constants, only : dp
  !
  integer nkpt, nbnd
  real(dp), allocatable :: kvec(:, :)
  real(dp), allocatable :: eig(:, :)
  real(dp), allocatable :: wgt(:, :)
  real(dp) ef
 !
 CONTAINS
 !
SUBROUTINE init_banddata()
  !
  IMPLICIT NONE
  !
  allocate(kvec(1:3, 1:nkpt), eig(1:nbnd, 1:nkpt), wgt(1:nbnd, 1:nkpt))
  !
END SUBROUTINE

SUBROUTINE finalize_banddata()
  !
  IMPLICIT NONE
  !
  deallocate(kvec, eig, wgt)
  !
END SUBROUTINE

SUBROUTINE calc_weight(tau)
  !
  USE constants, only : dp
  !
  IMPLICIT NONE
  !
  real(dp) tau
  !
  integer ik, ib
  real(dp) norm
  !
  norm=0.d0
  !
  do ik=1, nkpt
    wgt(:, ik)=exp(-((eig(:, ik)-ef)/tau)**2)
    norm=norm+sum(wgt(:,ik))
  enddo
  !
  wgt(:,:)=wgt(:,:)/norm
  !
END SUBROUTINE

SUBROUTINE read_banddata()
  !
  USE constants, only : dp, fin, stdout
  !
  IMPLICIT NONE
  !
  integer ik, ib, ii
  !
  open(unit=fin, file="fs.dat")
  read(fin, *) nkpt, nbnd, ef
  !
  CALL init_banddata()
  !
  do ik=1, nkpt
    do ib=1, nbnd
      read(fin, *) kvec(:, ik), ii, eig(ib, ik)
      if (ib.ne.ii) then
        write(stdout, *) " !!!! WRONG fs.dat"
      endif
    enddo
  enddo
  !
  close(unit=fin)
  !
END SUBROUTINE
 !
END MODULE
