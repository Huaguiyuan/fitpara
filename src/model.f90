#include 'lapack.f90'

MODULE model
  !
  USE constants,    ONLY: dp, twopi, cmplx_1, cmplx_0, cmplx_i, sqrt3
  !
  complex(dp), dimension(1:3, 1:3, 0:8) :: lambda(1:3, 1:3, 0:8)=reshape( &
    (/ cmplx_1, cmplx_0, cmplx_0 ,  cmplx_0, cmplx_1, cmplx_0 ,  cmplx_0, cmplx_0, cmplx_1,  &
       cmplx_0, cmplx_1, cmplx_0 ,  cmplx_1, cmplx_0, cmplx_0 ,  cmplx_0, cmplx_0, cmplx_0,  &
       cmplx_0,-cmplx_i, cmplx_0 ,  cmplx_i, cmplx_0, cmplx_0 ,  cmplx_0, cmplx_0, cmplx_0,  &
       cmplx_1, cmplx_0, cmplx_0 ,  cmplx_0,-cmplx_1, cmplx_0 ,  cmplx_0, cmplx_0, cmplx_0,  &
       cmplx_0, cmplx_0, cmplx_1 ,  cmplx_0, cmplx_0, cmplx_0 ,  cmplx_1, cmplx_0, cmplx_0,  &
       cmplx_0, cmplx_0,-cmplx_i ,  cmplx_0, cmplx_0, cmplx_0 ,  cmplx_i, cmplx_0, cmplx_0,  &
       cmplx_0, cmplx_0, cmplx_0 ,  cmplx_0, cmplx_0, cmplx_1 ,  cmplx_0, cmplx_1, cmplx_0,  &
       cmplx_0, cmplx_0, cmplx_0 ,  cmplx_0, cmplx_0,-cmplx_i ,  cmplx_0, cmplx_i, cmplx_0,  &
       cmplx_1/sqrt3, cmplx_0, cmplx_0 ,  cmplx_0, cmplx_1/sqrt3, cmplx_0 ,  cmplx_0, cmplx_0,-2.d0*cmplx_1/sqrt3 /), (/3,3,9/) )
  !
 CONTAINS
  !
SUBROUTINE calc_band_ek(ek, k, t)
  !
  USE constants, only : dp, sqrt3, cmplx_0
  USE lapack95,  only : heev
  !
  IMPLICIT NONE
  !
  real(dp),dimension(1:3) :: ek
  real(dp),dimension(1:3) :: k
  real(dp) t(:)
  !
  real(dp) cos_ka, cos_kb, cos_kc
  real(dp) sin_ka, sin_kb, sin_kc
  real(dp) cos_kz, cos_2kz
  !
  real(dp),dimension(0:8) :: ksi
  complex(dp),dimension(1:3,1:3) :: ham
  !
  integer ii, info;
  !
  cos_kz=cos(k(3))
  cos_2kz=cos(2*k(3))
  !
  cos_ka=cos(k(1))
  cos_kb=cos(k(2))
  cos_kc=cos(k(1)+k(2))
  !
  sin_ka=sin(k(1))
  sin_kb=sin(k(2))
  sin_kc=-sin(k(1)+k(2))
  !
  ksi(0)=t(1)+(t(2)+t(3)*cos_kz)*(cos_ka+cos_kb+cos_kc)+t(4)*cos_kz+t(5)*cos_2kz;
  ksi(8)=t(6)+(t(7)+t(8)*cos_kz)*(cos_ka+cos_kb+cos_kc)+t(9)*cos_kz+t(10)*cos_2kz;

  ksi(2)=-(t(11)+t(12)*cos_kz)*(sin_ka+sin_kb+sin_kc);

  ksi(1)=sqrt3*(t(13)+t(14)*cos_kz)*(cos_ka-cos_kb);
  ksi(3)=(t(13)+t(14)*cos_kz)*(cos_ka+cos_kb-2*cos_kc);

  ksi(4)=(t(15)+t(16)*cos_kz)*(cos_ka+cos_kb-2*cos_kc);
  ksi(6)=sqrt3*(t(15)+t(16)*cos_kz)*(cos_kb-cos_ka);

  ksi(5)=sqrt3*(t(17)+t(18)*cos_kz)*(sin_ka-sin_kb);
  ksi(7)=(t(17)+t(18)*cos_kz)*(sin_ka+sin_kb-2*sin_kc);
  !
  ham(:,:)=cmplx_0
  !
  do ii=0, 8
    ham=ham+ksi(ii)*lambda(:, :, ii)
  enddo
  !
  call heev(ham, ek, 'N', 'U', info)
  !
END SUBROUTINE
 !
END MODULE
