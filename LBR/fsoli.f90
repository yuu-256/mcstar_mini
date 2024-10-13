subroutine fsoli(wl1,wl2,fsol)
! wavelength integrated solar incidence.
!---history
! 94. 5.30 created by M.Tsukamoto
! 97. 4.19 tuned parameters by T. Nakajima
! 03. 9. 9 Changed for free source form (Fortran 90) by M. Sekiguchi 
!---input
! wl1: min. wavelength of the band (micron)
! wl2: max. wavelength of the band (micron)
!---output
! fsol: integrated solar irradiance (w/m2)
!---
  implicit none
  real(8),intent(in)::wl1,wl2
  real(8),intent(out)::fsol
  integer,parameter::knx=3000,knm=10
  real(8),parameter::dwwlc=1.d6
  integer::n,i
  real(8)::wll1,wll2,wl,dwll,wn,sunir
!--exec
  wll1=log(wl1); wll2=log(wl2)
  n=min(int((wll2-wll1)*dwwlc),knx); n=max(n,knm)
  dwll=(wll2-wll1)/dble(n)
  fsol=0.d0
  do i=1,n
     wl=exp(wll1+(dble(i)-0.5d0)*dwll); wn=1.d4/wl
     fsol=fsol+wl*sunir(wn)*dwll
  enddo
  return
end subroutine fsoli
