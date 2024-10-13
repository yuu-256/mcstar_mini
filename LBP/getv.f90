subroutine getv(wlc,pr,intvl,szp,vl)
! get total volume
!--- history
! 95. 9.14   created from
! 08. 9.19 modified Fortran 90 free-style form
!--- input
! wlc     r      wavelength in cm
! pr   r(6,4)    parameter packet of size distribution (see vlspc2)
!--- output
! vl     r       volume (cm3/cm3)
!---
  use paras, only: kintvl,pi
  implicit none

! input & output
  real(8),intent(in):: wlc       !! wavelength [cm]
  real(8),intent(inout):: pr(6,4)  !! parameter packet
  integer,intent(in):: intvl
  real(8),intent(in):: szp(kintvl) !! size parameter
  real(8),intent(out):: vl       !! volume

! work
  integer:: j
  integer:: nav
  real(8):: wvn
  real(8):: del,rszp2
  real(8):: rmin,rmin1,r1,r2,r3,dr,dlr
  integer:: intvl9,intvl8
  integer:: ir8,ir
  real(8):: v2,vlspc2
!--exec
  nav=10
  wvn=2.d0*pi/wlc
  del=log(szp(intvl)/szp(1))/(intvl-1)
  rszp2=exp(0.5d0*del)
  rmin=pr(1,3)
  rmin1=szp(1)/wvn/rszp2

  if(rmin1>rmin) then
     intvl9=int(log(rmin1/rmin)/del)+1
     r2=rmin1*exp(-intvl9*del)
  else
     intvl9=0
  endif

  vl=0.d0
  intvl8=intvl+intvl9
  do ir8=1,intvl8
     ir=ir8-intvl9
     if(ir<=0) then
        r1=r2
        r2=r1*rszp2**2
     else
        r1=szp(ir)/wvn/rszp2
        r2=szp(ir)/wvn*rszp2
     endif
     dr=(r2-r1)/dble(nav)
     r3=r1-0.5d0*dr
     do j=1,nav
        r3=r3+dr
        dlr=log(r3/(r3-dr))
        pr(1,1)=r3
        v2=vlspc2(pr)
        vl=vl+v2*dlr
     enddo
  enddo
  return
end subroutine getv
