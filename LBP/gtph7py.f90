subroutine gtph7py(iup,wlc,pr,ph,cext,cabs,cg,vl,nang,ang,err)
! get phase function and cross section 
! from table of interpolated refractive index
! with log-regular grid wavelengths.
!--- history
! 95. 9.14  created from getph3
! 96. 2.29  2nd order polynomial interpolation (t.y.nakajima)
! 96. 3.11  pollack and cuzzi method implemented (noguchi, tokai)
! 96. 4. 4  dcg(ir)=dcg(ir)+cg1 -> f(indr.le.1) then ... else...endif
! 08. 9.19  modified Fortran 90 free-style form
! 10. 3.30  Modified for Ping-Yang column scattering
! 22. 1.22  debuged.  
!--- input
! iup     i      read unit number of the kernel file for ir range.
! wlc     r      wavelength in cm
! pr   r(6,4)    parameter packet of size distribution (see vlspc2)
!--- output
! nang    i      nbr of scattering angle
! req  r(kreq)   equivalent-volume spherical radius
! ang  r(kang)   scattering angle in degree
! ph   r(knang)  volume phase function
! cext   r       extinction cross section (cm-1)
! cabs   r       absorption cross section (cm-1)
! cg     r       geometrical cross section (cm-1)
! vl     r       volume (cm3/cm3)
!---
  use paras, only: knang,knang2,kpol2, kwl,kreq
  implicit none

! for input
  integer,intent(in):: iup        !! device number for PY kernel
  real(8),intent(in):: wlc        !! wavelength [cm]
  real(8),intent(inout):: pr(6,4)   !! parameter packet

! for output
!  real(8),intent(out):: cr         !! refractive index (real part)
!  real(8),intent(out):: ci         !! refractive index (imaginary part)
  real(8),intent(out):: ph(knang,kpol2)   !! volume phase function
  real(8),intent(out):: cext       !! extinction cross section
  real(8),intent(out):: cabs       !! absorption cross section
  real(8),intent(out):: cg         !! geometrical cross section
  real(8),intent(out):: vl         !! volume
  integer,intent(out):: nang       !! number of angles
  real(8),intent(out):: ang(knang)  !! angles
  character,intent(inout):: err*64

! kernel parameters
  integer:: nwl,nreq,nang2,npol2
  integer:: nw
  real(8):: wl(kwl),req(kreq)
  real(8):: dwl,phs(knang2,kpol2,kreq)
  real(4):: phsf(knang2,kpol2,kreq,2)

  integer:: iw,ireq,j,ip
  integer:: nav

  real(8):: r1,r2,r3,dr,dlr,v1,v2,cg2
  real(8):: vlspc2
!--exec
  npol2=kpol2
  nwl=kwl
  nreq=kreq
  nang=knang
!  call pykrnl(wl(1:nwl),req(1:nreq),ang(1:knang))
  call pykrnl(wl,req,ang)
  
  if(wlc<wl(1)) then
!            1234567890123456789012345678901234567890
     err='Out of PY kernel range: Shorter than 0.199 micron.'
     return
  endif
  if(wlc>wl(nwl)) then
     err='Out of PY kernel range: Longer than 99 micron.'
     return
  endif
  
  do iw=2,nwl
     if(wlc<wl(iw)) then
        nw=iw
        dwl=(wl(iw)-wlc)/(wl(iw)-wl(iw-1))
        exit
     endif
  enddo

  nang2=nang+2
  do ireq=1,nreq
     read(iup,rec=(nw-2)*kreq+ireq) phsf(1:nang2,1:npol2,ireq,1)
     read(iup,rec=(nw-1)*kreq+ireq) phsf(1:nang2,1:npol2,ireq,2)
  enddo
  phs(1:nang2,1:npol2,1:nreq)=phsf(1:nang2,1:npol2,1:nreq,1)*dwl &
       +phsf(1:nang2,1:npol2,1:nreq,2)*(1.d0-dwl)

! P??/P11-->P?? (debuged)
  do ip=2,npol2
     phs(1:nang2,ip,1:nreq)=phs(1:nang2,ip,1:nreq)*phs(1:nang2,1,1:nreq)
  enddo
  
! number of averaging
  nav=10

! get interpolated optical constants
  cg=0.d0; vl=0.d0
  cext=0.d0; cabs=0.d0
  ph(1:nang,1:npol2)=0.d0

! integrate size distribution
  do ireq=1,nreq-1
     r1=req(ireq)
     r2=req(ireq+1)
     dr=(r2-r1)/dble(nav)
     r3=r1-dr*0.5d0
     v1=0.d0
     cg=0.d0
     do j=1,nav
        r3=r3+dr
        dlr=log(r3/(r3-dr))
        pr(1,1)=r3
        v2=vlspc2(pr)
        vl=vl+v2*dlr
        v1=v1+v2
!        cg2=v2*3.d0/4.d0/r3
!        cg =cg +cg2*dlr
        cg2=3.d0/(4.d0*r3)
        cg=cg+cg2*dlr
     enddo
     v1=v1/dble(nav)
     cg=cg/dble(nav)
     
     ph(1:nang,1:npol2)=ph(1:nang,1:npol2)+(phs(1:nang,1:npol2,ireq)+phs(1:nang,1:npol2,ireq+1))*cg*v1*0.5d0
     cext=cext+(phs(nang+1,1,ireq)+phs(nang+1,1,ireq+1))*cg*v1*0.5d0
     cabs=cabs+(phs(nang2,1,ireq)+phs(nang2,1,ireq+1))*cg*v1*0.5d0
  enddo
  return
end subroutine gtph7py
