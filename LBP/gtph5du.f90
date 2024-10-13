subroutine gtph5du(iuk,wlc,cr3,ci3,pr,nang,ipol,ang,ph, &
     cext,cabs,cg,vl,intvl,szp,err)
! get phase function and cross section 
! from table of interpolated refractive index
! with log-regular grid wavelengths.
! cr3,ci3 assigned; assume ipol=1
! pollack and cuzzi method implemented
!--- history
! 95. 9.14  created from getph3
! 96. 2.29  2nd order polynomial interpolation (t.y.nakajima)
! 96. 3.11  pollack and cuzzi method implemented (noguchi, tokai)
! 96. 4. 4  dcg(ir)=dcg(ir)+cg1 -> f(indr.le.1) then ... else...endif
! 08. 9.19  modified Fortran 90 free-style form
! 10. 5.17  Fix some bugs.
! 16. 6.28  extended for polarization.
! 22. 1.28  Fixed bugs.  
!--- input
! iuk     i      read unit number of the kernel file.
! wlc     r      wavelength in cm
! cr3     r      real part of refractive index
! ci3     r      imaginary part of refractive index (positive value)
! pr   r(6,4)    parameter packet of size distribution (see vlspc2)
! x0      r      critical size parameter for mie scattering (pollack-cuzzi)
!                give 1.0e9 for mie theory only calculations.
! gg      r      asymmetry parameter for transmitted ray (pollcack-cuzzi)
!                any value for mie theory only calculations.
! rp      r      surface area fraction to that of sphere (about 1.1-1.3)
!                any value for mie theory only calculations.
!--- output
! nang    i      nbr of scattering angle
! intvl   i      nbr of size parameter
! ipol    i      nbr of polarization component
! ang  r(nang)   scattering angle in degree
! ph   r(knang,ipol) volume phase function
! cext   r       extinction cross section (cm-1)
! cabs   r       absorption cross section (cm-1)
! cg     r       geometrical cross section (cm-1)
! vl     r       volume (cm3/cm3)
! intvl  i       number of size interval in the kernel file
! szp  r(intvl)
!---
  use paras, only: knang, kintvl=>kintvl_du,kpol,kpol2,pi,rad
  implicit none

! for input
  integer,intent(in):: iuk       !! device number
  real(8),intent(in):: wlc       !! wavelength [cm]
  real(8),intent(in):: cr3       !! refractive index (real part)
  real(8),intent(in):: ci3       !! refractive index (imaginary part)
  real(8),intent(inout):: pr(6,4)   !! parameter packet

! for output
  integer,intent(in):: nang     !! number of scattering angles
  integer,intent(in):: ipol     !! number of polarization component
  real(8),intent(inout):: ang(knang)   !! scattering angles [degree]
  real(8),intent(out):: ph(knang,kpol2) !! volume phase function
  real(8),intent(out):: cext         !! extinction cross section
  real(8),intent(out):: cabs         !! absorption cross section
  real(8),intent(out):: cg           !! geometrical cross section
  real(8),intent(out):: vl           !! volume
  integer,intent(inout):: intvl        !! number of size interval in kernel file
  real(8),intent(inout):: szp(kintvl)  !! size parameter
  character,intent(inout):: err*64

! initialization
  integer:: init=1
  real(8),save:: del,rszp2,phm(knang,kpol2)

! kernel parameters
  integer,parameter:: kr=15          !! number of real part in each section
  integer,parameter:: ki=15          !! number of imaginary part in each sec
  integer,parameter:: krf=kr*ki
  integer,parameter:: knang2=knang+2

  integer:: i,j,l,ip,jp,ipt,jpt,ir,ix,ii
  integer:: nang1,ipol1
  real(8),save:: cr(kr),ci(ki)
  real(8),save:: q(kintvl,kpol2,kr,ki,knang2)
  integer,save:: nav

  integer:: nintpr,nintpi
  real(8):: xr(3),xi(3)

  integer:: intvl9,intvl8,ir8
  integer:: indr
  real(8):: r1,r2,r3,dr,dlr,v1,v2,cg1,rmin1,wvn,rmin,rmax
  real(8):: vlspc2

  real(8):: cof
  real(8):: xi2(3),y(3),zi(3),z
  real(8):: bintp,pint

! for smllop 
  real(8):: x,qext9,qsca9,qabs9,cs(knang)

! kernel data parameters  
  integer:: nr = 15        !! refractive indices (real)
  integer:: ni = 15        !! refractive indices (imaginary)
!  real(8):: crmn = 1.33d0  !! minimum value of rfi (real) 
!  real(8):: crmx = 1.60d0  !! maximum value of rfi (real) 
!  real(8):: cimn = 5.d-4   !! minimum value of rfi (imaginary)
!  real(8):: cimx = 5.d-1   !! maximum value of rfi (imaginary)
!--exec
  if(init>0) then
     init=0
     rewind iuk
! get kernel ang, szp, ak, cext, cabs
     read(iuk) intvl, nang1, ipol1
     if (nang1 > knang) then
        err = 'nang > knang'; return
     end if
     if(nang1/=nang) then
        err='inappropriate number of angles in Dubovik KRNL'; return
     endif

     if (ipol1 > kpol2) then  !! mm220128
        err = 'ipol > kpol' ; return
     end if
     if(ipol1/=6) then        !! mm220128
        err='assume only ipol=6 case in gtph5du'; return
     endif
     
!mm220128 removed
!     if(ipol1/=ipol) then
!        err='inappropriate KRNL in gtph5b'; return
!     endif

     read(iuk) szp(1:intvl)
     read(iuk) ang(1:nang)

     del=log(szp(intvl)/szp(1))/dble(intvl-1)
     rszp2=exp(del*0.5d0)

!! read kernel
     do ir=1,kr; do ii=1,ki
        read(iuk)
        read(iuk) cr(ir),ci(ii)
        ci(ii)=abs(ci(ii))
        do i = 1, intvl
           do j = 1, ipol1       !! mm220128
              read(iuk) q(i,j,ir,ii,1:nang)
           enddo
        enddo
        read(iuk) q(1:intvl,1,ir,ii,nang+1)
        read(iuk) q(1:intvl,1,ir,ii,nang+2)
     enddo; enddo

! phase function for small particles
     cs(1:nang)=cos(ang(1:nang)*rad)
     phm(1:nang,1)=3.d0/(16.d0*pi)*(1.d0+cs(1:nang)**2)
     phm(1:nang,2)=phm(1:nang,1)
     phm(1:nang,3)=3.d0/(8.d0*pi)*cs(1:nang)
     phm(1:nang,4)=phm(1:nang,3)
     phm(1:nang,5)=3.d0/(16.d0*pi)*(cs(1:nang)**2-1.d0)
     phm(1:nang,6)=0.d0

! number of averaging
     nav=10
  endif
! initialization end

  wvn=2.d0*pi/wlc

! get grid points
!! 96.2.29
  nintpr=3
  nintpi=3

  do ir=1,kr
     if(cr3<cr(ir)) then
        ip=ir-1; exit
     endif
  enddo
  do ii=1,ki
     if(ci3<ci(ii)) then
        jp=ii-1; exit
     endif
  enddo

  ipt=ip; jpt=jp
  if(ip >= nr-1) ipt=nr-2
  if(jp >= ni-1) jpt=ni-2
  xr(1:nintpr)=cr(ipt:ipt+nintpr-1)
  xi(1:nintpi)=ci(jpt:jpt+nintpi-1)

! get interpolated optical constants
  cg=0.d0; vl=0.d0
  cext=0.d0; cabs=0.d0

  ph(1:nang,1:6)=0.d0   !! mm220128

  rmin=pr(1,3)
  rmax=pr(1,4)
  rmin1=szp(1)/(wvn*rszp2)
  if(rmin1 > rmin) then
     intvl9=int(log(rmin1/rmin)/del)+1
     r2=rmin1*exp(-dble(intvl9)*del)
  else
     intvl9=0
  endif
  intvl8=intvl+intvl9

  do ir8=1,intvl8
     ir=ir8-intvl9
     if(ir <= 0) then
        indr=1
        r1=r2
        r2=r1*rszp2**2
     else
!        dcg(ir)=0.d0
        indr=2
        r1=szp(ir)/(wvn*rszp2)
        r2=szp(ir)/wvn*rszp2
     endif
     dr=(r2-r1)/dble(nav)
     r3=r1-dr*0.5d0
     v1=0.d0

     do j=1,nav
        r3=r3+dr
        dlr=log(r3/(r3-dr))
        pr(1,1)=r3
        v2=vlspc2(pr)
        cg1=v2*3.d0/4.d0/r3*dlr
        cg=cg+cg1
        vl=vl+v2*dlr
        v1=v1+v2
     enddo
     v1=v1/dble(nav)

! interpolation
     if(indr==1) then
        x=(r1+r2)*0.5d0*wvn
        call smllop(x,cr3,ci3,qext9,qsca9)
        qabs9=qext9-qsca9
        cof=del*3.d0/4.d0/x*v1*wvn
        cext=cext+qext9*cof
        cabs=cabs+qabs9*cof

        ph(1:nang,1:6)=ph(1:nang,1:6)+phm(1:nang,1:6)*cof*qsca9
     else
        do l=1,nang+2
! interpolation to get corresponding q to cr3 and ci3
           xi2(1:nintpi)=log(xi(1:nintpi))

           ipol1=ipol
           if(ipol==4) ipol1=kpol2
           if(l>nang) ipol1=1
           do ip=1,ipol1
              if(ip<=2) then
                 do ix=1,nintpr
                    y(1:nintpi)= log(q(ir,ip,ix+ipt-1,jpt:jpt+nintpi-1,l))
                    pint=log(ci3)
                    zi(ix)= bintp(pint,nintpi,xi2,y)
                 enddo
                 z=bintp(cr3,nintpr,xr,zi)
                 z=exp(z)*v1*wvn
              else
                 do ix=1,nintpr
                    y(1:nintpi)=q(ir, ip, ix+ipt-1, jpt:jpt+nintpi-1,l)
                    pint=log(ci3)
                    zi(ix)= bintp(pint,nintpi,xi2,y)
                 enddo
                 z=bintp(cr3,nintpr,xr,zi)
                 z=z*v1*wvn
              endif

              if(l <= nang) then
                 ph(l,ip)=ph(l,ip)+z
              else if(l==nang+1) then
                 cext=cext+z
              else if(l==nang+2) then
                 cabs=cabs+z
              endif
           enddo
        enddo
     endif
  enddo
  err=' '
  return
end subroutine gtph5du
