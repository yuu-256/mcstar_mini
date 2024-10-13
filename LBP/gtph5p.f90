subroutine gtph5p(iuk,wlc,cr3,ci3,pr,x0,gg,rp,nang,ipol,ang,ph, &
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
! 22. 2.21  fixed bug.
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
  use paras, only: knang, kintvl,kpol,pi,rad
  implicit none

! for input
  integer,intent(in):: iuk       !! device number
  real(8),intent(in):: wlc       !! wavelength [cm]
  real(8),intent(in):: cr3       !! refractive index (real part)
  real(8),intent(in):: ci3       !! refractive index (imaginary part)
  real(8),intent(inout):: pr(6,4)   !! parameter packet
  real(8),intent(in):: x0        !! critical size parameter for mie scaattering
  real(8),intent(in):: gg        !! asymmetry parameter for transmittetd ray
  real(8),intent(in):: rp        !! surface area fraction to that of sphere

! for output
  integer,intent(inout):: nang     !! number of scattering angles
  integer,intent(in):: ipol     !! number of polarization component
  real(8),intent(inout):: ang(knang)   !! scattering angles [degree]
  real(8),intent(out):: ph(knang,kpol) !! volume phase function
  real(8),intent(out):: cext         !! extinction cross section
  real(8),intent(out):: cabs         !! absorption cross section
  real(8),intent(out):: cg           !! geometrical cross section
  real(8),intent(out):: vl           !! volume
  integer,intent(inout):: intvl        !! number of size interval in kernel file
  real(8),intent(inout):: szp(kintvl)  !! size parameter
  character,intent(inout):: err*64

! initialization
  integer:: init=1
  integer,parameter:: npol=4
  real(8),save:: del,rszp2,phm(knang,npol)

! kernel parameters
  integer,parameter:: krf=200        !! number of refractive index
  integer,parameter:: kr=20          !! number of real part in each section
  integer,parameter:: ki=20          !! number of imaginary part in each sec
  integer,parameter:: ksec=3         !! number of section
  integer,parameter:: knang2=knang+2

  integer:: k,i,j,l,ip,jp,kp,ipt,jpt,ir,ix,iy
  integer:: ierr,ipol0
  integer:: nrf
  character:: ch*1
  real(8),save:: cr(krf),ci(krf)
  real(8),save:: q(kintvl,kpol,krf,knang2)
  real(8),save:: dcr,crt(kr,ksec),cit(ki,ksec)
  integer,save:: ict(kr,ki,ksec)
  integer,save:: nav
  integer:: ns,nt

  integer:: nintpr,nintpi,irf(3,3)
  real(8):: xr(3),xi(3)

  integer:: intvl9,intvl8,ir8
  integer:: indr
  real(8):: r1,r2,r3,dr,dlr,v1,v2,cg1,rmin1
  real(8):: vlspc2

  real(8):: cof
  real(8):: xi2(3),y(3),zi(3),z,dci
  integer:: ic,irfi
  real(8):: bintp,pint

! for smllop 
  real(8):: x,qext9,qsca9,qabs9

! for nons
  integer:: ix0
  real(8):: cr1,ci1
  real(8):: dcgl,dcgs
  real(8):: rmin,rmax
  real(8):: wvn

  real(8):: dabs(kintvl)
  real(8):: dcg(kintvl)
  real(8):: dext(kintvl)
  real(8):: phs(knang)

! kernel data parameters
  integer:: nr(1:ksec)  =(/    3,     9,    18/)
  integer:: ni(1:ksec)  =(/   16,    10,     4/)
  real(8):: crmn(1:ksec)=(/1.3d0, 1.0d0, 1.0d0/)
  real(8):: crmx(1:ksec)=(/1.5d0, 1.8d0, 2.7d0/)
  real(8):: cimn(1:ksec)=(/1.d-9, 1.d-4, 1.d-1/)
  real(8):: cimx(1:ksec)=(/1.d-4, 1.d-1, 1.0d0/)

!--exec
  if(init>0) then
     init=0
     rewind iuk
! get kernel ang, szp, ak, cext, cabs
     read(iuk,*) intvl, nang, ipol0
     if (nang > knang) then
        err = 'nang > knang'; return
     end if

     if (ipol0 > kpol) then
        err = 'ipol > kpol' ; return
     end if
     if(ipol0/=4) then
        err='assume only ipol=4 case for PKRNL.OUT'; return
     endif
     
     read(iuk,*) szp(1:intvl)
     read(iuk,*) ang(1:nang)
     del=log(szp(intvl)/szp(1))/dble(intvl-1)
     rszp2=exp(del*0.5d0)

     phm(1:nang,1)=3.d0/8.d0/pi
     phm(1:nang,2)=3.d0/8.d0/pi*cos(ang(1:nang)*rad)**2
     phm(1:nang,3)=3.d0/8.d0/pi*cos(ang(1:nang)*rad)
     phm(1:nang,4)=0.d0

!! read kernel
     do ir=1,krf
        read(iuk,'(a1)',iostat=ierr) ch
        if(ierr/=0) then
           nrf=ir-1; exit
        endif
        read(iuk,*) cr(ir),ci(ir)
        ci(ir)=abs(ci(ir))
        do i = 1, intvl
           do j = 1, ipol0
              read(iuk,*) q(i,j,ir,1:nang)
           enddo
        enddo
        read(iuk,*) q(1:intvl,1,ir,nang+1:nang+2)
     enddo

!! get refractive index tables (3 tables)
     do k=1,ksec
        dcr=(crmx(k)-crmn(k))/dble(nr(k)-1)
        do i=1,nr(k)
           cr1=crmn(k)+dcr*dble(i-1)
           if(abs(cr1-1.0) <= 0.d0) cr1=1.01
           crt(i,k)=cr1
        enddo

        dci=log(cimx(k)/cimn(k))/dble(ni(k)-1)
        do i=1,ni(k)
           ci1=dci*dble(i-1)
           cit(i,k)=cimn(k)*exp(ci1)
        enddo
     enddo

!! find grids for kernel refractive indices
     do k=1,ksec
        ict(1:nr(k),1:ni(k),k)=0
     enddo

     do ir=1,nrf
        cr1=cr(ir); ci1=ci(ir)
        do k=1,ksec
           do i=1,nr(k)
              if(abs(crt(i,k)/cr1-1.d0) <= 0.01d0) exit
           enddo
           do j=1,ni(k)
              if(abs(cit(j,k)/ci1-1.d0) <= 0.01d0) exit
           enddo
           ict(i,j,k)=ir
        enddo
     enddo

     do k=1,ksec; do i=1,nr(k); do j=1,ni(k)
        if(ict(i,j,k)==0) then
           err='no matching grid'; return
        endif
     enddo; enddo; enddo

! number of averaging
     nav=10
  endif
! initialization end

  wvn=2.d0*pi/wlc
! write(*,*) wlc
! find the domain
  do kp=1,ksec
     if(ci3 <= cimx(kp).or.kp==ksec) then
        ns=nr(kp); nt=ni(kp)
        exit
     endif
  enddo

  do ip=1,ns-1
     if(cr3 <= crt(ip,kp)) exit
  enddo
  do jp=1,nt-1
     if(ci3.le.cit(jp,kp)) exit
  enddo

  ip=max(1, ip-1)
  jp=max(1, jp-1)

! get grid points
!! 96.2.29
  nintpr=3
  nintpi=3
  if(nintpr > nr(kp)) nintpr=nr(kp)
  if(nintpi > ni(kp)) nintpi=ni(kp)

  ipt=ip
  jpt=jp
  if(ip >= nr(kp)-1) ipt=nr(kp)-2
  if(jp >= ni(kp)-1) jpt=ni(kp)-2
  xr(1:nintpr)=cr(ict(ipt:ipt+nintpr-1,jpt             ,kp))
  xi(1:nintpi)=ci(ict(ipt             ,jpt:jpt+nintpi-1,kp))

  irf(1:nintpr,1:nintpi)=ict(ipt:ipt+nintpr-1,jpt:jpt+nintpi-1,kp)

! get interpolated optical constants
  cg=0.d0; vl=0.d0
  cext=0.d0; cabs=0.d0

  dcgl=0.d0; dcgs=0.d0
  do ix0=1,intvl
     if(szp(ix0) > x0) exit
  enddo
  phs(1:nang)=0.d0; ph(1:nang,1:ipol)=0.d0

  rmin=pr(1,3)
  rmax=pr(1,4)
  rmin1=szp(1)/wvn/rszp2
  if(rmin1 > rmin) then
     intvl9=int(log(rmin1/rmin)/del)+1
     r2=rmin1*exp(-intvl9*del)
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
        dcg(ir)=0.d0
        indr=2
        r1=szp(ir)/wvn/rszp2
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

        if(ir <= ix0)then
           dcgs=dcgs+cg1
        else
           dcgl=dcgl+cg1
        endif

        if(indr <= 1) then
           dcg(1)=dcg(1)+cg1
        else
           dcg(ir)=dcg(ir)+cg1
        endif
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
!        ph(1:nang,1:ipol)=ph(1:nang,1:ipol)+phm(1:nang,1:ipol)*cof*qsca9
        ph(1:nang,1)=ph(1:nang,1)+(phm(1:nang,2)+phm(1:nang,1))*0.5d0*cof*qsca9
        if(ipol/=1) then
!!           ph(1:nang,2)=ph(1:nang,2)+(phm(1:nang,2)-phm(1:nang,1))*0.5d0*cof*qsca9
!           ph(1:nang,2)=ph(1:nang,2)+phm(1:nang,2)*cof*qsca9
!ms230625    ph(1:nang,3:npol)=ph(1:nang,3:npol)+phm(1:nang,3:npol)*cof*qsca9
           ph(1:nang,2:npol)=ph(1:nang,2:npol)+phm(1:nang,2:npol)*cof*qsca9           
        endif
     else
        do l=1,nang+2
! interpolation to get corresponding q to cr3 and ci3
           xi2(1:nintpi)=log(xi(1:nintpi))

! I (ip=1)
           do ix=1,nintpr
              ic=0
              do iy=1,nintpi
                 ic=ic+1
                 irfi=irf(ix,iy)
                 if(l<=nang) then
                   y(ic)=log((q(ir,2,irfi,l)+q(ir,1,irfi,l))*0.5d0)
                 else
                   y(ic)=log(q(ir,1,irfi,l))
                 endif
              enddo
              pint=log(ci3)
              zi(ix)= bintp(pint,nintpi,xi2,y)
           enddo
           z=bintp(cr3,nintpr,xr,zi)
           z=exp(z)*v1*wvn
           if(l <= nang) then
              ph(l,1)=ph(l,1)+z
           else if(l==nang+1) then
              cext=cext+z
!write(*,*) ir8,z,cext
              dext(ir)=z
              cycle
           else if(l==nang+2) then
              cabs=cabs+z
              dabs(ir)=z
              exit
           endif
           if(ipol==1) cycle
! Q (ip=2)
           do ix=1,nintpr
              ic=0
              do iy=1,nintpi
                 ic=ic+1
                 irfi=irf(ix,iy)
                 y(ic)=log(q(ir,2,irfi,l))
              enddo
              pint=log(ci3)
              zi(ix)= bintp(pint,nintpi,xi2,y)
           enddo
           z=bintp(cr3,nintpr,xr,zi)
           z=exp(z)*v1*wvn
           if(l<=nang) then
              ph(l,2)=ph(l,2)+z
           endif
! U,V (ip=3,4)
           do ip=3,ipol
              do ix=1,nintpr
                 ic=0
                 do iy=1,nintpi
                    ic=ic+1
                    irfi=irf(ix,iy)
                    y(ic)=q(ir, ip, irfi, l)
                 enddo
                 zi(ix)= bintp(ci3,nintpi,xi,y)
              enddo
              z=bintp(cr3,nintpr,xr,zi)
              z=z*v1*wvn
              ph(l,ip)=ph(l,ip)+z
           enddo
        enddo
     endif
  enddo
  err=' '
  return
end subroutine gtph5p
