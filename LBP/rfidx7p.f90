subroutine rfidx7p(iuk,iup,iukd,ipol,ins0,wlc,rh,mp0,ispcvp,rfracp,asphr, &
     rop,dryap,nawp,awcrp,nv,nwlv,wlv,rfi,crw1,ciw1,&
     voldp,volp,cextp,cabsp,nang,ang,php,err)
! Get parameters for wet particle polydispersions
!--- history
! 95. 2.22 created from getpar
! 95. 9.16 modified
! 96. 2.29  2nd order polynomial interpolation (t.y.nakajima)
! 96. 3.11 allow non-spherical parameters
! 96. 7.19 debugged (takayabu, t. y. nakajima)
! 96. 7.19 log for ci and cr interpolation.(t.y.nakajima)
!     5. 5 awcrp(*,*,2)
! 97. 3.17 crw1,ciw1 on subroutine rfidxb by t.y.nakajima
! 97. 5.6  same debug as
! 08. 9.19 modified Fortran 90 free-style form
! 10. 3.31 modified for Ping-Yang column scattering
! 10. 5.17 Remove duplicated save
! 20.12.31 add input ipol
! 21.12.28 add polarized parameters (Ping-Yang Hex column)
! 22. 1.22 debuged.  
!--- input
! iuk      i          device unit number for reading a kernel file
! iup      i          device unit number for reading a kernel file (PY)
! iukd     i          device unit number for reading a kernel file (Dubovik)
! ipol     i          flag of polarization
! ins      i          flag of nonspherical calculation
! wlc      r          wavelength (cm)
! rh       r          relative humidity
! mp       i          aerosol number
! ispcvp i(3,kptc)    fundamental materials for 3 comp. internal mixture (1-8)
! rfracp r(3,kptc)    dry volume fraction of the dry mixture
! asphr  r(3,kptc)    non-spherical parameters (x0, g, r)
! rop    r(kptc)      paricle density relative to water
! dryap  r(6,4,kptc)  dv/dlnr parameters for the dry mixture
!                     see vlspc2
!                     c-values (coefficients of volume spectrum) are ralative
! nawp   i(kptc)      number of aw
! awcrp r(kaw,kptc,2) 1: aw      water activity (see shettle and fenn.)
!                     2: rmmd    rmmd
! nv        i         number of fundamental species (1-8)
! nwlv      i         number of wavelengths for refractive index tables
! wlv    r(kwl)       wavelengths (micron)
! rfi   r(kwl,knv,2)  refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                     with log-regular wavelength interval
!--- output
! voldp    r          total dry volume
! volp     r          total volume
! cextp    r          extinction cross section
! cabsp    r          absorption cross section
! nang     i          number of scattering angles
! ang    r(nang)      scattering angle in degrees
! php    r(knang)     volume scattering phase function
! err      c*64       error indicater. if ' ' then normal end.
!---
  use paras, only: kptc,kaw,kwlv,knv,knang,kintvl,kpol,kpol2,kintvl_du
  implicit none

! input
  integer,intent(in):: iuk       !! device number
  integer,intent(in):: iup       !! device number
  integer,intent(in):: iukd      !! device number
  integer,intent(in):: ipol      !! number of polarization component
  integer,intent(in):: ins0       !! flag of nonspherical calculation
  real(8),intent(in):: wlc       !! wavenumber [cm]
  real(8),intent(in):: rh        !! relativev humidity
  integer,intent(in):: mp0        !! aerosol number
  integer,intent(in):: ispcvp(3,kptc) !! internal mix of fundamental materials
  real(8),intent(in):: rfracp(3,kptc) !! dry component volume fractions
  real(8),intent(in):: asphr(3,kptc)  !! non-spherical parameters
  real(8),intent(in):: rop(kptc)     !! particle density of dry mixture [g/cm3]
  real(8),intent(in):: dryap(6,4,kptc) !! params to define volume size dist.
  integer,intent(in):: nawp(kptc)    !! number of RH to define particle growth
  real(8),intent(in):: awcrp(kaw,kptc,2) !! growth parameters
  integer,intent(in):: nv            !! number of fundamental species
  integer,intent(in):: nwlv          !! number of wavelenggths for rfi tables
  real(8),intent(in):: wlv(kwlv)     !! wavelengths [micron] for rfi tables
  real(8),intent(in):: rfi(kwlv,knv,2)  !! refractive index of materials

! output
  real(8),intent(out):: crw1     !! relative reflactive index (real part)
  real(8),intent(out):: ciw1     !! relative reflactive index (imaginary part)
  real(8),intent(out):: voldp    !! total dry volume of aerosol
  real(8),intent(out):: volp     !! total volume of aerosol
  real(8),intent(out):: cextp    !! extinction cross section
  real(8),intent(out):: cabsp    !! absorption cross section
  integer,intent(out):: nang     !! number of scaattering angles
  real(8),intent(out):: ang(knang)    !! scattering angles [degree]
  real(8),intent(out):: php(knang,kpol2) !! volume scattering phase function
  character,intent(inout):: err*64    !! error code

! for work (aerosol parameters)
  integer:: nmode                !! number of aerosol mode
! for rawb and wet grouth of aerosol
  real(8):: rmd                  !! dry aerosol radius relative to water
  real(8):: rmw                  !! wet aerosol radius relative to water
  integer:: itp
  real(8):: aw(kaw), rmmd(kaw)   !! water activity
  real(8):: rawb                 !! mode radius of wet aerosols [cm]

! for wavelength
  real(8):: wl1,dwl,pint
  integer:: iwl,iwlt

! for refractive index
  real(8):: cr0,ci0,cr01,ci01,cr1,ci1
  integer:: nintp
  real(8):: x(3),yr(3),yi(3)
  real(8):: bintp

! for gtph7
  real(8):: pr(6,4)              !! parameter packet
  real(8):: x0                   !! critical size parameter for mie scaattering
  real(8):: gg                   !! asymmetry parameter for transmittetd ray
  real(8):: rp                   !! surface area fraction to that of sphere

  real(8):: ph(knang,kpol2)      !! volume phase function
  real(8):: ph0(knang,kpol)      !! volume phase function
  real(8):: phpy(knang,kpol2)      !! volume phase function
  real(8):: cext                 !! extinction cross section
  real(8):: cabs                 !! absorption cross section
  real(8):: cg                   !! geometrical cross section
  real(8):: vl                   !! volume
  integer,save:: intvl                !! number of size interval in kernel file
  real(8),save:: szp(kintvl)          !! size parameter

! for gtph4du
  integer,save:: intvl_sp
  real(8),save:: szp_sp(kintvl_du)

! range of wavelength (Ping Yang)
  real(8):: wlmn = 0.225d0  !! minimum wavelength (micron)
  real(8):: wlmx = 99.9d0   !! maximum wavelength (micron)

! range of refractive index (Dubovik Kernel)
  real(8):: crmn = 1.33d0  !! minimum value of rfi (real) 
  real(8):: crmx = 1.60d0  !! maximum value of rfi (real) 
  real(8):: cimn = 5.d-4   !! minimum value of rfi (imaginary)
  real(8):: cimx = 5.d-1   !! maximum value of rfi (imaginary)

! for work
  integer:: ins           !! flag of nonspherical scattering
  integer:: mp            !! aerosol type
  integer:: m,isp,is
!--exec
  err=' '

! changable nonspherical flag  
  ins=ins0
! changable aerosol type
  mp=mp0

! loop for polydispersion
  nmode=int(dryap(1,2,mp))
  pr(1,2)=1.d0
  pr(1,3)=dryap(1,3,mp)

! better to change maxmum radius from 30 micron (to 100 micron)
  pr(1,4)=1.d-2

! These parameters are used in case of ins=3.
  if(ins==3) then
     x0=asphr(1,mp)
     gg=asphr(2,mp)
     rp=asphr(3,mp)
  else
     x0=1.d9
  endif

! loop for size distribution mode
  cextp=0.d0; cabsp=0.d0
  volp =0.d0; voldp=0.d0
  php(1:knang,1:kpol2)=0.d0

  wl1=wlc*1.0d4

  if(ins==2.and.(wl1<wlmn.or.wl1>wlmx)) then
     ins=0; mp=2
  endif

  do m=1,nmode
     pr(2:6,1) =dryap(2:6,m,mp)

     if(nawp(mp)<=0) then
!! dry aerosols
        rmd=1.d0; rmw=1.d0
     else
!! wet aerosols
        itp=int(pr(2,1))
        if(itp/=2) then
           err='wet aerosols assume only log-normal size distribution'
           return
        endif
        rmd=dryap(5,m,mp)
        aw  (1:nawp(mp))=awcrp(1:nawp(mp),mp,1)
        rmmd(1:nawp(mp))=awcrp(1:nawp(mp),mp,2)
        rmw=rawb(rh,rmd,rop(mp),nawp(mp),aw,rmmd)
        pr(5,1)=rmw
        pr(3,1)=pr(3,1)*(rmw/rmd)**3
     endif

! search wavelength in refractive index table
     wl1=wlc*1.0e4
     dwl=log(wl1/wlv(1))/log(wlv(nwlv)/wlv(1))*(nwlv-1)+1
     iwl=int(dwl)
     nintp=3
     if(nintp>nwlv) nintp=nwlv
     iwlt=iwl
     if(iwl>=nwlv-1) iwlt=nwlv-2
     x (1:nintp)=log(wlv(iwlt:iwlt+nintp-1))
        
! water refractive index:
     yr(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,1,1) )
     yi(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,1,2) )
     pint=log(wl1)

     crw1=exp( bintp(pint,nintp,x,yr) )
     ciw1=exp( bintp(pint,nintp,x,yi) )

     if(ins/=2) then
! 3 dry aerosol components refractive index:
        cr0=0.d0; ci0=0.d0
        do isp=1,3
           is=ispcvp(isp,mp)
           if(is==0) exit
        
! interpolation to get corresponding cr1 and ci1 to wl1
           nintp=3
           if(nintp>nwlv) nintp=nwlv
           iwlt=iwl
           if(iwl>=nwlv-1) iwlt=nwlv-2
           yr(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,is,1) )
           yi(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,is,2) )
           pint=log(wl1)
           
           cr01=exp( bintp(pint,nintp,x,yr) )
           ci01=exp( bintp(pint,nintp,x,yi) )
           
! dry aerosol refractive index:
           cr0=cr0+cr01*rfracp(isp,mp)
           ci0=ci0+ci01*rfracp(isp,mp)

! wet aerosol refractive index:
           cr1=crw1+(cr0-crw1)*(rmd/rmw)**3
           ci1=ciw1+(ci0-ciw1)*(rmd/rmw)**3
        enddo

        if(cr0<=0.d0) then
           err='particle species does not exist'
           return
        endif
     endif

! check range
     if(ins==1.and.(cr1<crmn.or.cr1>crmx.or.ci1<cimn.or.ci1>cimx)) then
        ins=0
        mp=mp0-1
     endif

! get optical parameters (Dubovik spheroid)
     if(ins==1) then
        call gtph5du(iukd,wlc,cr1,ci1,pr,nang,ipol,ang,ph, &
             cext,cabs,cg,vl,intvl_sp,szp_sp,err)
        if(err/='') return
        volp = volp +vl
        
! total dry volume
        pr(3,1)=dryap(3,m,mp)
        pr(5,1)=dryap(5,m,mp)
        call getv(wlc,pr,intvl_sp,szp_sp(1:intvl_sp),vl)
        voldp=voldp+vl

! get optical parameters (Ping-Yang Hex-column)
     elseif(ins==2) then

        call gtph7py(iup,wlc,pr,phpy,cext,cabs,cg,vl,nang,ang,err)
        if(err/='') return
        volp = volp +vl
        voldp = volp

! input scattering matrix parameters
! P11:ph(1), P22:ph(2), P33:ph(3), P44:ph(4), P12:ph(5), P34: ph(6)
        ph(1:knang,1)=phpy(1:knang,1)
        ph(1:knang,2)=phpy(1:knang,3)
        ph(1:knang,3)=phpy(1:knang,4)
        ph(1:knang,4)=phpy(1:knang,6)
        ph(1:knang,5)=phpy(1:knang,2)
        ph(1:knang,6)=phpy(1:knang,5)
     else
! get optical parameters
        call gtph5p(iuk,wlc,cr1,ci1,pr,x0,gg,rp,nang,ipol,ang,ph0, &
             cext,cabs,cg,vl,intvl,szp,err)
        if(err/='') return
        volp =volp +vl
        
        ph(1:knang,1)=ph0(1:knang,1)
        !           ph(1:knang,1)=(ph0(1:knang,1)+ph0(1:knang,2))*0.5d0
        ph(1:knang,2)=ph(1:knang,1)
        ph(1:knang,3)=ph0(1:knang,3)
        ph(1:knang,4)=ph0(1:knang,3)
        ph(1:knang,5)=ph0(1:knang,2)-ph0(1:knang,1)
        !           ph(1:knang,5)=(ph0(1:knang,2)-ph0(1:knang,1))*0.5d0
        ph(1:knang,6)=ph0(1:knang,4)
        
! total dry volume
        pr(3,1)=dryap(3,m,mp)
        pr(5,1)=dryap(5,m,mp)
        call getv(wlc,pr,intvl,szp,vl)
        voldp=voldp+vl
     endif
     cextp=cextp+cext
     cabsp=cabsp+cabs
     php(1:knang,1:kpol2)=php(1:knang,1:kpol2)+ph(1:knang,1:kpol2)
  enddo
  return
end subroutine rfidx7p
