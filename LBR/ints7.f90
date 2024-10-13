subroutine ints7(init,iua,ius,matm,co2ppm,nl,alt,prs,tmp,nmol,cng,nptc,ins, &
  cnp,ispcvp,rfracp,asphr,rop,dryap,nawp,awcrp,nv,nwlv,wlv,rfi)
! initilization of star5b
! assume non-spherical parameter input
!--- history
! 95. 9.20  created
! 96. 3.11  kptc=20 for more aerosol models and non-spherical parameters
!           getpr1 included
!     5. 5  drop npoly and get all the polydispersion information
! 01.10.14  changed for fortran90 version
! 08.09.20  modified free-style form
! 12.02.28  include nonspherical flag
! 22. 6.26  debuged.
!--- input
! iua        i        device number for atmospheric model file
! ius        i        d!evice number for particle model file
! matm       i        model number of atmosphere (afgl)
!                      1: tropical
!                      2: mid-latitude summer    3: mid-latitude winter
!                      4: high-lat summer        5: high-lat winter
!                      6: us standard
!--- output
! nl         i        number of atmospheric layers
! alt     r(knl)      height at interfaces of layers (km), top to bottom <-- should be bottom to top (tr220724)
! prs     r(knl)      pressure    at interfaces of layers (mb)
! tmp     r(knl)      temperature at interfaces of layers (k)
! nmol       i        number of gases
! cng    r(knl,knm0)  gas concentration (ppmv)
! nptc      i         number of particle polydispersions
! ins    i(kptc)      nonspherical flag
! iver   i(kptc)      type of vertical profile
! cnp  r(knl,kver)    dry volume concentration profile (relative unit)
! ispcvp  i(3,kptc)   fundamental materials for 3 comp. internal mixture (1-8)
! rfracp  r(3,kptc)   dry volume fraction of the dry mixture
! asphr r(3,kptc)     non-spherical parameters (x0, g, r)
! rop    r(kptc)      paricle density relative to water
! dryap r(6,4,kptc)   dv/dlnr parameters for the dry mixture
!                     see vlspc2
!                     c-values (coefficients of volume spectrum) are ralative
! nawp    i(kptc)      number of aw
! awcrp  r(kaw,kptc,2) 1: aw      water activity (see shettle and fenn.)
!                     2: rmmd    rmmd
! nv        i         number of fundamental species (1-8)
! nwlv      i         number of wavelengths for refractive index tables
! wlv    r(kwl)       wavelengths (micron)
! rfi   r(kwl,knv,2)  refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                     with log-regular wavelength interval
!---Declarations
  use paras,only: kver,kptc,knl,knm0,kaw,kwlv,knv
  implicit none

! input
  integer,intent(inout):: init
  integer,intent(in):: iua,ius   !! device number
  integer,intent(in):: matm      !! atmospheric number
  real(8),intent(in):: co2ppm
  
! output
  integer,intent(inout):: nl       !! number of atmospheric layer
  real(8),intent(out):: alt(knl) !! altitude [km]
  real(8),intent(out):: prs(knl) !! pressure [mb]
  real(8),intent(out):: tmp(knl) !! temperature [K]
  integer,intent(out):: nmol          !! number of gases
  real(8),intent(inout):: cng(knl,knm0) !! gas concentration [ppmv]

  integer,intent(inout):: nptc          !! number of particle polydispersions
  integer,intent(inout):: ins(kptc)     !! nonspherical flag
  real(8),intent(inout):: cnp(knl,kptc) !! dry volume concentration
  integer,intent(inout):: ispcvp(3,kptc) !! internal mixture of fundamental
  real(8),intent(inout):: rfracp(3,kptc) !! dry component volume fractions
  real(8),intent(inout):: asphr(3,kptc)  !! non-spherical parameters
  real(8),intent(inout):: rop(kptc)     !! particle density of dry mixture [g/cm3]
  real(8),intent(inout):: dryap(6,4,kptc) !! param to define volume size dist.
  integer,intent(inout):: nawp(kptc)    !! number of RH to define particle growth
  real(8),intent(inout):: awcrp(kaw,kptc,2) !! growth parameters

  integer,intent(inout):: nv            !! number of fundamental species
  integer,intent(inout):: nwlv          !! number of wavelenggths for rfi tables
  real(8),intent(inout):: wlv(kwlv)     !! wavelengths [micron] for rfi tables
  real(8),intent(inout):: rfi(kwlv,knv,2)  !! refractive index of fundamental

! for mlatm
  integer,parameter:: knm1=7     !! number of major gas
  integer,parameter:: knm2=21    !! number of minor gas
  integer,parameter:: katm=6     !! number of standard atmospheric condition
  integer,parameter:: knm=knm1+knm2
  character(8)::idm(knm)
  integer::idms(knm,10)
  real(8)::wmol(knm,10),rams(knm,10)
  real(8),save::pmatm(knl,katm),tmatm(knl,katm),amol(knl,knm1,katm)
  real(8)::dnsty(knl,katm),trac(knl,knm2)

! for getpar
  integer::natm,nm1,nm2,iptc
  integer:: iver(kptc)
  real(8)::airm,ratio,cnp1(knl,kver)
!---
  if(init > 0) then
     init=0
! model atmospheres
     call mlatm(iua,nl,natm,nm1,nm2,idm,idms,wmol,rams,airm, &
          alt,pmatm,tmatm,dnsty,amol,trac)

!  get aerosol parameters
     call gtpar7(ius,nl,nptc,ins,iver,cnp1,ispcvp,rfracp,asphr,rop, &
          dryap,nawp, awcrp,nv,nwlv,wlv,rfi)
  endif

  prs(1:knl)=pmatm(1:knl,matm)
  tmp(1:knl)=tmatm(1:knl,matm)
  cng(1:knl,1:knm1)=amol(1:knl,1:knm1,matm)
  cng(:,8:28)=trac(:,1:21)

  ratio=co2ppm/cng(1,2)
  cng(1:knl,2)=cng(1:knl,2)*ratio

  nmol=knm1

! set vertical profiles of polydispersion
  do iptc=1,nptc 
     cnp(1:nl,iptc)=cnp1(1:nl,iver(iptc))
  enddo
  
  return
end subroutine ints7

subroutine mlatm(iu,nl,natm,nm1,nm2,idm,idms,wmol,rams,airm, &
     alt,pmatm,tmatm,dnsty,amol,trac)
! set model profiles
!--- history
! 90. 2.24 modified by t. nakajima
! 92. 3.30 idms*8 -> idms for sun-sparc read-statement trouble
! 94. 3.29 idms(i,j)=' ' -> 0
! 94. 6.20 add return
!    10.30 9 format(a8) for ibm
! 95. 9.15 rewind iu
! 01.10.14  changed for fortran90 version
! 08. 9.20  modified free-style form
!--- input
! iu        i         read unit number
!--- output
! nl        i         number of layers = 50
! natm      i         number of model atmospheres = 6
! nm1       i         number of molecules of first kind = 7
! nm2       i         number of molecules of second kind = 21
! idm     c(knm)*8    molecular id
! idms   i(knm,10)    isotope id (0 means end)
! wmol   r(knm,10)    molecular weights
! rams   r(knm,10)    relative abundance of isotopes (0 means end)
! airm      r         molecular weight of air
! alt     r(knl)      altitude (k)
! pmatm  r(knl,katm)  model pressure profiles (mb)
!                     1: tropical,
!                     2: mid-latitude summer,  3: mid-latitude winter
!                     4: high-lat summer    ,  5: high-lat winter
!                     6: us standard
! tmatm  r(knl,katm)  model temperature profiles
! dnsty  r(knl,katm)  density (air molecules / cm3) = p*na/r/t
! amol   r(knl,knm1   molecular profiles (ppmv)
!            ,katm)   1: h2o    2: co2     3: o3      4: n2o     5: co
!                     6: ch4    7: o2
! trac   r(knl,knm2)  trace gase frofiles (ppmv)
!                     8: no     9: so2    10: no2    11: nh3    12: hno3
!                    13: oh    14: hf     15: hcl    16: hbr    17: hi
!                    18: clo   19: ocs    20: h2co   21: hocl   22: n2
!                    23: hcn   24: ch3cl  25: h2o2   26: c2h2   27: c2h6
!                    28: ph3
!
! err      c*64       error indicater
!
!--- data file
! MLATMD
!--
  use paras , only : knl

! parameters in this routine
  integer,parameter:: knm1=7     !! number of major gas
  integer,parameter:: knm2=21    !! number of minor gas
  integer,parameter:: katm=6     !! number of standard atmospheric condition
  integer,parameter:: knm=knm1+knm2

! input
  integer,intent(in)::iu

! output
  character(8),intent(out)::idm(knm)
  integer,intent(out)::idms(knm,10)
  real(8),intent(out)::wmol(knm,10),rams(knm,10)
  real(8),intent(out)::pmatm(knl,katm),tmatm(knl,katm)
  real(8),intent(out)::amol(knl,knm1,katm),trac(knl,knm2)
  real(8),intent(out)::alt(knl),dnsty(knl,katm)
!--- work area
  character ch*1
  integer::nm1,nm2,natm,nm,i,j,k,ns,nl
  real(8)::airm
!--exec
! read numbers
  read(iu,'(a1)') ch
  read(iu,*) nm1,nm2,natm,nl
  nm=nm1+nm2

! read molecular information
  read(iu,'(a1)') ch
  read(iu,*) airm
  read(iu,'(a1)') ch
  idms(:,:)=0
  rams(:,:)=0.0d0
  do i=1,nm
     read(iu,'(a8)') idm(i)
     read(iu,*) ns
     read(iu,*) idms(i,1:ns)
     read(iu,*) wmol(i,1:ns)
     read(iu,*) rams(i,1:ns)
  enddo

! read altitude
  read(iu,'(a1)') ch
  read(iu,*) alt(1:nl)

! read pressure
  do j=1,natm
     read(iu,'(a1)') ch
     read(iu,*) pmatm(1:nl,j)
  enddo

! read temperature
  do j=1,natm
     read(iu,'(a1)') ch
     read(iu,*) tmatm(1:nl,j)
  enddo

! read molecular profiles
  do k=1,nm1
     do j=1,natm
        read(iu,'(a1)') ch
        read(iu,*) amol(1:nl,k,j)
     enddo
  enddo

  do j=1,natm
     read(iu,'(a1)') ch
     read(iu,*) dnsty(1:nl,j)
  enddo

! read trace gase profiles
  do  k=1,nm2
     read(iu,'(a1)') ch
     read(iu,*) trac(1:nl,k)
  enddo
  return
end subroutine mlatm
