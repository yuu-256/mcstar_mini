!mm220221
subroutine mcst2(init0,icont,iuk,iup,iukd,iug,indg,inda,isol,ipol,imthd,nda,&
     fis,na0,th0,na1,th1,nfi,fi,npv,llu,llr,rpvw,ngas,wl, &
     icn,wlcn,npoly,ncomp,mptc,vptc,cnpt, &
     nl,alt,prs,tmp,gtmp,cng,ins,cnp,ispcvp,rfracp,asphr, &
     rop,dryap,nawp,awcrp,nwlv,wlv,rfi,nln,ipbf, &
     icycl,dseed,np0,nb,bx,cconc,ag,jatm,fsol,tauav,tausav, &
     taurv,tckdav,gav,flxs1,flxs0,flxe1,flxuxy1,flxd001,ais1,aie1,flxsp1,flxep1,nsos,err)

!-- Monte-Carlo radiative transfer program
!  assuming Lambert surface
! this routine is for broad-band calculation with use mstrn-ckd with abgasx routine
!-- history
!2019. 9. 4 subcal2 (M. Okata)
!tr211026 debugged, added itray, itkd
!tr211104 dbxsf3 with knng,mcstr1
!tr211120 wlcn should be in cm unit and given by the main program, not needed to define in this routine
!tr211123 icont, gtmp
!2021.12.6 shipped (T. Nakajima)
!mm211218 add tms with tmsbxb1
!mm211218 srdbxb1 => srdbxb2 to calc. u1t
!mm211218 debegged in lseg = 3
!mm211220 tmsbxb1 => tmsbxb2
!mm211220 add standard error calculation
!mm211223 add p2ims for monte carlo system
!mm220131 add 3d-ims for monte carlo system
!mm231002 add errors in ais1, aie1, flxsp1, flxep1: 1: values, 2: random errors, 3: truncated error
!-- input variables
! init0            give 1 to itialize rstr6cms
! icont            (1) : on/off Rayleigh (1/0)
!                  (2) : on/off CKD gas absorption (1/0)
!                  (3) : on/off 0th Fourier L(0) output (0/1)
!                  (4) : on/off photon info output to outw-file (0/1)
! iuk              input unit for Mie kernel
! iug              input CKD gas absorption table
! indg             -1: no ground surface
!                   0: lambert surface
!                   1: ocean surface initialization
!                   2: ocean surface with no initialization
!                     when indg>0 and imthd>0 then cabspsingle scattering correction
!                     for ocean surface reflection
!                   3: flat ocean surface tr201216
! inda              0: flux; 1: radiance (not usesd)
! isol              0: night
!                   1: day
! imthd            This routine supports only -1 and 2 in
!                  the rstar definition as follows
!                  -1: nt,  0: dms-method  for intensity/flux
!                   1: ms,  2:tms,  3:ims-method for intensity.
!                     when indg>0 and imthd>0 then single scattering correction
!                     for ocean surface reflection
! fis              solar azimuthal angle relative to x-axis
! nda              nunber of nadir-quadrature angles in the hemisphere
! na0              number of solar incidences
! th0(na0)         solar zenith angles 0-90 degrees
! na1              number of emergent nadir angles in sphere
! th1(kna1u)       emergent nadir angles
!                    0-90 for downward; 90-180 for upward
!  (On the other hand, zenith angle is adopted inside this MC program adopts )
!                   90-180 for downward; 0-90 for upward
! nfi              number of azimuthal angles
! fi(knfi)         azimuthal angles in degrees relative to fis
! npv              number of user-defined viewing (receiver) positions for lseg=1, 2 (backward)
! llu(npv)         side numbers of viewing positions
!                   1 xy+z (toa), 2 xy-z (boa), 3 yz+x, 4 yz-x, 5 zx+y, 6 zx-y
! llr              0 prescribed rpvw by  user; 1 dondom positioning
! rpvw(kpv,3)      user-derined viewing positions when llr>0
! ngas(1)             wavelength number wl for broadband, number of quadrature for narrowband
! ngas(2)             dintvl for narrowband
! wl          CKD band wavelengths
! icn              0: cnpt is the total column dry volume for each polydispersion
!                  1: cnpt is the total column volume for each polydispersion
!                  2: cnpt is optical thickness at wlcn
!                  3: cnpt is optical thickness at each wavelength.
! wlcn             scaling wavelength for icn=2,3
! npoly            number of external mixtures
! ncomp(npoly)     number of aerosol types for the external mixture
! mptc(kpoly,4)    aerosol types for the external mixture
! vptc(kpoly,4)    dry volume ratio for the mixture
! cnpt(npoly)      total volume (icn=0,1) or optical thickness (icn=2,3)
!                    for each polydispersion.
! nl               number of atmospheric layers to define model atmosphere
! alt(knl)         height at interfaces of layers (km), top to bottom
! prs(knl)         pressure    at interfaces of layers (mb)
! tmp(knl)         temperature at interfaces of layers (k)
! gtmp             ground temperature defined by user, standard set is tmp(1)
! cng(knl,knm0)    gas concentration (ppmv)
! cnp(knl,kptc)    dry volume concentration profile (relative unit)
! ispcvp(3,kptc)   fundamental materials for 3 comp. internal mixture (1-8)
! rfracp(3,kptc)   dry volume fraction of the dry mixture
! asphr(3,kptc)    non-spherical parameters (x0, g, r)
! rop(kptc)        paricle density relative to water
! dryap(6,4,kptc)  dv/dlnr parameters for the dry mixture
!                   see vlspc2
!                   c-values (coefficients of volume spectrum) are ralative
! nawp(kptc)       number of aw
! awcrp(kaw,kptc,2) 1: aw      water activity (see shettle and fenn.)
!                   2: rmmd    rmmd
! nwlv             number of wavelengths for refractive index tables
! wlv(kwl)         wavelengths (micron)
! rfi(kwl,knv,2)   refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                    with log-regular wavelength interval
! nln              number of atmospheric layers for radiative transfer.
! ipbf(nln1)       interface number to define sublayers to construct
!                    the transfer atmosphere in mlatm.
!                    ipbf(1)=50, ipbf(nln+1)=1
!                    if >0 then field calculations
!                       <0 then no field calculations
! icycl            1 for cyclic x-y boudary
! dseed            Initial seed for random generator
! np0              Number of photons for forward solar routine (lseg=3)
! nb(3)            box numbers along x, y, z axes
! bx(kb1,3)        box boundary coordinates along x, y, z axes
! cconc(kx,ky,kz,kpoly)
!                  relartive concentrations of boxes relative to userdefined atmosphere for each
!                     aerosol type
! ag(kx,ky)        user defined ground albedo distribution along x-y axes

!-- output variables
! init0     0
! fsol             extraterrestrial solar irradiance (W/m2)
! tauav(nz)        area-averaged optical thickness of an atmospheric sub-layer
!                   defined by ipbf. (1) at bottom layer and (nz=nb(3)) at top layer
! tausav(nz)       same as tauav but for scattering optical thickness
! taurz(nz)        same as tauav but for Rayleigh scattering optical thickness
! tckdav(nz)       same as tauav but for gas absorption optical thickness by CKD approximation
! gav(klgt1,kz)    area-averaged Legendre polynomial coefficients up to klgt1-th order
! flxs1(kna0,6,3)  CKD-averaged spectral solar fluxes  per horizontal area by forward MC (lseg=3)
!                  for solar zenith angles th0, and six system boundaries (ll) and directions (idu)
!                  ll: 1 xy+z (toa), 2 xy-z (boa), 3 yz+x, 4 yz-x, 5 zx+y, 6 zx-y
!                  idu: downward direct, downward diffuse, upward
! flxs0(kna0,6)    Solar insolation fluxes at 6 illuminating surface per horizontal area
! flxe1(6,3)      Same as flxs1, but for thermal emission fluxes per horizontal area by backward MC (lseg=4)
! flxuxy1(kna0,kx,ky)  Boa horizontal upward solar fluxes
! flxd001(kna0)    Donward direct solar flux at boa by plane-parallel approximation
! ais1(kna1,kna0,knfi,kpv,2) CKD-averaged spectral solar diffuse radiances by backward solar MC (lseg=1)
!                   at user-defined viewring positions given by rpvw and for downward/upward
! aie1(kna1u,knfi,kpv,2) Same as ais1, but for thermal emission radiance by lseg=2 (backward thermal)
! nsos(ksos)       number of photons with multiple scattering scattering (1 for no scattering, 2 for
!                     single scattering) in the forward solar MC calculation (lseg=3)
! err*64           error code, ' ' if normal termination

!$end
  use paras

  implicit none

! control parameter of the programs
! 1: zero for no-Rayleigh (itray),2: zero for no-ckd (itkd),3: zero for L(0) radiance (immax)
  integer:: icont(10)

!c file name
  integer,intent(in) :: iuk,iug
  integer,intent(in) :: iup,iukd !mm220221 from rp
  character(len=64) :: err

  integer :: nch,isol,inda,indg(kx,ky),nfi,ngas(2),nda,npoly,icn,nwlv,llu(kpv),llr
  integer :: nww,init0,init,init1,init2,init3,iww,ich,nln,na0,na1,npv
  real(8) :: wlcn,fsol,wgt(kch),ang(knang)

!c rstar parameter
  real(8) ::alt(knl),prs(knl,katm),tmp(knl,katm),cng(knl,knm0,katm), &
       cnp(knl,kptc),rfracp(3,kptc), &
       rop(kptc),dryap(6,4,kptc), &
       awcrp(kaw,kptc,2),asphr(3,kptc),wlv(kwlv),&
       rfi(kwlv,knv,2),ai00(kna1u,kna0,knfi,kntau)
  integer :: ispcvp(3,kptc),nawp(kptc)
  real(8),intent(inout)  :: th0(kna0),th1(kna1u),fi(knfi),wl
  real(8) :: flxd0(kna0,kntau),flxd00(kna0,kntau),flxu0(kna0,kntau),  &
    flxs0(kna0,6),flxs1(kna0,6,3),flxe(kpv,3),flxe1(6,3),ais(kna1u,kna0,knfi,kpv,2,kpol,3), &
    ais1(kna1u,kna0,knfi,kpv,2,kpol,3),aie(kna1u,knfi,kpv,2,kpol,3),aie1(kna1u,knfi,kpv,2,kpol,3), &
    thk0(knln,10),cnpt(kpoly),vptc(kpoly,4),rpvw(kpv,2)
    
  integer :: ipbf(knln1),ncomp(kpoly),mptc(kpoly,4)

!c monte carlo parameter
  real(8) :: drs(3),dvw(3),pvw(3),bx(kb1,3),ag(kx,ky),ps(3),dbx(kb1,3,2)

  integer :: nb(3),mlim(2),mlim1(2)
!c file name
!c data patrameter
  integer :: icycl,np,np0
!c work area
  integer :: nang,l
  real(8) :: sths
  integer :: li,lj,lk,i,j,k,lp
  real(8) :: u1(kpol)
  real(8) :: fis,angr(na4),sum1
  integer :: ind1,ipoly
!tr  type1-3 variables
  real(8) :: phsfr(knang,knln,kpol2), dpf(knln)
  real(8) :: phsfrz(knang,kz,katm,kpol2), dpfz(kz,katm)

!tr210901
integer:: iw99,np1,nl,k1,initd,mxlgn1,mxlgt1,m,nx,ny,nz, &
  lim1,nns,nng

real(8):: wgt1,dz
!mm211218: p2 => p2t
real(8):: ssaa(knln),tkd(knln,kch,kww),taur(knln),taua(knln), &
  ssapoly(knln,kpoly),epoly(knln,kpoly),ppoly(knang,knln,kpoly,kpol2), &
  altz(kz),pz(kz),tz(kz),cconc(kx,ky,kz,kpoly),  &
  ssaaz(kz,katm),tkdz(kz,kch,kww,katm),taurz(kz,katm),tauaz(kz,katm), &
  ssapolyz(kz,kpoly,katm),epolyz(kz,kpoly,katm),ppolyz(knang,kz,kpoly,katm,kpol2), &
  da,a1,p1,pd1,p2t(knang,kpol2),cs2(knang),cs(na4),ff(kz,kpoly,katm), &
  w1,f1,tau1,taut1,taust1,taut3,taust3,cext(kx,ky,kz),ssat(kx,ky,kz),  &
  tautp(kx,ky,kz),taustp(kx,ky,kz),g(1:klgn1,kpol2),gt(1:klgn1,kpol2),g1, &
  ppolyzt(knang,kz,kpoly,katm,kpol2),cs1,  &
  angdp(na4,kx,ky,kz),dx,x1,aintp, &
  wg0,taurv(kz),tauav(kz),tausav(kz),gpoly(klgt1,kz,kpoly,katm), &
  gav(klgt1,kz),taus1,gmol(klgt1,kz,katm),tckdav(kz)
  !mm211218: p3 => p3t
  real(8):: pint4c,plgd,fb(kpol),am0(kna0),p3t(na4,kx,ky,kz,kpol2),am1(kna1)
  real(8):: cplkz(kplk1,kz,katm)  !! coefficient of planck function
  real(8):: b1,b2,bplnk,bgnd(katm),gtmp(katm),u,dseed,rndv(3),galb2
  real(8):: bgnd2,bgnd3,bgndz(kx,ky),cplkzz(kplk1,kx,ky,kz)
  integer:: imthd,na12,lz,ll,ltau,nplk1,lseg,ierr,ibxg(knng),ibyg(knng), &
      lw,nsos(ksos),ii,jj,nfi12,nrnd,nnsg(knng),na0se,ibs(3),idui,ll1, &
      ll2,i1,i2,i3,j1,j2,j3,lpv
  real(8):: fbf1(kna0,kpv,3),fbf(6,3),fdg,fdg0,fug,  &
      flxuxy(kx,ky),flxuxyb(kx,ky),flxuxy1(kna0,kx,ky),c,ug(knng,2),  &
      dsg(knng,3,2),psg(knng,3),flxd001(kna0),ds(3), &
      fbe(3),galb3,earea,earea0
  integer i33(3,3),li2,lli
!mm211214 for tms
  integer :: initt
  real(8) :: u1t(kpol), pow
  real(8) :: taup(kx,ky,kz), tausp(kx,ky,kz), p2(knang,kpol2), p3(knang,kx,ky,kz,kpol2), fp(kx,ky,kz)
  real(8) :: tau3, taus3
!mm211220 for ims
  real(8) :: fims, p3s(na4,kx,ky,kz,kpol2), p32(na4,kx,ky,kz,kpol2), &
  angdps(na4,kx,ky,kz), ssaf(kx,ky,kz)
  integer :: np2
!mm220110 for pims
  real(8) :: pbnd
  integer :: mlimp(2)
  real(8) :: reb, rea
  integer :: kims(2) = (/1,1/)

  character(len=64) ::  erc
  
!mm220221
  integer :: ipol
  integer :: ins(kptc)
  
!mm220416 for point flux module and reflesh
  integer :: lsp, nsp ! delta-m and ims space
  real(8) :: flxep1(kpv,3,3), flxsp1(kna0,kpv,3,3)
  
!mm220417
  integer :: jatm(kx,ky), iatm
  
  !mm231002 for truncated radiance
  real(8) :: utrn(kpol), fbtrn(kpol), fb2(kpol)
  real(8) :: flxep(kpv,3,3), flxsp(kna0,kpv,3,3)
  
  !mm231003 for vector
  integer :: npol
  real(8) :: eps
  
  !mm240206 for satellite sensor
  real(8) :: pb(3)
  
  !mm220221
  call cpcon(ccp)
  eps = ccp(1)*100.0
  
  !mm231114 box parameters
  do i = 1, 3
    dbx(1:nb(i),i,1) = bx(2:nb(i)+1,i) - bx(1:nb(i),i)
    dbx(1:nb(i),i,2) = bx(2:nb(i)+1,i) + bx(1:nb(i),i)
    dbx(1:nb(i),i,1:2) = dbx(1:nb(i),i,1:2) * 0.5d0
  end do
  
  npol = ipol
  if(ipol>1) npol = 6

  do i=1,3
    j1=i-1
    do j=1,3
      j1=j1+1
      if(j1>3) j1=1
      i33(i,j)=j1
    enddo
  enddo

  nx=nb(1); ny=nb(2); nz=nb(3)
  na12=na1/2
  nfi12=nfi/2
  init=init0
  init0=0
  
  flxd001(1:na0)=0.
  
  do iatm = minval(jatm(1:nx,1:ny)),maxval(jatm(1:nx,1:ny))
    if( maxval(jatm(1:nx,1:ny)) /= minval(jatm(1:nx,1:ny)) ) init = 1
    if( count(jatm(1:nx,1:ny)==iatm)==0 ) cycle
    
    call rstr7cms(init,icont,iuk,iup,iukd,iug,indg(1,1),inda,isol,ipol,na0,th0,na1,th1,nfi,  &
      wl,icn,wlcn,npoly,ncomp,mptc,  &
      vptc,cnpt,nl,alt,prs(:,iatm),tmp(:,iatm),cng(:,:,iatm),ins,cnp,ispcvp,rfracp,asphr,rop,dryap, &
      nawp,awcrp,nwlv,wlv,rfi,nln,ipbf,ngas,thk0,fsol,flxd0,flxd00,flxu0,  &
      ai00,nch,nww,wgt,ang,ssaa,tkd,taur,taua,ssapoly,epoly,ppoly,phsfr,dpf,nang,err)

    if(err /= ' ') return

!tr reorder layer-number, bottom to top from top to bottom (rstr),
! other than alt: 1 at ground
! assumed increment of ipbf =1 from top to bottom for rstr
! ppoly: phase function is normalized by Csca for each ipoly
! nz=nln
    do  lz=1,nz+1
      ltau=nz+2-lz
      l=iabs(ipbf(ltau))
      altz(lz)=bx(lz,3)
      pz(lz)=prs(l,iatm)
      tz(lz)=tmp(l,iatm)
    enddo

    lw=0
    do l=1,nln+1
      if(ipbf(l)> 0) lw=lw+1
    enddo
    flxd001(1:na0)=flxd001(1:na0)+flxd00(1:na0,lw)*count(jatm(1:nx,1:ny)==iatm)

! Plank function

    if(wl<2.d0) then
      nplk1=0
     else
      nplk1=2
    endif

    bgnd(iatm)=bplnk(wl,gtmp(iatm))

    if(nplk1>0) then
      b2=bplnk(wl,tz(1))
      do l=1,nz
        b1=b2
        b2=bplnk(wl,tz(l+1))
        cplkz(1,l,iatm)=b1
        cplkz(2,l,iatm)=(b2-b1)/(altz(l+1)-altz(l))
      enddo
    endif
!
    mxlgn1=2*nda+1
    mxlgt1=2*nda
    cs2(1:nang)=cos(ang(1:nang)*rad)

    da=180.0/(na4-1)
    do i=1,na4
      angr(i)=(i-1)*da
    enddo
    cs(1:na4)=cos(angr(1:na4)*rad)

    do k=1,nz
      k1=nln+1-k
      ssaaz(k,iatm)=ssaa(k1)
      tkdz(k,1:nch,1:nww,iatm)=tkd(k1,1:nch,1:nww) !tr assumed nww=1
      taurz(k,iatm)=taur(k1)
      tauaz(k,iatm)=taua(k1)
      ssapolyz(k,1:npoly,iatm)=ssapoly(k1,1:npoly)
      epolyz(k,1:npoly,iatm)=epoly(k1,1:npoly)
      ppolyz(1:nang,k,1:npoly,iatm,1:npol)=ppoly(1:nang,k1,1:npoly,1:npol)
      phsfrz(1:nang,k,iatm,1:npol)=phsfr(1:nang,k1,1:npol)
      dpfz(k,iatm)=dpf(k1)
    enddo

! Legendre expansion and truncation (imthd>=0)
!    m=3
!mm220416
!    call lgndf3(m,nang,cs2,phsfr(:,1),gmol)
!    g1=gmol(1)
!    gmol(1:m)=gmol(1:m)/g1
    do k=1,nz
      g1=2.d0*(1.d0-dpfz(k,iatm))/(2.d0+dpfz(k,iatm))
      gmol(1:3,k,iatm) = (/1.,0.,g1*0.1/)
!memo from rpstar
! greek constants are not used in MC
!      g1=2.d0*(1.d0-dpfz(k))/(2.d0+dpfz(k))
!      gmol(1:3,1,k) = (/1.,0.,g1*0.1/)
!      gmol(3,2,k) = 0.6*g1
!      g1=2.d0*(1.d0-3.d0*dpfz(k))/(2.d0+dpfz(k))
!      gmol(2,4,k) = 0.5*g1
!      gmol(3,5,k) = sqrt(6)*gmol(3,1,k)
    end do

    do ipoly=1,npoly
    do k=1,nz
      p2(1:nang,1:npol)=ppolyz(1:nang,k,ipoly,iatm,1:npol)
!      call lgndf3(mxlgn1,nang,cs2,p2(:,1),g(:,1))
      call phasa2mc(ipol,mxlgn1,nang,ang,p2,g,gt,f1,g1)
      !mm231004 g1=g(1,1)
!tr211026 debugged, cs2 is descending order series, so g becomes negative
      if(abs(g1)<=0) then
        g(1:mxlgn1,1)=0
      else
        !mm231004 g(1:mxlgn1,1)=g(1:mxlgn1,1)/g1
        !mm220224
        if(istar(2) == 1) then
!          call lgndf3(1,nang,cs2,ppolyz(1:nang,k,ipoly,iatm,1),g(1:1,1))
!          g1=g(1,1)
!          g(1,1)=1.d0
          ppolyz(1:nang,k,ipoly,iatm,1:npol) = ppolyz(1:nang,k,ipoly,iatm,1:npol)/4/pi/abs(g1)
        end if
      endif

      !mm211218: mxlgt1 => mxlgn1
      gpoly(1:mxlgn1,k,ipoly,iatm)=g(1:mxlgn1,1)

!truncation
      if(imthd < 0) then
        ff(k,ipoly,iatm)=0
        ppolyzt(1:nang,k,ipoly,iatm,1:npol)=ppolyz(1:nang,k,ipoly,iatm,1:npol)
      else
!        f1=g(mxlgt1+1,1)
!        ff(k,ipoly,iatm)=f1
!        gt(1:mxlgt1,1)=(g(1:mxlgt1,1)-f1)/(1-f1)
!        init1=1
!        do i=1,nang
!          cs1=cs2(i)
!          init1=1
!          sum=0
!          do l=1,mxlgt1
!            sum=sum+(2*l-1)*gt(l,1)*plgd(init1,cs1)
!          enddo
!          p2(i,1)=sum/4/pi
!        enddo
!        ppolyzt(1:nang,k,ipoly,iatm,1:npol)=p2(1:nang,1:npol)
        ff(k,ipoly,iatm)=f1
        ppolyzt(1:nang,k,ipoly,iatm,1:npol)=p2(1:nang,1:npol)
      endif
    enddo; enddo
  end do ! atmosphere
  flxd001(1:na0)=flxd001(1:na0)/(nx*ny)
! averaged ground albedo and bgnd
    galb3=0
    bgnd3=0
    k=0
    do i=1,nx; do j=1,ny
      if(indg(i,j)==0) then
        k = k + 1
        galb3=galb3+ag(i,j)
        bgnd3=bgnd3+log(bgnd(jatm(i,j)))
      end if
    enddo; enddo
    if(k/=0) then
      galb3=galb3/k
      bgnd3=exp(bgnd3/k)
    end if
! construction of 3D optical field
    taurv(1:nz)=0
    tauav(1:nz)=0
    tausav(1:nz)=0
    gav(1:mxlgt1,1:nz)=0

    do i=1,nx
    do j=1,ny
      ! thermal
      bgndz(i,j)=bgnd(jatm(i,j))
      cplkzz(:,i,j,1:nz)=cplkz(:,1:nz,jatm(i,j))
    do k=1,nz
      !mm211218: for tms
      tau3=0
      taus3=0
      p2(1:nang,1:npol)=0
      !end for tms
      taut3=0
      taust3=0
      !mm211218: p2 => p2t
      p2t(1:nang,1:npol)=0
      do ipoly=1,npoly
        w1=ssapolyz(k,ipoly,jatm(i,j))
        f1=ff(k,ipoly,jatm(i,j))
        tau1=epolyz(k,ipoly,jatm(i,j))*cconc(i,j,k,ipoly)
        taus1=w1*tau1
        !mm211218 for tms
        tau3=tau3 + tau1
        taus3=taus3 + taus1
        p2(1:nang,1:npol)=p2(1:nang,1:npol)+taus1*ppolyz(1:nang,k,ipoly,jatm(i,j),1:npol)
        !tauav(k)=tauav(k)+tau1
        !tausav(k)=tausav(k)+taus1
        !end for tms
        gav(1:mxlgt1,k)=gav(1:mxlgt1,k)+gpoly(1:mxlgt1,k,ipoly,jatm(i,j))*taus1
        
        taut1=(1-w1*f1)*tau1
        taust1=(1-f1)*taus1
        taut3 =taut3 +taut1
        taust3=taust3+taust1
        !mm211218: p2 => p2t
        p2t(1:nang,1:npol)=p2t(1:nang,1:npol)+taust1*ppolyzt(1:nang,k,ipoly,jatm(i,j),1:npol)
      enddo
      !mm211218 for tms
      taup(i,j,k)=tau3+taurz(k,jatm(i,j))
      tausp(i,j,k)=taus3+taurz(k,jatm(i,j))
      p2(1:nang,1:npol)=p2(1:nang,1:npol)+phsfrz(1:nang,k,jatm(i,j),1:npol)*taurz(k,jatm(i,j))
      !scaling: important for TMS!!
      if(istar(2)==0) then ! Rstar
        call lgndf3(1,nang,cs2,p2(:,1),g(1:1,1))
        p2(1:nang,1:npol) = p2(1:nang,1:npol)/4/pi/abs(g(1,1))
        g(1,1) = 1.0
      else if(istar(2)==1) then ! Rpstar
        p2(1:nang,1:npol)=p2(1:nang,1:npol)/tausp(i,j,k)
      end if
      !end scaling
      taurv(k) =taurv(k) +taurz(k,jatm(i,j))
      tauav(k) =tauav(k) +taup(i,j,k)
      tausav(k)=tausav(k)+tausp(i,j,k)
      !tauav(k) =tauav(k) +taurz(k)
      !tausav(k)=tausav(k)+taurz(k)
      !end for tms
      gav(1:3,k)=gav(1:3,k)+gmol(1:3,k,jatm(i,j))*taurz(k,jatm(i,j))
      
      tautp(i,j,k)=taut3+taurz(k,jatm(i,j))
      taustp(i,j,k)=taust3+taurz(k,jatm(i,j))
      !mm211218: p2 => p2t
      p2t(1:nang,1:npol)=p2t(1:nang,1:npol)+phsfrz(1:nang,k,jatm(i,j),1:npol)*taurz(k,jatm(i,j))
      !p2(1:nang)=(p2(1:nang)+phsfr(1:nang)*taurz(k))/taustp(i,j,k)
      !end p2 => p2t
      !mm211218 truncation fraction
      fp(i,j,k) = (tausp(i,j,k)-taustp(i,j,k))/tausp(i,j,k)

! phase function table
      !mm211223 move to pflut
      call pflut2(ipol,0,nang,npol,ang(1:nang),p2t(1:nang,1:npol),na4,angr(1:na4),cs(1:na4),&
      p3t(1:na4,i,j,k,1:npol),angdp(1:na4,i,j,k),err)
      !mm211218 for tms
      if( imthd == 2 ) p3(1:nang,i,j,k,1:npol) = p2(1:nang,1:npol)
      if( imthd >= 3 .and. fp(i,j,k) .ne. 0. ) then
        !mm211220 for ims
! phase function table: total function
        call pflut2(ipol,1,nang,npol,ang(1:nang),p2(1:nang,1:npol),na4,angr(1:na4),cs(1:na4),&
        p32(1:na4,i,j,k,1:npol),angdps(1:na4,i,j,k),err)
! phase function table: forward peak
        p3s(1:na4,i,j,k,1:npol) = (p32(1:na4,i,j,k,1:npol)-(1-fp(i,j,k))*p3t(1:na4,i,j,k,1:npol))/fp(i,j,k)
        
! phase function table: total function
        call pflut2(ipol,0,na4,1,angr(1:na4),abs(p3s(1:na4,i,j,k,1)),na4,angr(1:na4),cs(1:na4),&
        p32(1:na4,i,j,k,1),angdps(1:na4,i,j,k),err)
        !call lgndf3(1,na4,cs,p3s(1:na4,i,j,k),g)
        !p3s(1:na4,i,j,k) = p3s(1:na4,i,j,k)/4/pi/abs(g(1))
        
        !end for ims
      end if

    enddo; enddo; enddo
! area averaged optical quantities
    taurv (1:nz)=taurv (1:nz)/(nx*ny)
    tauav (1:nz)=tauav (1:nz)/(nx*ny)
    tausav(1:nz)=tausav(1:nz)/(nx*ny)
    do k=1,nz
       g1=gav(1,k)
       gav(1:mxlgt1,k)=gav(1:mxlgt1,k)/g1
    enddo

! flxs1(1:na0,i,j) i=1: top, 2: boa, 3-6 sides; j=1 dn, j=2 up
    wg0=0

! CKDs
    if(icont(4)>0) then
      write(iuow,'(2i4,a)') nww,nch,' nww nch'
    endif

    nsos(1:ksos)=0
    ais1(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)=0
    aie1(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)=0
    flxs1(1:na0,1:6,1:3)=0
    flxe1(1:6,1:3)=0
    !mm220416 point flux
    flxsp1(1:na0,1:npv,1:3,1:3)=0
    flxep1(1:npv,1:3,1:3)=0

    do iww=1,nww
    do ich=1,nch

      wgt1=wgt(ich)
      wg0=wg0+wgt1

      do iatm = minval(jatm(1:nx,1:ny)), maxval(jatm(1:nx,1:ny))
        if( count(jatm(1:nx,1:ny)==iatm)==0 ) cycle
        tckdav(1:nz)=tckdav(1:nz)+tkdz(1:nz,ich,iww,iatm)*wgt1*count(jatm(1:nx,1:ny)==iatm)/(nx*ny)
      end do

! construction of 3D optical field
      do i=1,nx; do j=1,ny; do k=1,nz
        dz=abs(altz(k)-altz(k+1))
        cext(1:nx,1:ny,k)=(tautp(1:nx,1:ny,k)+tkdz(k,ich,iww,jatm(i,j)))/dz
        ssat(1:nx,1:ny,k)=taustp(1:nx,1:ny,k)/dz/cext(1:nx,1:ny,k)
        if( imthd == 2 ) then
          !mm231007 for user defined flux
          ssaf(1:nx,1:ny,k)=(tausp(i,j,k)-taustp(i,j,k))/dz
        else if( imthd >= 3 ) then
          !mm211220 for ims
          ssaf(1:nx,1:ny,k)=ssat(1:nx,1:ny,k)*fp(1:nx,1:ny,k)/(1-fp(1:nx,1:ny,k))
          !end for ims
        end if
      enddo; end do; end do
 
      write(*,*) "Radiance calculation"
! radiance calculation mode for vis and nir
      mlim(1) =klim; mlim(2) =klim
      mlim1(1)=klim; mlim1(2)=klim
      am0(1:na0)=cos(th0(1:na0)*rad)
      am1(1:na1)=cos(th1(1:na1)*rad)

! solar (lseg=1) and thermal (loseg=2) backward MC
      np=max(1000,np0/kback)
      if( imthd > 2 ) kims(1) = 2 !mm220213
      ais(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)=0  ! lpv-position; dn/up
      aie(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)=0
      flxsp(1:na0,1:npv,1:3,1:3)=0
      flxep(1:npv,1:3,1:3)=0  ! toa boa side, down-direct -diffuse up

      do lseg=1,2
        if(lseg==1) then
          if(isol<=0 .or. wl>6.0) cycle !mm220214
          na0se=na0
          nsp = 1
          if(imthd==3) nsp = 2 ! ims space
         elseif(lseg==2) then
          if(nplk1<=0) cycle
          na0se=1
          nsp = 1
        endif
! user-defined relative viewing positions (rpvw) at a boundary plane:
!   llu: 1 xy+z (toa), 2 xy-z (boa), 3 yz+x, 4 yz-x, 5 zx+y, 6 zx-y'

        do lpv=1,npv
          ll=llu(lpv)
          if(ll==0) ll=2
          if(icycl>0 .and. ll>2) cycle
          ll1=(ll+1)/2   ! boundary plane number
          ll2=ll-(ll1-1)*2  ! +boundary or -boundary
          i1=i33(ll1,1); i2=i33(ll1,2); i3=i33(ll1,3) ! coordinate axis numbers
          j1=nb(i1)+1; j2=nb(i2)+1  ! coordinates of +the boundary plane intercepts along i1, i2, i3 axis
          if(ll2==1) then
            j3=nb(i3)+1
           else
            j3=1
          endif
          if(ll==2 .and. llu(lpv)==0) ll=0

          do lj=1,na0se
! Rstar (th1) adops nadir angle mu1=cos(th1)>0 downward, <0 upward
! On the other hand, MC adopts zenith anlge   ds(3)<0 downward, >0 upward
!  mu1>0 upward, <0 downward
! solar direction
            sths=sin(th0(lj)*rad)
            drs(1)=sths*cos(fis*rad); drs(2)=sths*sin(fis*rad)
            drs(3)=am0(lj)
          
            write(*,*) "Loop in :", lj
          do lsp=1,nsp ! 1: delta-m, 2: ims spaces
          
          if(inda==1) then !mm220425
          
          do li=1,na1
            li2=mod(li-1,na12)+1

          do lk=1,nfi
! -viewing direction, for backward photon shooting from view point (receiver)
            sths=sin(th1(li)*rad)
            dvw(1)=sths*cos((fi(lk)+fis)*rad); dvw(2)=sths*sin((fi(lk)+fis)*rad)
            dvw(3)=am1(li) ! upward radiance in nadir polar coorinate (rstar) -> th1>90
                           !   -> dvw(3)<0 -> downward shooting in backward MC in zenith polar coord.
            if(j3>1  .and. dvw(i3)>0) cycle
            if(j3==1 .and. dvw(i3)<0 .and. i3/=3) cycle
            if(ll==0 .and. dvw(3)>0) cycle

            if(dvw(3)>0) then
              idui=1    ! down-going
             else
              idui=2    ! up-going
            endif
! photon loop
            if(icont(4)>0) then
              write(iuow,'(2i4,a)') lseg,np,' lseg np'
              write(iuow,'(6i4,a)') li,lj,lk,ll,ich,iww,' i j k ll ich iww'
            endif

            fb(1:ipol)=0
            fb2(1:ipol)=0.
            fbtrn(1:ipol)=0 !mm231005
            init2=1 ! initialize rdbx
            init3=1 ! initialize pvw set
            np1=0    ! valid photon number
            do lp=1,int(np/kims(lseg))
! backward solar ray tracing in the box system
              if(init3>0) then
                if(llr==0) then
                  pvw(i1)=bx(1,i1)+(bx(j1,i1)-bx(1,i1))*rpvw(lpv,1)
                  pvw(i2)=bx(1,i2)+(bx(j2,i2)-bx(1,i2))*rpvw(lpv,2)
                  init3=0
                else
                  nrnd=2
                  call rdmu2(dseed,nrnd,rndv)
                  pvw(i1)=bx(1,i1)+(bx(j1,i1)-bx(1,i1))*rndv(1)
                  pvw(i2)=bx(1,i2)+(bx(j2,i2)-bx(1,i2))*rndv(2)
                endif
                pvw(i3)=bx(j3,i3)
                if(ll==0) then
                  !write(*,'(A,1P3E12.4)') "before update: ", pvw(1:3)
                  call dsbnd1(icycl,-dvw,pvw,nb,bx,dbx,pb,cext,pbnd,mlim)
                  pvw(1:3) = pb(1:3)
                  !write(*,'(A,1P3E12.4)') "after update: ", pvw(1:3)
                end if
                call fbox(pvw,bx,nb,ibs,erc)
                galb2=ag(ibs(1),ibs(2))
                if(lsp==2) call dsbnd1(icycl,dvw,pvw,nb,bx,dbx,pb,cext,pbnd,mlim) ! for ims
              endif

              ind1=98
              if(lseg==1.and.lsp==1) then
                !mm211218 add tms
                call srdbxb2(init2,ipol,icycl,indg,dseed,drs,dvw,pvw,nb,bx,dbx,ag, &
                 cext,ssat,p3t,angdp,nns,nng,u1,u1t,utrn,lim1,mlim,ind1)
               else if(lseg==1.and.lsp==2) then ! ims
                 if( pbnd.ge.1. .and. init3.eq.0 ) exit
                 call imsbxb1(init2,ipol,icycl,dseed,pbnd,drs,dvw,pvw,nb,bx,dbx,&
                 cext,ssaf,p3s,p32(:,:,:,:,1),angdps,u1,utrn,lim1,mlim,ind1)
                 u1 = u1 * (1.-pbnd) ! retrieve
                 utrn = utrn * (1.-pbnd)
               else
                 call erdbxb2(init2,ipol,icycl,indg,dseed,dvw,pvw,nb,bx,dbx,nplk1,cplkzz,bgndz,ag, &
                   cext,ssat,p3t,angdp,nns,nng,u1,utrn,lim1,mlim,ind1)
              endif

! invalid photon
              if(iabs(ind1)<=3 .or. (itrn==0.and.ind1==7)) then
! photon output
                np1=np1+1
                !mm211218 for tms
                !mm220109: imthd==2
                !only short-wave
                if(lseg==1 .and. imthd==2) then
                  fb(1:ipol)=fb(1:ipol)+(u1(1:ipol)-u1t(1:ipol))
                  fb2(1:ipol)=fb2(1:ipol)+(u1(1:ipol)-u1t(1:ipol))**2.
                  if(utrn(1).ne.u1t(1)) fbtrn(1:ipol)=fbtrn(1:ipol)+utrn(1:ipol)
                else
                  fb(1:ipol)=fb(1:ipol)+u1(1:ipol)
                  fb2(1:ipol)=fb2(1:ipol)+u1(1:ipol)**2.
                  fbtrn(1:ipol)=fbtrn(1:ipol)+utrn(1:ipol)
                end if
              endif
            enddo  ! lp

            if(icont(4)>0) then
              write(iuow,'(i3,f10.5)') 99,0.0
            endif

            if(np1<=0.and.lsp==1) then
              write(*,'(a,4i4,2f6.1,3f8.3,a)') 'no photon',lpv,li,lj,lk,th1(li),fi(lk),dvw(1:3),' lpv i j k th1 fi dvw(1:3)'
              cycle
            endif

            if(lseg==1) then
              !mm211218 tms
              !substract single-scattering when imthd=2
              fb(1:ipol) = fb(1:ipol)*fsol/max(np1,1)
              fb2(1:ipol) = fb2(1:ipol)*(fsol**2.)/max(np1,1)
              ais(li2,lj,lk,lpv,idui,1:ipol,1)=ais(li2,lj,lk,lpv,idui,1:ipol,1) + fb(1:ipol)
              if(np1>1) ais(li2,lj,lk,lpv,idui,1:ipol,2)=ais(li2,lj,lk,lpv,idui,1:ipol,2) &
                + (fb2(1:ipol)-fb(1:ipol)**2.)/(max(np1,1)-1)
              ais(li2,lj,lk,lpv,idui,1:ipol,3)=ais(li2,lj,lk,lpv,idui,1:ipol,3) + fbtrn(1:ipol)*fsol/max(np1,1)
              !ais(li2,lj,lk,lpv,idui)=fb*fsol/np1
              !add single-scattering when imthd=2
              if( imthd == 2 ) then
                initt = 1
                !mm211220 add p
                !call tmsbxb1(initt,icycl,drs,dvw,pvw,nb,bx,cext,ssat,fp,nang,ang,p3,u1t,mlim,ind1)
                pow = -1.d0
                call tmsbxb2(initt,ipol,icycl,drs,dvw,pvw,nb,bx,dbx,cext,ssat,fp,&
                nang,knang,ang,p3,u1t,10,pow,mlim,ind1)
                ais(li2,lj,lk,lpv,idui,1:ipol,1)=ais(li2,lj,lk,lpv,idui,1:ipol,1)+u1t(1:ipol)*fsol
              end if
              !end tms
             else
              aie(li2,lk,lpv,idui,1:ipol,1)=fb(1:ipol)/np1
              if(np1>1) aie(li2,lk,lpv,idui,1:ipol,2)=(fb2(1:ipol)/np1-(fb(1:ipol)/np1)**2.)/(np1-1)
              aie(li2,lk,lpv,idui,1:ipol,3)=fbtrn(1:ipol)/np1
            endif

          enddo  ! lk
          enddo  ! li
          
          end if !inda
          
          ! point flux module
          do li=1,3 ! direct, up and downward
            fb(1:ipol)=0.
            fb2(1:ipol)=0.
            fbtrn(1:ipol)=0.
            init2=1
            init3=1
            np1=0    ! valid photon number
            if(lseg==2 .and. li==1) cycle
            do lp=1,int(np/kims(lseg))
! backward solar ray tracing in the box system
              if(init3>0) then
                if(llr==0) then
                  pvw(i1)=bx(1,i1)+(bx(j1,i1)-bx(1,i1))*rpvw(lpv,1)
                  pvw(i2)=bx(1,i2)+(bx(j2,i2)-bx(1,i2))*rpvw(lpv,2)
                  init3=0
                else
                  nrnd=2
                  call rdmu2(dseed,nrnd,rndv)
                  pvw(i1)=bx(1,i1)+(bx(j1,i1)-bx(1,i1))*rndv(1)
                  pvw(i2)=bx(1,i2)+(bx(j2,i2)-bx(1,i2))*rndv(2)
                endif
                pvw(i3)=bx(j3,i3)
                call fbox(pvw,bx,nb,ibs,erc)
                galb2=ag(ibs(1),ibs(2))
                if(indg(ibs(1),ibs(2))==0) then
                  galb2=ag(ibs(1),ibs(2))
                  bgnd2=bgndz(ibs(1),ibs(2))
                else
                  galb2=0.; bgnd2=0.
                end if
              endif
              
              call ggdir1(dseed,ds)
              if(li==2) then
                dvw(i1)=+ds(1); dvw(i2)=+ds(2); dvw(i3)=+ds(3)
                ds(i1)=0.; ds(i2)=0.; ds(i3)=1.
              else
                dvw(i1)=-ds(1); dvw(i2)=-ds(2); dvw(i3)=-ds(3)
                ds(i1)=0.; ds(i2)=0.; ds(i3)=-1.
              endif

              ind1=98
              if( lseg==1 .and. li==1 .and. lsp==1 ) then ! direct solar insolation
                call dsbnd1(icycl,drs,pvw,nb,bx,dbx,pb,cext,pbnd,mlim)
                if( imthd==2 ) then
                  call dsbnd1(icycl,drs,pvw,nb,bx,dbx,pb,ssaf,u1(1),mlim)
                  pbnd = u1(1)*pbnd
                end if
                if( init3 == 0 ) then
                  fb(1) = pbnd/2/pi*abs(drs(i3)); np1 = 1; exit
                else
                  u1 = 0.
                  u1(1) = pbnd/2/pi*abs(drs(i3)); dvw(i3) = 1.; ind1 =1
                end if
              else if( lseg==1 .and. lsp==1 ) then ! delta-M
                if( imthd==2 .and. (sum(drs(1:3)*ds(1:3))>0) ) then
                  call dsbnd1(icycl,drs,pvw,nb,bx,dbx,pb,cext,pbnd,mlim)
                  call dsbnd1(icycl,drs,pvw,nb,bx,dbx,pb,ssaf,u1(1),mlim)
                  if(abs(dvw(i3))>1.d-8) then
                    pbnd = (1-u1(1))*pbnd/2/pi*abs(drs(i3)/dvw(i3))
                  else
                    pbnd = 0.
                  end if
                end if
                call srdbxb2(init2,ipol,icycl,indg,dseed,drs,dvw,pvw,nb,bx,dbx,ag, &
                cext,ssat,p3t,angdp,nns,nng,u1,u1t,utrn,lim1,mlim,ind1)
                if(imthd==2 .and. (sum(drs(1:3)*ds(1:3))>0)) u1 = u1 + pbnd
              else if( lseg==1 .and. lsp==2 ) then ! IMS
                call dsbnd1(icycl,dvw,pvw,nb,bx,dbx,pb,cext,pbnd,mlim)
                if( pbnd .ge. 1. ) then
                  u1 = 0.; ind1 = 1; utrn =0.
                else
                  call imsbxb1(init2,ipol,icycl,dseed,pbnd,drs,dvw,pvw,nb,bx,dbx,&
                  cext,ssaf,p3s,p32(:,:,:,:,1),angdps,u1,utrn,lim1,mlim,ind1)
                  u1 = u1 * (1.-pbnd)
                  utrn = utrn * (1.-pbnd)
                end if
              else
                call erdbxb2(init2,ipol,icycl,indg,dseed,dvw,pvw,nb,bx,dbx,nplk1,cplkzz,bgndz,ag, &
                cext,ssat,p3t,angdp,nns,nng,u1,utrn,lim1,mlim,ind1)
                if( ll==2 .and. li==3 ) u1(1) = u1(1) + (1-galb2)*bgnd2 *0.5
              end if

! invalid photon
              if(iabs(ind1)<=3 .or. (itrn==0.and.ind1==7)) then
! photon output
                np1=np1+1
                fb(1:ipol)=fb(1:ipol)+u1(1:ipol)*abs(dvw(i3))
                fb2(1:ipol)=fb2(1:ipol)+(u1(1:ipol)*abs(dvw(i3)))**2.
                fbtrn(1:ipol)=fbtrn(1:ipol)+utrn(1:ipol)*abs(dvw(i3))
              endif

            end do  ! lp
            if(lseg==1) then
              fb(1:ipol) = fb(1:ipol)*fsol*2*pi/max(np1,1)
              fb2(1:ipol) = fb2(1:ipol)*((fsol*2*pi)**2.)/max(np1,1)
              flxsp(lj,lpv,li,1) = flxsp(lj,lpv,li,1) + fb(1)
              if(np1>1) flxsp(lj,lpv,li,2) = flxsp(lj,lpv,li,2) + (fb2(1)-fb(1)**2.)/(max(np1,1)-1)
              flxsp(lj,lpv,li,3) = flxsp(lj,lpv,li,3) + fbtrn(1)/max(np1,1)*fsol*2*pi
              ! delta-M space to real space via IMS space : energy conservation
              if( lsp==2 .and. li > 1 ) then
                if(dvw(i3)<0) then ! upward
                  flxsp(lj,lpv,1,1) = flxsp(lj,lpv,1,1) +fb(1)
                  flxsp(lj,lpv,1,3) = flxsp(lj,lpv,1,3) + fbtrn(1)/max(np1,1)*fsol*2*pi
                else ! downward
                  flxsp(lj,lpv,1,1) = flxsp(lj,lpv,1,1) - fb(1)
                  flxsp(lj,lpv,1,3) = flxsp(lj,lpv,1,3) - fbtrn(1)/max(np1,1)*fsol*2*pi
                end if
                if(np1>1) flxsp(lj,lpv,1,2) = flxsp(lj,lpv,1,2) + (fb2(1)-fb(1)**2.)/(max(np1,1)-1)
              end if
            else
              fb(1:ipol) = fb(1:ipol)*2*pi/max(np1,1)
              fb2(1:ipol) = fb2(1:ipol)*((2*pi)**2.)/max(np1,1)
              flxep(lpv,li,1) = flxep(lpv,li,1) + fb(1)
              if(np1>1) flxep(lpv,li,2) = flxep(lpv,li,2) + (fb2(1)-fb(1)**2.)/(max(np1,1)-1)
              flxep(lpv,li,3) = flxep(lpv,li,3) + fbtrn(1)/max(np1,1)*2*pi
            end if
          end do ! li
          end do ! lsp
          enddo  ! lj
        enddo  ! lpv
      enddo  ! lseg

      if(inda==1) then
      ais1(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)=ais1(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)   &
                                   +ais(1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)*wgt1

      aie1(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)=aie1(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)   &
                                   +aie(1:na12,1:nfi,1:npv,1:2,1:ipol,1:3)*wgt1
      end if !inda

      flxsp1(1:na0,1:npv,1:3,1:3)=flxsp1(1:na0,1:npv,1:3,1:3)+flxsp(1:na0,1:npv,1:3,1:3)*wgt1

      flxep1(1:npv,1:3,1:3)=flxep1(1:npv,1:3,1:3)+flxep(1:npv,1:3,1:3)*wgt1

      write(*,*) "Flux calculation with Forward MC"
! Forward solar MC (lseg=3)
      lseg=3
      flxs0(1:na0,1:6)=0   ! solar insolation flux for each illuminating surface lli
      fbf1(1:na0,1:6,1:3)=0
      flxuxy1(1:na0,1:nx,1:ny)=0

      if(lseg==3 .and. isol>0 .and. wl<6.0) then

! flux calculation mode for vis and nir
        do  lj=1,na0
! Polar coordinate with positive z-direction toward zenith
!  mu1>0 upward, <0 downward
! solar direction along x-axis, fis>0 for
          sths=sin(th0(lj)*rad)
          drs(1)=sths*cos(fis*rad); drs(2)=sths*sin(fis*rad)
          drs(3)=-am0(lj)
          earea0=abs((bx(nb(1)+1,1)-bx(1,1))*(bx(nb(2)+1,2)-bx(1,2)))  ! target area

          do lli=1,6 ! toa, boa, four lateral sides for solar illumination
            if(lli==2) cycle
            if(icycl>0 .and. lli>2) cycle
            if(lli==1) then
              np=max(1000,np0/2)
             else
              np=max(1000,np0/4)
            endif
            ll1=(lli+1)/2   ! boundary plane number
            ll2=lli-(ll1-1)*2  ! +boundary or -boundary
            i1=i33(ll1,1); i2=i33(ll1,2); i3=i33(ll1,3) ! coordinate axis numbers
            j1=nb(i1)+1; j2=nb(i2)+1  ! coordinates of +the boundary plane intercepts along i1, i2, i3 axis
            if(ll2==1) then
              j3=nb(i3)+1
             else
              j3=1
            endif

            if(j3>1  .and. drs(i3)>0) cycle
            if(j3==1 .and. drs(i3)<0 .and. i3/=3) cycle

            earea=abs(drs(i3)*(bx(j1,i1)-bx(1,i1))*(bx(j2,i2)-bx(1,i2)))

            if(icont(4)>0) then
              write(iuow,'(7i4,a)') lseg,np,' lseg np'
              write(iuow,'(6i4,a)') li,lj,lk,lli,ich,iww,' i j k lli ich iww'
            endif

! photon loop
            fbf(1:6,1:3)=0
            flxuxyb(1:nx,1:ny)=0
 
            np1=0
            init2=1
            do lp=1,np
              nrnd=2
              call rdmu2(dseed,nrnd,rndv)
              pvw(i1)=bx(1,i1)+(bx(j1,i1)-bx(1,i1))*rndv(1)
              pvw(i2)=bx(1,i2)+(bx(j2,i2)-bx(1,i2))*rndv(2)
              !pvw(3)=bx(nb(3)+1,3) !mm211218
              pvw(i3)=bx(j3,i3)

! ray tracing in the box system
              ind1=98
              call srdbxf3(init2,icycl,indg,dseed,drs,pvw,nb,bx,dbx,ag, &
                cext,ssat,angdp,nns,nng,u,ds,ps,ii,jj,fdg,fdg0,fug, &
                nnsg,ug,dsg,psg,ibxg,ibyg,lim1,mlim,ind1)
! photon output
              if(icont(4)>0) then
                write(iuow,'(3i4,7f11.5)') ind1,nns,nng,u,fdg,fdg0,fug,ds(1:3)
                if(nng>0) then
                  do i=1,nng
                    write(iuow,'(4i6,11f11.5)') i,ibxg(i),ibyg(i),nnsg(i),ug(i,1:2),dsg(i,1:3,1:2),psg(i,1:3)
                  enddo
                endif
              endif

! invalid photon
!!          if(iabs(ind1)>3) cycle ! lp

! extincted photon with u=0
              if(ind1==0 .or. (itrn==0.and.ind1==7)) then
                np1= np1+1
!tr211213 (no-need)           fbf(1,1)=fbf(1,1)+1
                cycle
              endif

! fluxes and radiance
! flux  fbf(i,x): i=1 toa, 2 boa, 3 lateral sides
! toa and lateral sides
              ierr=0
              i=0; j=0; k=iabs(ind1)
              if(ind1==3) then   ! toa
                if(ds(3)<=0) ierr=1
                i=1; j=3
               else if(k<=2) then   ! lateral sides
                i=2*(k+1)
                if(ind1>0) i=i-1
                if(ds(3)>0) then
                  if(nns<=0) then
                    j=1  ! downward direct from lateral sides
                   else
                    j=2  ! downward diffuse from sides
                  endif
                 else
                  j=3  ! upward from lateral sides
                endif
              endif
              if(ierr>0) cycle
! flux sum
              if(i>0 .and. j>0) then
                fbf(i,j)=fbf(i,j)+u
                if(i==1) flxuxyb(ii,jj)=flxuxyb(ii,jj)+u
              endif

! boa flux and radiances at ground
              i=2
              fbf(i,1)=fbf(i,1)+fdg0
              fbf(i,2)=fbf(i,2)+fdg
              fbf(i,3)=fbf(i,3)+fug

! valid scattered photon counting
              np1=np1+1
              fbf(1,1)=fbf(1,1)+1

! scattering order statistics
              k=nns+1
              if(k<=ksos) nsos(k)=nsos(k)+1

! np-loop end
            enddo  ! lp=1,np
            if(icont(4)>0) then
              write(iuow,'(3i4,7f11.5)') 99,0,0,0.0,0.0,0.0,0.0,0.0
            endif

            if(np1<=0) then
              write(*,*) 'no photon'
              cycle
            endif

            c=fsol*earea/earea0/np1
            fbf1(lj,1:6,1:3)=fbf1(lj,1:6,1:3)+fbf(1:6,1:3)*c
            flxuxy(1:nx,1:ny)=flxuxy(1:nx,1:ny)+flxuxyb(1:nx,1:ny)*c
            flxs0(lj,lli)=flxs0(lj,lli)+fsol*earea/earea0*wgt1

          enddo ! lli
! fluxes
!! fb: calculated by the present photon counting

! ckd integration
! flxs1(lk,ll,lud)  ll=toa, boa, side,  lud=direct down,diffuse down, up
            flxs1(lj,1:6,1:3)=flxs1(lj,1:6,1:3)+fbf1(lj,1:6,1:3)*wgt1
            flxuxy1(lj,1:nx,1:ny)=flxuxy1(lj,1:nx,1:ny)+flxuxy(1:nx,1:ny)*wgt1
        enddo  ! lj

      endif ! lseg=3

! Thermal fluxes by randome sampling (lseg=4)
      lseg=4
      np=max(1000,np0/kback)
      flxe(1:6,1:3)=0  ! toa boa side, down-direct -diffuse up
      flxe(2,3)=(1-galb3)*bgnd3*pi

      if(lseg==4 .and. nplk1>0) then
! photon loop
        earea0=abs((bx(nb(1)+1,1)-bx(1,1))*(bx(nb(2)+1,2)-bx(1,2)))  ! target area

        do ll=1,6 ! toa, boa, four lateral sides
          if(icycl>0 .and. ll>2) cycle
          ll1=(ll+1)/2
          ll2=ll-(ll1-1)*2
          i1=i33(ll1,1); i2=i33(ll1,2); i3=i33(ll1,3)
          j1=nb(i1)+1; j2=nb(i2)+1
          if(ll2==1) then
            j3=nb(i3)+1
           else
            j3=1
          endif
 
          nrnd=2
 
          fbe(1:3)=0
          earea=abs(drs(i3)*(bx(j1,i1)-bx(1,i1))*(bx(j2,i2)-bx(1,i2)))

          if(icont(4)>0) then
            write(iuow,'(7i4,a)') lseg,np,' lseg np'
            write(iuow,'(6i4,a)') li,lj,lk,ll,ich,iww,' i j k ll ich iww'
          endif

          init2=1
          np1=0    ! valid photon number
          do lp=1,np
! backward solar ray tracing in the box system
            call rdmu2(dseed,nrnd,rndv)
            pvw(i1)=bx(1,i1)+(bx(j1,i1)-bx(1,i1))*rndv(1)
            pvw(i2)=bx(1,i2)+(bx(j2,i2)-bx(1,i2))*rndv(2)
            pvw(i3)=bx(j3,i3)

            call ggdir1(dseed,ds)
            if(j3==1) then
              dvw(i1)=+ds(1); dvw(i2)=+ds(2); dvw(i3)=+ds(3)
             else
              dvw(i1)=-ds(1); dvw(i2)=-ds(2); dvw(i3)=-ds(3)
            endif

            ind1=98
            call erdbxb2(init2,1,icycl,indg,dseed,dvw,pvw,nb,bx,dbx,nplk1,cplkzz,bgndz,ag, &
                 cext,ssat,p3t,angdp,nns,nng,u1,utrn,lim1,mlim,ind1)

            if(icont(4)>0) then
              write(iuow,'(i3,f10.5)') ind1,u1
            endif

! invalid photon
            if(iabs(ind1)<=3 .or. (itrn==0.and.ind1==7)) then
! photon output
              np1=np1+1
              if(dvw(3)<0) then
                fbe(3)=fbe(3)+u1(1)*abs(dvw(i3))
               else
                fbe(2)=fbe(2)+u1(1)*abs(dvw(i3))
              endif
            endif

          enddo  ! lp

          if(icont(4)>0) then
            write(iuow,'(i3,f10.5)') 99,0.0
          endif

          if(np1<=0) then
            write(*,*) 'no photon'
            cycle
          endif

          flxe(ll,1:3)=flxe(ll,1:3)+fbe(1:3)/np1*2*pi*earea/earea0
        enddo ! ll

        flxe1(1:6,1:3)=flxe1(1:6,1:3)+flxe(1:6,1:3)*wgt1

      endif  ! lseg=4

! CKD-loop end
    enddo ! ich
    enddo  ! iww

    tckdav(1:nz)=tckdav(1:nz)/wg0
    if(inda==1) then
    ais1 (1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)=ais1 (1:na12,1:na0,1:nfi,1:npv,1:2,1:ipol,1:3)/wg0
    aie1(1:na12,      1:nfi,1:npv,1:2,1:ipol,1:3)=aie1(1:na12,      1:nfi,1:npv,1:2,1:ipol,1:3)/wg0
    end if
    flxsp1(1:na0,1:npv,1:3,1:3)=flxsp1(1:na0,1:npv,1:3,1:3)/wg0 !mm220416
    flxep1(1:npv,1:3,1:3)=flxep1(1:npv,1:3,1:3)/wg0 !mm220416
    flxe1(1:6,1:3)=flxe1(1:6,1:3)/wg0
    flxs1(1:na0,1:6,1:3)=flxs1(1:na0,1:6,1:3)/wg0
    flxs0(1:na0,1:6)=flxs0(1:na0,1:6)/wg0
    flxuxy1(1:na0,1:nx,1:ny)=flxuxy1(1:na0,1:nx,1:ny)/wg0

  return
end subroutine mcst2