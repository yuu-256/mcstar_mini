!mm220221
subroutine rstr7cms(init,icont,iuk,iup,iukd,iug,indg,inda,isol,ipol,na0, &
     th0,na1,th1,nfi,wl,icn,wlcn,npoly, &
     ncomp,mptc,vptc,cnpt,nl,alt,prs,tmp,cng,ins,cnp, &
     ispcvp,rfracp,asphr,rop,dryap,nawp,awcrp,nwlv,wlv,rfi, &
     nln,ipbf,ngas,thk0,fsol,flxd0,flxd00,flxu0,ai0,nch,nww, &
     wgt,ang,ssaa,tkd,taur,taua,ssapoly,epoly, &
     ppoly,phsm,dpf,nang,err)

! general subroutine for calculating optical properties at
! wavelength wl and CDK, dropped nw and rx and rf; nww=1
!--- history
! 95. 9.15 created
! 95.12.28 corrected the comments for icn=1, and added icn=3 case.
! 95.12.30 modified the relative humidity profile calculation
! 96. 3.11 added non-spherical parameters
!     5. 5 introduce ncomp and some bug-fix for cnprf1
!     5. 9 cextp() -> cextp2, cabsp() -> cabsp2, voldp()->voldp2, volp()->volp2
!     5.24 change the structure of do 18
! 96.12.27 bug in initialization with do 17
! 96.12.28 change the meaning of rx
! 97. 1.14  20.0->200.0
! 97. 3.17  get scr and sci from subroutine rfidxb.
!           call rtrn21 with scr and sci.
! 98.1.7
!   u10 > 0.01 for rstr and pstr.  u10 < 0.01 will have significant
!   error in radiances due to angular integration errors for surface
!   reflection.  recommend that u10 is set to be 0.01 for flat surface
!   case
! 06.06.25 version up to 6. (miho)
! 07.03.28 fixed some bugs. (miho)
! 07.11.30 dynamical => statical dimension (nln=>knln, nch=>kch) (miho)
! 08. 2.19 flxd00 for direct solar flux
! 08. 2.19 added notes for thk0
! 08. 2.19 deleted ifrh and trh
! 08. 8.14 rtrn2 --> rtrn21, add dintvl; ssaa1 -> ssaa in subroutine intefrace
!tr210905 rstr6bb, dropped tau, ssapoly1, ssapoly, ppolyt-> ppoly,ttkd, ttaur
!tr  rd in use parameter
!tr211026 itray, itkd
!tr211123 rename as rstr6cms, icont, gtmp
!--- input
! init       i        if 1 initilize the routine
!                     0 skip wavelength-independent part
! icont   i(10)       1: zero for no-Rayleigh (itray); 2: zero for no-ckd (itkd)
!                     3: zero for L(0) radiance (immax)
! iuk        i        device number for a kernel file
! indg       i       -1: no ground surface
!                     0: lambert surface
!                     1: ocean surface initialization
!                     2: ocean surface with no initialization
!                     when indg>0 and imthd>0 then cabspsingle scattering correc
!                      for ocean surface reflection
! inda       i        0: flux; 1: radiance
! imthd      i       -1: nt,  0: dms-method  for intensity/flux
!                     1: ms,  2:tms,  3:ims-method for intensity.
!                     when indg>0 and imthd>0 then single scattering correction
!                      for ocean surface reflection
! isol       i        0: night
!                     1: day
! nda        i        nunber of nadir-quadrature angles in the
!                     hemisphere.
! na0        i        number of solar incidences.
! th0     r(na0)      solar zenith angles 0-90 degrees
! na1       i        number of emergent nadir angles in sphere.
! th1      r(kna1u)    emergent nadir angles.
!                      0-90 for downward; 90-180 for upward
! nfi        i        number of azimuthal angles.
! fi     r(knfi)      azimuthal angles in degrees.
! wl         r        center wavelengths (micron).
! dw         r        scaling factor to define scaled wavelengths rx.
!                     give 0 for monochromatic calculations.
! galb       r        ground albedo if indg=0
!                     u10 (m/sec)   if indg>0; u10>0.01
! icn        i        0: cnpt is the total column dry volume for each polydisper
!                     1: cnpt is the total column volume for each polydispersion
!                     2: cnpt is optical thickness at wlcn
!                     3: cnpt is optical thickness at each wavelength.
! wlcn       r        scaling wavelength for icn=2,3 (cm)
! npoly      i        number of external mixtures
! ncomp   i(npoly)    number of aerosol types for the external mixture
! mptc    i(kpoly,4)  aerosol types for the external mixture
! vptc    r(kpoly,4)  dry volume ratio for the mixture
! cnpt    r(npoly)    total volume (icn=0,1) or optical thickness (icn=2,3)
!                      for each polydispersion.
! nl         i        number of atmospheric layers to define model atmosphere
! alt     r(knl)      height at interfaces of layers (km), top to bottom
! prs     r(knl)      pressure    at interfaces of layers (mb)
! tmp     r(knl)      temperature at interfaces of layers (k)
! gtmp       r        ground temperature
! nmol       i        number of gases
! cng    r(knl,knm0)  gas concentration (ppmv)
! nptc       i        number of defined aerosol types
! cnp    r(knl,kptc)  dry volume concentration profile (relative unit)
! ispcvp i(3,kptc)    fundamental materials for 3 comp. internal mixture (1-8)
! rfracp r(3,kptc)    dry volume fraction of the dry mixture
! asphr  r(3,kptc)    non-spherical parameters (x0, g, r)
! rop    r(kptc)      paricle density relative to water
! dryap  r(6,4,kptc) dv/ƒaadlnr parameters for the dry mixture
!                     see vlspc2
!                     c-values (coefficients of volume spectrum) are ralative
! nawp   i(kptc)      number of aw
! awcrp r(kaw,kptc,2)1: aw      water activity (see shettle and fenn.)
!                     3: rmmd    rmmd
! nv         i        number of fundamental species (1-8)
! nwlv       i        number of wavelengths for refractive index tables
! wlv    r(kwl)       wavelengths (micron)
! rfi   r(kwl,knv,2)  refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                     with log-regular wavelength interval
! nln        i        number of atmospheric layers for radiative transfer.
! ipbf    i(nln1)     interface number to define sublayers to construct
!                      the transfer atmosphere in mlatm.
!                     ipbf(1)=50, ipbf(nln+1)=1
!                     if >0 then field calculations
!                        <0 then no field calculations
!--- output
! thk0    r(nln,10)   (l,1): total optical thickness for extinction (tau)
!                     (l,2): particle optical thickness for extinction (taua)
!                     (l,3): optical thickness for rayleigh scattering (taur)
!                     (l,4): single scattering albedo (w)
!                     (l,5): g1 (1st legendre polynomial moment)
!                     (l,6): g2 (2nd legendre polynomial moment)
!   these values are not for monochromatic but for wavelength-mean values
!   other important parameters can be obtained as follows:
!   optical thickness for gas absorption: taug= tau - taua - taur
!   single scattering albedo of particles: wa = (w*tau-taur)/taua
!   legendre poly. moments of particles: gan= (w*tau*gn - grn*taur)/(wa*taua)
!    because  tau= taua + taur + taug
!             w*tau = wa*taua + taur
!             w*tau*gn = wa*taua*gan + grn*taur
!    where gr is the phase function of rayleigh scattering phase function
!   good approximation is: gr1= 0 and gr2= 0.1
! sol        r        solar incident irradiance (w/m2/micron)
! flxd0 r(kna0,kntau) downward flux (w/m2/micron)
! flxd00 r(kna0,kntau) downward direct flux (horizontal plane flux)
! flxu0 r(kna0,kntau) upward   flux (w/m2/micron)
! ai0   r(kna1u,kna0,knfi,kntau)
!                     radiance (w/m2/micron/str)
! err       c*64      error code.  if ' ' then normal termination.
! add output data 120513
! ssaa  single scattering albedo parameter
! ttkd  gas absorption mass
!---
!120520
use paras, only: kna0,kna1u,knfi,kpoly,knl,knm0,kptc, &
     kaw,kwlv,knv,knln1,knln,kntau,kmol1,kch,kww,knang,kplk1,kmol, &
     klgn1,pstd,pi,rad,kpol2, istar
implicit none

! control parameter of the programs
! 1: zero for no-Rayleigh (itray),2: zero for no-ckd (itkd),3: zero for L(0) radiance (immax)
  integer:: icont(10)

  integer, intent(inout) :: init, inda
  integer, intent(in) :: iuk,iug,indg
  integer, intent(in) :: iup,iukd !mm220221
  integer, intent(in) :: isol,ipol,na0,na1,nfi
  real(8), intent(in) :: th0(kna0),th1(kna1u)
  real(8), intent(in) :: wl!,dw
!  real(8), intent(in) :: galb
  integer, intent(in) :: icn
  real(8), intent(in) :: wlcn
  integer, intent(in) :: npoly
  integer, intent(in) :: ncomp(kpoly),mptc(kpoly,4)
  real(8), intent(in) :: vptc(kpoly,4),cnpt(kpoly)
  real(8), intent(in) :: alt(knl),prs(knl),tmp(knl)
  real(8), intent(in) :: cng(knl,knm0)
  integer, intent(in) :: ins(kptc) !mm2202221
  real(8), intent(in) :: cnp(knl,kptc)
  integer, intent(in) :: ispcvp(3,kptc)
  real(8), intent(in) :: rfracp(3,kptc),asphr(3,kptc),rop(kptc), &
                      dryap(6,4,kptc)
  integer, intent(in) :: nawp(kptc)
  real(8), intent(in) :: awcrp(kaw,kptc,2)
  integer, intent(in) :: nwlv
  real(8), intent(in) :: wlv(kwlv),rfi(kwlv,knv,2)
  integer, intent(in) :: nln,ipbf(knln1)
  integer  :: nch
  real(8), intent(inout)  ::thk0(knln,10)
  real(8) :: sol
  real(8), intent(inout) :: flxd0(kna0,kntau),flxd00(kna0,kntau), &
                         flxu0(kna0,kntau), ai0(kna1u,kna0,knfi,kntau)
  integer, intent(inout) :: nww
  real(8), intent(inout) :: wgt(kch),ang(knang)
!  real(8), intent(inout) :: prb(knang,kpoly),prbr(knang)
  real(8), save  :: scaa(knln)
  real(8)  :: tkd(knln,kch,kww),taur(knln)
  integer, intent(inout) :: nang
  character(len=64), intent(inout) :: err
!  save
!work area
  integer, save :: indp,nln1,ntau,l,indt
  real(8), save :: epsp,epsu,vptc2,rhl,pp,tt,dz,dzcm0,p
  integer, save :: j,k,ll,l1,ipoly,i,ik,iww
  real(8), save :: t,ppmv,rh1,gm3,e,gm3s,scr,sci
  integer, save :: m,ind,initrf
  real(8), save :: voldp2,volp2,cextp2,cabsp2
  real(8), save :: wg0,es,p2,t2,rh2,wgt1,wlc,wn
  real(8), intent(out) :: fsol
  integer, save :: nplk1
  real(8), save :: tray,cext1,cabs1,csca1,fac,sum
  integer, save :: nang0
  real(8), save :: ang0(knang)
!tr
real(8):: ssaa(knln),spoly(knln,kpoly),epoly(knln,kpoly),  &
        ssapoly(knln,kpoly),ppoly(knang,knln,kpoly,kpol2)!,bplnk,gtmp

! local
! for scald
! for abgask
! 07.11.30 nln => knln, nch => kch (miho)
  real(8), save :: amtpb(kmol1,knln)
! for qgausn
! for sdp2
! 07.03.28 nln => knln (miho, pointed out by fukuda)
  real(8), save :: pl(knln),tl(knln),dp(knln)
  integer, save :: np(knln)
! for rtrn21
  real(8),save :: thk(knln),omg(knln), &
       am0(kna0),am1u(kna1u)
  integer, save :: ipb(knln1)
!
  real(8), save :: rh(knln),cnml(kpoly,4),cnprf1(knln,kpoly,4), &
       vptc1(kpoly,4),php(knang,kpol2) !mm220221 php <= kpol2
  real(8)::taua(knln)
  real(8) :: phsm(knang,knln,kpol2)
  real(8), save :: volt(kpoly),cextpn(kpoly)
  real(8), save :: taua3(knln),scaa3(knln)
  real(8) :: p1,t1
!  real(8) :: tau(kx,ky,kz)
  real(8) :: taumax,dpt,b2,b1,wgt2,thks,bgnd
  integer :: lf
  integer,intent(in):: ngas(2)
  integer,save:: initg=1 !mm220617 from Dai
!  real(8) :: cot(kz), cf(kz)
!tr
  integer:: nl,iw99
  real(8):: sunir,fx
!mm220221
  integer :: nv
  integer,save :: npol
!mm220222
! for qgausn
  real(8):: gwt(kch)             !! Gaussian weights 
  real(8):: gmu(kch)             !! Gaussian points
! for abgask
  real(8) :: dintvl
!mm220222 for rayleigh on rpstar
  real(8) :: dzcm(knln), co2(knln), dpf(knln),g
!---
  if(init.eq.1) then
     init=0
     !initg=1 !mm220617 moved
! 96.12.27       ntau   = 2
     if(inda.ge.1) inda   = 1
     indp   = 1
     epsp   = 0
     epsu   = 0
!mm220221     ipol   = 1
     nln1=nln+1
! level checking
     ntau=0
! raylegh scattering coefficient
     do  l=1,nln1
        ipb(l)=abs(ipbf(l))
        if(ipbf(l).gt.0) ntau=ntau+1
     enddo
! 96.12.27 added t. nakajima
     if(ntau.le.0) then
        err='no interface for calculating radiation field'
        return
     endif
! 96.12.27 end
     if(ntau.gt.kntau) then
        err='number of field calculation layers larger than kntau'
        return
     endif
     if(ntau.eq.nln1) then
        indt=2
     else
        indt=1
     endif

! level merging (from top to bottom)
!     gtmp=tmp(1)  ! gtmp is given by user
!!
     do  ipoly=1,npoly
        do  j=1,ncomp(ipoly)
           vptc2=vptc2+vptc(ipoly,j)
           cnml(ipoly,j)=0
           do ll=1,nln
              cnprf1(ll,ipoly,j)=0
           enddo
        enddo
        do  j=1,ncomp(ipoly)
           vptc1(ipoly,j)=vptc(ipoly,j)/vptc2
        enddo
     enddo

     amtpb(1:kmol1,1:nln) = 0.

     do  ll = 1 ,nln
        rhl=0.d0
        pp=0.d0
        tt=0.d0
        co2(ll)=0.d0 !mm220222
        do ipoly = 1 ,npoly
           do j=1,ncomp(ipoly)
              cnprf1(ll,ipoly,j)=0
           enddo
        enddo

!        do 5  k = 1 ,53
!    5   amtpba(k,ll) = 0
!
! 07.03.28 bug fixed. (miho, pointed out by dr. kurino)
!! merging

        do l=ipb(ll+1),ipb(ll)-1
           l1=l+1
           dz=abs(alt(l)-alt(l1))
           dzcm0=dz*1.0e5
           if( istar(1) == 0 ) then ! Rstar
             p =(prs(l)+prs(l1))/2
           else if( istar(1) == 1 ) then ! Rpstar
             p =sqrt(prs(l)*prs(l1)) !mm220222 from rpstar
           end if
           t =(tmp(l)+tmp(l1))/2
! 06.03.17 miho
           if( istar(1) == 0 ) then ! Rstar
             pp=pp+p
           else if( istar(1) == 1 ) then ! Rstar
             pp=pp+p*dz !mm220222 from rpstar
           end if
           tt=tt+t
           dz=8.31451d0*273.15d0*(prs(l)-prs(l1))/(1013.25d0*9.8*28.8)

! 06.03.17 end
! gas amount
           do m=1,kmol
!    9   wm(m)=(cng(l,m)+cng(l1,m))/2
! 1.d-1=1.d-6*1.d5
              amtpb(m,ll)=amtpb(m,ll) &
                   +0.5d0*(cng(l,m)+cng(l1,m))*dz*1.d-1
           enddo

           amtpb(kmol1,ll)=amtpb(1,ll)**2/(amtpb(1,ll)+dz*1.d5)

! relative humidity (95.12.30)
           ind=2
           ppmv=cng(l,1)
           p1=prs(l)
           t1=tmp(l)
           call wvcal(ind,p1,t1,rh1,ppmv,gm3,e,gm3s,es)
           ppmv=cng(l1,1)
           p2=prs(l1)
           t2=tmp(l1)
           call wvcal(ind,p2,t2,rh2,ppmv,gm3,e,gm3s,es)
           rhl=rhl+(rh1+rh2)/2

!! scaled gas amount
!        call scald(p,t,wm,ibmol,densty)
!        do 8 k  = 1 ,53
!    8   amtpba(k,ll) = amtpba(k,ll) + densty(k)*dz
!! total dry concentration
           do ipoly=1,npoly
              do j=1,ncomp(ipoly)
                 m=mptc(ipoly,j)
                 cnprf1(ll,ipoly,j)=cnprf1(ll,ipoly,j) &
                      +(cnp(l,m)+cnp(l1,m))/2*dzcm0
              enddo
           enddo
           
           co2(ll)=co2(ll)+(cng(l,2)+cng(l1,2))*0.5d0*abs(alt(l)-alt(l1)) !mm220222
        enddo

        do ipoly=1,npoly
           do j=1,ncomp(ipoly)
              cnml(ipoly,j)=cnml(ipoly,j)+cnprf1(ll,ipoly,j)
           enddo
        enddo
        
        dzcm(ll)=alt(ipb(ll))-alt(ipb(ll+1))
! 06.03.17 miho
        if( istar(1) == 0 ) then ! Rstar
          pl(ll)=pp/abs(ipb(ll)-ipb(ll+1))
        else if( istar(1) == 1 ) then ! Rpstar
          pl(ll)=pp/dzcm(ll)
        end if
!end
        tl(ll)=tt/abs(ipb(ll)-ipb(ll+1))
! 06.03.17 end
        rh(ll)=rhl/abs(ipb(ll)-ipb(ll+1))
!mm220222 for king factor
! King factor
        co2(ll)=co2(ll)/abs(alt(ipb(ll))-alt(ipb(ll+1)))
        dzcm(ll)=dzcm(ll)*1.d5
!end king factor
     enddo
     
! extinction at a scaling wavelength, wlcn
     do ipoly=1,npoly
        volt(ipoly)=0
        cextpn(ipoly)=0
        do j=1,ncomp(ipoly)
           m=mptc(ipoly,j)
! dry volume profiles
           initrf=1
           do l=1,nln

              rhl=rh(l)
              if(cnprf1(l,ipoly,j).gt.0 .and. initrf.gt.0) then
                 !call rfidxb(iuk,wlcn,rhl,m,ispcvp,rfracp,asphr, &
                 !     rop,dryap,nawp,awcrp,nwlv,wlv,rfi,scr,sci, &
                 !     voldp2,volp2,cextp2,cabsp2,nang,ang,php,err)
                 call rfidx7p(iuk,iup,iukd,ipol,ins(m),wlcn,rhl,m,&
                      ispcvp,rfracp,asphr,rop,dryap,nawp,awcrp, &
                      nv,nwlv,wlv,rfi,scr,sci, &
                      voldp2,volp2,cextp2,cabsp2,nang,ang,php,err)
                 if(err.ne.' ') return
                 if(nawp(m).le.0) initrf=0
              endif
              if(voldp2*cnml(ipoly,j).gt.0)then
! 96.5.5 cnprf1(l,ipoly)=cnprf1(l,ipoly)/cnml(ipoly)/voldp2
                 cnprf1(l,ipoly,j)=cnprf1(l,ipoly,j)/cnml(ipoly,j) &
                      *vptc1(ipoly,j)
                 cextpn(ipoly)=cextpn(ipoly)+cextp2/voldp2 &
                      *cnprf1(l,ipoly,j)
                 volt(ipoly)=volt(ipoly)+volp2/voldp2*cnprf1(l,ipoly,j)
              else
                 cnprf1(l,ipoly,j)=0
              endif
           enddo
        enddo
     enddo
!mm220416 move to rayleigh scattering loop
! rayleigh phase function
!     do i = 1 ,nang
!        phsm(i,1) = 3.0/16.0/pi*(1+cos(ang(i)*rad)**2)
!        phsm(i,2) = phsm(i,1)
!        phsm(i,3) = 3.0/8.0/pi*cos(ang(i)*rad)
!        phsm(i,4) = phsm(i,3)
!        phsm(i,5) = 3.0/16.0/pi*(cos(ang(i)*rad)**2-1)
!        phsm(i,6) = 0.
!     enddo
! 06.03.17 miho
     call sdp2(nln,pl,dp,np)
! 06.03.17 end
     !mm220221
     npol = ipol
     if( ipol > 1 ) npol = 6
     !mm231005
     nang0 = nang
     ang0 = ang
  endif
  nang = nang0
  ang = ang0
! initialization end
  do i=1,na0
     am0(i)=cos(th0(i)*rad)
  enddo
  do i=1,na1
     am1u(i)=cos(th1(i)*rad)
  enddo
! wavelength loop start
  sol   =0
! 96.12.27 takayabu
  wg0=0
  do l=1,nln
     do i=1,10
        thk0(l,i)=0
     enddo
  enddo
  do l=1,ntau
     do j=1,na0
        flxd0(j,l)=0
        flxd00(j,l)=0
        flxu0(j,l)=0
        do i=1,na1
           do k=1,nfi
              ai0(i,j,k,l)=0
           enddo
        enddo
     enddo
  enddo
  wgt1=1
  wlc=wl*1.0e-4
  wn=10000/wl
  if(isol>0) then
     fsol=sunir(wn)
   else
     fsol=0
  endif
  if(wl.le.2.0) then
     nplk1=0
   else
     nplk1=2
  endif

!mm220222 gas absorption
  nww = 1
  if(icont(2)==0) then ! no gas absorption
    nch = 1
    tkd(1:nln,1:nch,1) = 0
    wgt(1:nch) = 1.0
!!tr211026
!  nww=1
!  if(icont(2)==0) tkd(1:nln,1:nch,1:nww)=0
  else if(icont(2)==1) then ! broadband
    call abgasx(initg,iug,nln,nch,ngas(1),dp,np,amtpb,tl,wgt,tkd)
  else if(icont(2)==2) then ! narrowband for Rstar standard LUT
    nch = ngas(1)
    dintvl = dble(ngas(2))
    call abgask(iug,nln,nch,dintvl,nww,wl,dp,np,amtpb,tl,tkd,err)
    if(err.ne.'') return 
    call qgausn(gwt,gmu,nch)
    do ik=1,nch
       wgt(ik)=2.d0*gwt(ik)*gmu(ik)
    enddo
  end if
!  nww=1;nch=1;wgt=1.d0;tkd=0.d0
!end gas absorption

! polydispersions
!mm220222
  if(icont(1)==0) then ! no rayleigh
    taur(1:nln) = 0.
    dpf(1:nln)=0. !mm220416
!  if(icont(1)==0) tray=0
  else if(icont(1)==1) then ! rayleigh for rstar
  ! Rayleigh (Frohlich and Shaw, 1980)
    tray=0.00864/wl**(3.916+0.074*wl+0.05/wl)
    do l=1,nln
       taur(l)=abs(prs(ipb(l))-prs(ipb(l+1)))/pstd*tray
  ! added rayleah scattering tau =0 ;tayr=0
  ! to get only aerosol phase function
  ! 2012.11.25
  ! phsf->hsmphsfr chnged
  ! 2012.11.25 by okata
    enddo
    dpf(1:nln)=0. !mm220416
  else if(icont(1)==2) then ! rayleigh for rpstar
  ! Rayleigh (Bodhaine et al.,1999)
    call raysca(nln,wl,co2,dpf,taur)
    do l=1,nln
      taur(l) = taur(l) * (pl(l)/pstd) * (288.15d0/tl(l)) * dzcm(l)
    enddo
  end if
!mm220416
  do l = 1, nln
    g = (1.0-dpf(l))/(1.0+.5*dpf(l)) !mm220707
! rayleigh phase function from rp (currect??)
     do i = 1 ,nang
        phsm(i,l,1) = 3.0/16.0/pi*((1.0+dpf(l))/(1.0-dpf(l))+cos(ang(i)*rad)**2)*g !mm220623
        phsm(i,l,2) = 3.0/16.0/pi*(1+cos(ang(i)*rad)**2)*g
        phsm(i,l,3) = 3.0/8.0/pi*cos(ang(i)*rad)*g
        phsm(i,l,4) = 3.0/8.0/pi*cos(ang(i)*rad)*g*(1.0-2.0*dpf(l))/(1.0-dpf(l)) !tr220627
        phsm(i,l,5) = 3.0/16.0/pi*(cos(ang(i)*rad)**2-1)*g
        phsm(i,l,6) = 0.
     enddo
  end do
  
  taua(1:nln) = 0.
  scaa(1:nln) = 0.
  
! rayleig 5km-120km tayr=0
! 2013.03.12 by okata
! open file phase function !
  do ipoly=1,npoly
     do l=1,nln
        taua3(l)=0
        scaa3(l)=0
        spoly(l,ipoly)=0
        epoly(l,ipoly)=0
!! added spoly,epoly,ppoly
        ppoly(1:nang,l,ipoly,1:npol)=0
     enddo
     do j=1,ncomp(ipoly)
        m=mptc(ipoly,j)
        initrf=1
        do l=1,nln
!080219       if(ifrh.eq.1) then
!080219        rhl=trh
!080219       else
           rhl=rh(l)
!080219       endif
           if(cnprf1(l,ipoly,j).gt.0 .and. initrf.gt.0) then
              !call rfidxb(iuk,wlc,rhl,m,ispcvp,rfracp,asphr, &
              !rop,dryap,nawp,awcrp,nwlv,wlv,rfi,scr,sci, &
              !voldp2,volp2,cextp2,cabsp2,nang,ang,php,err)
              call rfidx7p(iuk,iup,iukd,ipol,ins(m),wlc,rhl,m,&
                    ispcvp,rfracp,asphr,rop,dryap,nawp,awcrp, &
                    nv,nwlv,wlv,rfi,scr,sci, &
                    voldp2,volp2,cextp2,cabsp2,nang,ang,php,err)
              if(err.ne.' ') return
              if(nawp(m).le.0) initrf=0
           endif
           cext1=cextp2/voldp2*cnprf1(l,ipoly,j)
           cabs1=cabsp2/voldp2*cnprf1(l,ipoly,j)
           csca1=cext1-cabs1
           taua3(l)=taua3(l)+cext1
           scaa3(l)=scaa3(l)+csca1
! added spoly,epoly,ppoly
           spoly(l,ipoly)=spoly(l,ipoly)+csca1
           epoly(l,ipoly)=epoly(l,ipoly)+cext1
! added ppoly!!! 120722
           if (cnpt(ipoly).ne.0.) then
              ppoly(1:nang,l,ipoly,1:npol)=ppoly(1:nang,l,ipoly,1:npol)+  &
              php(1:nang,1:npol)/(cextp2-cabsp2)*csca1
            else
              ppoly(1:nang,l,ipoly,1:npol)=0
           endif
        enddo
     enddo
     if(icn.le.0) then
        fac=cnpt(ipoly)
      else if(icn.eq.1) then
        fac=cnpt(ipoly)/volt(ipoly)
      else if(icn.eq.2) then
        fac=cnpt(ipoly)/cextpn(ipoly)
      else
        sum=0
        do l=1,nln
           sum=sum+taua3(l)
        enddo
        if (cnpt(ipoly).ne.0.) then
           fac=cnpt(ipoly)/sum
         else
           fac=0
        endif
     endif

!tr bottom to top
     do l=1,nln
! added spoly,epoly,ppoly
        taua(l)=taua(l)+taua3(l)*fac
        scaa(l)=scaa(l)+scaa3(l)*fac
        if (cnpt(ipoly).ne.0.) then
           ssaa(l)=scaa(l)/taua(l)
         else
           ssaa(l)=0
        endif
!           ssaa1=ssaa(l)
! bug spoly and epoly each ipoly !
! 2012.12.14

        if(epoly(l,ipoly)>0) then
           ssapoly(l,ipoly)=spoly(l,ipoly)/epoly(l,ipoly)
         else
           ssapoly(l,ipoly)=1
        endif
        if(spoly(l,ipoly)>0) then
           ppoly(1:nang,l,ipoly,1:npol)=ppoly(1:nang,l,ipoly,1:npol)/spoly(l,ipoly)
         else
           ppoly(1:nang,l,ipoly,1:npol)=0
        endif
        epoly(l,ipoly)=epoly(l,ipoly)*fac
     enddo
  enddo
! gas absorption
! 06.03.17 changed for the k-distribution using hitran 2004.
!     nww=int(log10(1.d4/wla)*dintvl+1.d0) &
!          -int(log10(1.d4/wlb)*dintvl+1.d0)+1
!     call abgask(iug,nln,nch,dintvl,nww,wlb,dp,np,amtpb,tl,tkd,err)
!     if(err.ne.'') return
!     call qgausn(gwt,gmu,nch)
!     do ik=1,nch
!        wgt(ik)=2.d0*gwt(ik)*gmu(ik)
!     enddo

!tr210905 tau out of this routine

! initialization

! aerdb_rstar distribute at only 1km
! 06.03.24 => changed to n channnels by miho.
!start ik roop
!
! gas absorption
! read from parag

! for cloudsat
  nww=1
  do iww=1,nww
     do ik=1,nch
!!! gas absorption reject
        do l=1,nln
           thk(l)=tkd(l,ik,iww)+taua(l)+taur(l)
           omg(l)=(scaa(l)+taur(l))/thk(l)
! 97. 1.14  20.0->200.0
           taumax=200.0/(1.04-omg(l))
           thk(l)=max(thk(l),1.0e-7)
           thk(l)=min(thk(l),taumax)
        enddo
! indt=1 case
        if(indt.eq.1) then
           lf=0
           dpt=0
           do l=1,nln1
              if(ipbf(l).gt.0) then
                 lf=lf+1
              endif
              dpt=dpt+thk(l)
           enddo
        endif
!mm220222 purge
!! plank function fitting (nplk1=3)
!        bgnd=0
!        if(nplk1.gt.0) then
!!        xx(1)=0
!!        yy(1)=bplnk(wl1,tmp(1))
!!        yy(1)=bplnk(wl1,tmp(ipb(1)))
!!        do 37 l=1,nln
!!        xx(l+1)=xx(l)+thk(l)
!!   37   yy(l+1)=bplnk(wl1,tmp(l+1))
!!   37   yy(l+1)=bplnk(wl1,tmp(ipb(l+1)))
!!        call cspl1(nln1,xx,yy,aa,bb,cc,dd)
!!        do 25 l=1,nln
!!        l1=l+1
!!        cplk(1,l)=aa(l1)+bb(l1)*xx(l)+cc(l1)*xx(l)**2+dd(l1)*xx(l)**3
!!        cplk(2,l)=bb(l1)+2*cc(l1)*xx(l)+3*dd(l1)*xx(l)**2
!!        cplk(3,l)=cc(l1)+3*dd(l1)*xx(l)
!!   25   cplk(4,l)=dd(l1)
!           b2=bplnk(wl,tmp(ipb(1)))
!           do l=1,nln
!              b1=b2
!              b2=bplnk(wl,tmp(ipb(l+1)))
!!                 cplk(1,l)=b1
!!                 cplk(2,l)=(b2-b1)/thk(l)
!           enddo
!           if(indg.eq.0) then
!              bgnd=(1-galb)*bplnk(wl,gtmp)
!            else if(indg.gt.0) then
!              bgnd=bplnk(wl,gtmp)
!           endif
!        endif
!end purge
!
!1205120
!!           call rtrn21(indg,inda,indt,indp,imthd,nda,na1,am1u,na0,am0, &
!!                nfi,fi,nln,thk,omg,nlgn1,g,nang,ang,phsf,epsp,epsu,galb,fsol, &
!!                nplk1,cplk,bgnd,ntau,utau,scr,sci,flxd,flxu,ai,err)
           if (err .ne. ' ') return
!c 06.03.29 miho
        sol    = sol+            fsol*wgt(ik)*wgt1
        wgt2=wgt(ik)*wgt1
!v     wgt average calculation
           wg0=wg0+wgt2
!v     do l=1,nln
!      l9=min(l,nln)
!v     wgt2=wgt(ik)*wgt1
!v     wg0(l) = wg0(l)+              wgt2
! 06.03.29 end
!              thk0(l,1)=thk0(l,1)+  thk (l)*wgt2
!              thk0(l,2)=thk0(l,2)+  taua(l)*wgt2
!              thk0(l,3)=thk0(l,3)+  taur(l)*wgt2
!              sca2=omg(l)*thk(l)*wgt2
!              thk0(l,4)=thk0(l,4)+sca2
!              thk0(l,5)=thk0(l,5)+g(2,l)*sca2
!              thk0(l,6)=thk0(l,6)+g(3,l)*sca2
!           enddo
! k distribution correct ssa and tau input data in monte carlo
! ssa
        lf=0
        thks=0
        do l=1,nln1
           if(l.ge.2) thks=thks+thk(l-1)
!      l9=min(l,nln)
           if(ipbf(l).gt.0) then
              lf=lf+1
              do j=1,na0
!                    flxd0(j,lf)=flxd0(j,lf)+flxd(j,lf)*wgt2
                 fx=fsol*am0(j)*exp(-thks/am0(j))
iw99=0
if(iw99>0) then
write(*,'(3i4,4f12.4,a)') l,lf,j,thks,am0(j),fsol,fx,' l lf j thks am0 fsol fx'
endif
                 flxd00(j,lf)=flxd00(j,lf)+fx*wgt2
!                    flxu0(j,lf)=flxu0(j,lf)+flxu(j,lf)*wgt2
!!                    if(inda.ne.0) then
!                       do i=1,na1
!                          do k=1,nfi
!                          ai0(i,j,k,lf) = ai0(i,j,k,lf) + ai(i,j,k,lf)*wgt2
!                          enddo
!                       enddo
!                    endif
              enddo
!42      continue
           endif
        enddo
!44               continue
!     enddo
!   13 continue
     enddo
  enddo
!  enddo drop iws-loop

!ccenddo
!    7 continue
  sol=sol/wg0
! 06.03.02 miho
!      if(dw.gt.0.) then
!        sol=sol*dw
!    wg0=wg0/dw
!      endif
! 06.03.02 end
  lf=0
  do l=1,nln1
!      l9=min(l,nln)
     if(ipbf(l).gt.0) then
        lf=lf+1
        do j=1,na0
!          flxd0(j,lf)=flxd0(j,lf)/wg0
           flxd00(j,lf)=flxd00(j,lf)/wg0
!           write(*,*) j,lf,flxd00(j,lf)
!           flxu0(j,lf)=flxu0(j,lf)/wg0
!           if(inda.ne.0) then
!              do i=1,na1
!                 do k=1,nfi
!           ai0(i,j,k,lf)=ai0(i,j,k,lf)/wg0
        enddo
!              enddo
     endif
  enddo

  return
end subroutine rstr7cms
