!mm220221
subroutine tmsbxb2(init,ipol,icycl,drs,dvw,pvw,nb,bx,dbx, &
    ce,omg,fp,nang,kna,ang,phs,u,ng0,p,mlim,ind)
! single scattering in box system.
! inverse monte-carlo. one photon case.
!--- history
!mm211216 created from srdbxb1 by terry
!mm211218 validated with ss-approx in Rstar
!mm211218 ibound => icycl, rd => rad
!mm211220 add p for ss correction in actual space (-1: tms)
!--- parameter
! kx,ky,kz            declared numbers of boxes in x/y/z directions
! kb1                 kb1=max(kx, ky, kz) + 1
! knang               declared number of scattering angles
  use paras
  implicit none
  ! I/O
  character(len=64) :: erc
  integer :: init, mlim(2), ipol
  integer :: icycl
  integer :: nb(3)
  real(8) :: bx(kb1,3), dbx(kb1,3,2)
  real(8) :: ce(kx,ky,kz), omg(kx,ky,kz), fp(kx,ky,kz)
  integer :: nang, kna
  real(8) :: ang(kna)
  real(8) :: phs(kna,kx,ky,kz,kpol2)
  real(8) :: u(kpol)
  integer :: ind
  integer :: ng0
  real(8) :: p
  !integer :: nx,ny, nz
  real(8) :: pvw(3), dvw(3), drs(3)
  ! work
  integer, save :: mlim1, mlim2
  real(8), save :: cpc(3), eps
  real(8) :: ps(3), ds(3), pe(3), pe2(3), ps2(3)
  integer :: ibs(3), ibe2(3), ibx, iby, ibz, is
  real(8) :: ang1, phss(kpol2)
  real(8) :: cet, w0t
  real(8) :: dst, dtau, taus, trne, trns, tau
  integer :: initsc = 1
  ! function
  real(8) :: pint4c
  ! gaussian quadrature in box
  integer, parameter :: kgs = 200 !mm220223 10 => 200
  integer :: initg = 1, ng = kgs, ig, jpol, npol
  real(8), save :: gwt(kgs), gmu(kgs)
  real(8) :: u1(kpol2), phmx(kpol,kpol)
  
  npol = ipol
  if(ipol>1) npol = 6

  if( init == 1 ) then
    init = 0
    mlim1 = mlim(1)
    mlim2 = mlim(2)
  end if

  if( initg == 1 .or. ng /= ng0 ) then
    initg = 0
    ng = ng0
    call qgausn( gwt, gmu, ng )
    call cpcon(cpc)
    eps=cpc(1)*100.0
  end if

! start
! initialize
  ind=0
  u1=0
  tau = 0.d0
! backward tracing from the receiver
  ps(1:3)=pvw(1:3)
  ds(1:3)=dvw(1:3)
! find box, ibs: found position
  call fbox(ps,bx,nb,ibs,erc)
  if (erc.ne.' ') then
    ind=9
    u(1:ipol) = 0.0
    return
  endif
  
  ! scattering angle
  if( nang > 0 ) then
    ang1 = acos( sum(ds(1:3)*drs(1:3)) )/rad
    initsc = 1
  end if

  do10 :  do
    ibx=ibs(1)
    iby=ibs(2)
    ibz=ibs(3)
    call boxss(dbx,ds,ibs,ps,is,pe,eps)
!    if (erc(1:1).ne.' ') then
    if (is.eq.0) then
      ind=9
      if(ipol==1) then
        u(1) = u1(1)
      else
        call scatmat(ipol,u1(1:kpol2),ds,drs,phmx,eps)
        u(1:ipol) = phmx(1:ipol,1)
      end if
      return
    endif
!
    cet=ce(ibx,iby,ibz)
    w0t=omg(ibx,iby,ibz)
    if( nang > 0 ) then
      do jpol = 1, npol
        phss(jpol) = pint4c(initsc,ang1,nang,ang,phs(1:nang,ibx,iby,ibz,jpol))
      end do
    else
      phss(1:npol) = phs(iabs(nang),ibx,iby,ibz,1:npol)
    end if

    if(cet .le. 0) goto 9
    dst=sqrt((pe(1)-ps(1))**2 &
         +(pe(2)-ps(2))**2+(pe(3)-ps(3))**2)
    dtau=cet*dst
    
! quadrature in box
    do ig = 1, ng
      ps2(1:3)=dst*gmu(ig)*ds(1:3)+ps(1:3)
      ! scat - sun => trns
      call strnm(icycl,drs,ibs,ps2,ibe2,pe2,nb,bx,dbx, &
             ce,taus,trns,eps,mlim2,ind)
      trne=exp(-(tau+dtau*gmu(ig)))
      u1(1:npol)=u1(1:npol)+gwt(ig)*dtau*w0t*phss(1:npol)*trne*trns*((1.0-fp(ibx,iby,ibz))**p)
    end do

! next box for continued ray tracing
! if icycl.gt.0 then cyclic boundary for lateral sides
    tau=tau+dtau
9   call chksd(icycl,is,ibs,pe,ibe2,pe2,nb,bx,ind)

    if (iabs(ind).ge.4) then
      ind=9
      if(ipol==1) then
        u(1) = u1(1)
      else
        call scatmat(ipol,u1(1:kpol2),ds,drs,phmx,eps)
        u(1:ipol) = phmx(1:ipol,1)
      end if
      return
    endif
    ps(1:3)=pe2(1:3)
    ibs(1:3)=ibe2(1:3)

    if (iabs(ind).le.2) then
      if ((iabs(ind).eq.0).or.(icycl.gt.0)) cycle
    end if
    if(ipol==1) then
      u(1) = u1(1)
    else
      call scatmat(ipol,u1(1:kpol2),ds,drs,phmx,eps)
      u(1:ipol) = phmx(1:ipol,1)
    end if
    return
  enddo do10
! photon quench
  return
end subroutine tmsbxb2