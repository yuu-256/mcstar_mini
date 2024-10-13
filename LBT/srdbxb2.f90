!mm220221
subroutine srdbxb2(init,ipol,icycl,indg,dseed,drs,dvw,pvw,nb,bx,dbx,ag, &
    ce,omg,phs,angdp,nns,nng,u1,u1t,utrn0,lim1,mlim,ind)

! radiative transfer in box system.
! inverse monte-carlo. one photon case.
!--- history
!tr211111 created
!tr211212 deleted if(ind1 /= 3)
!mm211218 added output u1t for tms
!mm211219 debugged
!mm220118 interpolation of angdp
!--- input
! init      i          if <>0 then initialization
! ipol      i          polarization state: 1: scalar, 3: semi vector, 4: full vector
! icycl     i          if >0 then cyclic boundary condition
! dseed     d          initial value for random generator
! drs      r(3)        solar direction vector from eyes
! dvw      r(3)        viewing direction
! pvw      r(3)        receiver position  (x, y, z)
! nb       i(3)        number of boxes along each axis
! bx      r(kb1,3)     coordinates of box interfaces
!                      increasing box number should correspond to increasing coor
! ag     r(kx,ky)      ground albedo, give zero for the no-ground problem
! ce    r(kx,ky,kz)    extinction cross section (/length)
! omg   r(kx,ky,kz)    single scattering albedo
! phs   r(na4,kx,ky,kz) phase function
! angdp r(na4,kx,ky,kz) inverse distribution from accumulated(P) to Ang
!                      0 (forward)  1 (backward)
! mlim     i(2)        1: muximum loop number for scattering order
!                      2: muximum loop number for intercepting walls
!--- output
! init     i          0
! ind                 0: photon is extincted as u=0
!              +-1, +-2: exitted from a lateral side
!                     3: exitted from toa
!                    -3: exitted from boa when ag=0
!                   >=4:  error in chksd
!                     7: too many scattering
!                     8: too many ray tracing
!                     9: fatal error
! nns                 number of scatterings
! nng                 number of ground hit
! u        r          photon intensity (insolation photon=1)
! ug(knng,2)          dn/up diffuse photons hit ground at nng-th interaction between atmosphere and ground
! lim1                scattering order
!--- parameter
! kx,ky,kz            declared numbers of boxes in x/y/z directions
! kb1                 kb1=max(kx, ky, kz) + 1
! knang               declared number of scattering angles

  use paras
  implicit none
  character(len=64) :: erc
  real(8), save ::  eps,u1min,da
  integer, save :: mlim1,mlim2
!  double precision dseed
  real(8) :: dseed
  real(8) :: drs(3),dvw(3),pvw(3),bx(kb1,3),dbx(kb1,3,2),ag(kx,ky)
  integer :: nb(3),mlim(*), ipol
!c working area
  real(8) :: ps(3),pe(3),ds(3),pe2(3), &
       cpc(3)
  integer :: ibs(3),ibe2(3)
  integer :: init,ind,ibx,iby,ibz,i, &
       icycl,indg(kx,ky),lim1,lim2,is
  real(8) :: u(kpol2,kpol2),u1(kpol),trne,agg,tau, &
       cet,taue,dtau,w0t
  real(8) :: angdp(na4,kx,ky,kz),phs(na4,kx,ky,kz,kpol2)
  integer :: nns,nng,nrnd,nx,ny,nz
  real(8) :: taus1,trns1,ce(kx,ky,kz),omg(kx,ky,kz),scat
  real(8) :: rndv(2)
  integer:: ind1
  !mm211218 for tms
  real(8) :: u1t(kpol)
  !mm231002 for truncated order
  real(8) :: utrn(kpol), utrn0(kpol)
  integer :: jpol
  real(8) :: ds2(3), phmx(kpol,kpol), rag


! initialization
  if(init.ne.0) then
     init=0
     call cpcon(cpc)
     eps=cpc(1)*100.0
     u1min=10.0**(cpc(2)*0.8)
     mlim1=mlim(1)
     mlim2=mlim(2)
     da=180.0/(na4-1)
  endif

! start
  ind=0
  nx=nb(1); ny=nb(2); nz=nb(3)
  u(1:ipol,1:ipol)=0.
  do jpol = 1, ipol
    u(jpol,jpol)=1
  end do
  u1=0
  u1t=0 !mm211218 for tms
  utrn=0.
  utrn0=0.
  nns=0
  nng=0
! backward tracing from the receiver
  ps(1:3)=pvw(1:3)
  ds(1:3)=dvw(1:3)
! find box, ibs: found position
  call fbox(ps,bx,nb,ibs,erc)
  if (erc.ne.' ') then
    ind=9
    return
  endif

! extinction of the ray
  do10 :  do lim1 = 1, mlim1
!c     write(11,*) lim1,'-th order scattering'
    call gttau(dseed,trne,taue)
    if(trne.le.0) then
      cycle do10
    endif

! ray tracing along the incident direction
! ray trace in a box (number of intercepting walls)
    tau=0
    lim2=0
    do11 : do
      if(u(1,1).le.u1min) then
         utrn0(1:ipol) = utrn(1:ipol)
         return
      end if

      lim2=lim2+1
!     write(11,*) lim2,'-th wall intercepting'
      if(lim2.gt.mlim2) then
        ind=8
        return
      endif
      ibx=ibs(1)
      iby=ibs(2)
      ibz=ibs(3)

      call boxss(dbx,ds,ibs,ps,is,pe,eps)

!      if (erc(1:1).ne.' ') then
      if (is.eq.0) then
        ind=9
        return
      endif
!
      cet=ce(ibx,iby,ibz)
      w0t=omg(ibx,iby,ibz)

      if(cet .gt. 0) then
         dtau=cet*sqrt(sum((pe(1:3)-ps(1:3))**2))
! exit to next scattering process
         if(taue<=tau+dtau) exit do11

! next box for continued ray tracing
! if icycl.gt.0 then cyclic boundary for lateral sides
         tau=tau+dtau
      end if
      call chksd(icycl,is,ibs,pe,ibe2,pe2,nb,bx,ind)

        if (ind.ge.4) then
           ind=9
           return
        endif
        ps(1:3)=pe2(1:3)
        ibs(1:3)=ibe2(1:3)

        if (iabs(ind).le.2 .and. ind .ne. 0 .and. icycl.le.0 ) return
        if (ind.eq.3) return
! hit ground surface
        if (ind.eq.-3) then

!tr downward photon being hitted by ground
!tr test           fdg=fdg+u*abs(ds(3))  !tr for incident horizontal flux at ground
           nng=nng+1
           if(nng>knng) then
             utrn0(1:ipol) = utrn(1:ipol)
             ind=7
             return
           endif
           if(indg(ibs(1),ibs(2))==0) then !lambertian surface
             agg=ag(ibs(1),ibs(2))
             rag = agg*abs(drs(3))/pi
           else
             !set ag and rag for other case
             continue
           end if
           if(agg <= 0) return
! backward tracing to the solar direction
           call strnm(icycl,drs,ibs,ps,ibe2,pe2,nb,bx,dbx, &
             ce,taus1,trns1,eps,mlim2,ind1)
!tr211212  if(ind1 /= 3) then
!             ind=9
!             return
!           endif

           utrn(1:ipol)=trns1*u(1:ipol,1)*rag
!           if(ipol/=1) then ! when lambertian surface
!             utrn(2:ipol) = 0
!             utrn(1:ipol)=matmul(trns1*u(1:ipol,1:ipol)*abs(drs(3))/pi,fsol(1:ipol))
!           end if

!tr211116           u1=u1+trns1*u/pi
           u1(1:ipol)=u1(1:ipol)+utrn(1:ipol)

           call ggdir(dseed,ds)
           if(indg(ibs(1),ibs(2))==0) then !lambertian surface
             u(1:ipol,1)=u(1:ipol,1)*agg
             if(ipol/=1) u(1:ipol,2:ipol) = 0.0 ! when lambertian surface
           else
             ! weight function with lambertian
             continue
           end if
           ibs(3)=1
        endif  ! ind branch end
!  u=u*agg   ! no need in the randome reflection method

     enddo do11
!! scattering point
!tr
!     w0t=omg(ibx,iby,ibz)

     ps(1:3)=(taue-tau)/cet*ds(1:3)+ps(1:3)
     nns=nns+1
     u(1:ipol,1:ipol)=w0t*u(1:ipol,1:ipol)
! ray tracing to the solar direction
     call strnm(icycl,drs,ibs,ps,ibe2,pe2,nb,bx,dbx, &
         ce,taus1,trns1,eps,mlim2,ind1)

!tr211212  if(ind1 /= 3) then
!       ind=9
!       return
!     endif

     i=acos(sum(drs(1:3)*ds(1:3)))/rad/da+1.5
     if(i<1 .or. i>na4) then
       ind=9
       return
     endif
     if(ipol==1) then
       utrn(1)=trns1*u(1,1)*phs(i,ibx,iby,ibz,1)
     else
       ! changed to matrix operation in phs
       call scatmat(ipol,phs(i,ibx,iby,ibz,1:kpol2),ds,drs,phmx,eps)
       phmx(1:ipol,1:ipol)=trns1*matmul(u(1:ipol,1:ipol),phmx(1:ipol,1:ipol))
       utrn(1:ipol)=phmx(1:ipol,1)
     end if
     u1(1:ipol)=u1(1:ipol)+utrn(1:ipol)
     !mm211214 for tms
     if(nng+nns==1) u1t(1:ipol)=utrn(1:ipol)
     if(nns+1>knns) then
        utrn0(1:ipol) = utrn(1:ipol)
        ind=7
        return
     end if
! next scattering direction
     nrnd=2
     call rdmu2(dseed,nrnd,rndv)
     !mm220118 interpolation
     !!mm211219 debugged
     !!i=rndv(1)*(na4-1)+1
     !i=rndv(1)*(na4-1)+1.5
     !if(i<1 .or. i>na4) then
     !  ind=9
     !  return
     !endif
     !scat=angdp(i,ibx,iby,ibz)
     call angintp(rndv(1),na4,angdp(1:na4,ibx,iby,ibz),ind,scat)
     if(ind.eq.9) return
     !end interpolation

     ds2(1:3) = ds(1:3)
     call atrns(scat,rndv(2)*360.0,ds,eps)
     ! change polarization matrix after here
     if(ipol==1) cycle
     i=int(scat/da+1.5)
     if(i>=na4 .or. i<=1) cycle
     call scatmat(ipol,phs(i,ibx,iby,ibz,1:kpol2)/phs(i,ibx,iby,ibz,1),ds2,ds,phmx,eps)
     !u(1:ipol,1:ipol) = u(1:ipol,1:ipol)*phmx(1,1) !matmul(u(1:ipol,1:ipol),phmx(1:ipol,1:ipol))
     u(1:ipol,1:ipol) = matmul(u(1:ipol,1:ipol),phmx(1:ipol,1:ipol))

!     do jpol = 1, ipol
!     write(*,'(A,1p4e12.4)') "phmx:", u(jpol,1:ipol)
!     end do
!     return
  enddo do10
  ind=8
  return
! photon quench
  return
end subroutine srdbxb2