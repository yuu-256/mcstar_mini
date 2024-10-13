!mm220221
subroutine erdbxb1(init,icycl,dseed,dvw,pvw,nb,bx,nplk1,cplkz,bgnd,ag, &
    ce,omg,angdp,nns,nng,u1,lim1,mlim,ind)

! radiative transfer in box system.
! inverse Monte-Carlo. one photon case.
! thermal emssion
!--- history
!tr211111 created
!mm211219 debugged
!--- input
! init      i          if <>0 then initialization
! icycl     i          if >0 then cyclic boundary condition
! dseed     d          initial value for random generator
! dvw      r(3)        viewing direction
! pvw      r(3)        receiver position  (x, y, z)
! nb       i(3)        number of boxes along each axis
! bx      r(kb1,3)     coordinates of box interfaces
!                      increasing box number should correspond to increasing coor
! cplkz               B explansion coeff with z
! bgnd                 B at ground
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
  real(8) :: dvw(3),pvw(3),bx(kb1,3),ag(kx,ky)
  integer :: nb(3),mlim(*)
!c working area
  integer :: ibs(3),ibe2(3)
  integer :: init,ind,ibx,iby,ibz,i,icycl,lim1,lim2,is
  integer :: nns,nng,iexit,nrnd,nx,ny,nz,nplk1
  real(8) :: ps(3),pe(3),ds(3),pe2(3),cpc(3)
  real(8) :: u,u1,trne,agg,tau, &
       cst,cet,dst,taue,dtau,w0t
  real(8) :: angdp(na4,kx,ky,kz)
  real(8) :: ce(kx,ky,kz),omg(kx,ky,kz),scat,sfi
  real(8) :: rndv(2),bgnd,cplkz(kplk1,kz),bp
!tr
integer::iw99


! initialization
  if(init.ne.0) then
     init=0
     call cpcon(cpc)
     eps=cpc(1)*100.0
!     u1min=10.0**(cpc(2)*0.8)
     u1min=1.0e-6
     mlim1=mlim(1)
     mlim2=mlim(2)
     da=180.0/(na4-1)
  endif

! start

  ind=0
  iexit=0
  nx=nb(1); ny=nb(2); nz=nb(3)
  u=1
  u1=0
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
  lim1=0
  do10 :  do
    lim1=lim1+1
!c     write(11,*) lim1,'-th order scattering'

    if (lim1.gt.mlim1) then
      ind=8
      return
    endif
    call gttau(dseed,trne,taue)

    if(trne.le.0) then
      cycle do10
    endif



! ray tracing along the incident direction
! ray trace in a box (number of intercepting walls)
    tau=0
    lim2=0
    do11 : do
      if(u.le.u1min) return

      lim2=lim2+1
!     write(11,*) lim2,'-th wall intercepting'
      if(lim2.gt.mlim2) then
        ind=8
        return
      endif
      ibx=ibs(1)
      iby=ibs(2)
      ibz=ibs(3)

      call boxss(bx,ds,ibs,ps,is,pe,eps,erc)

      if (erc(1:1).ne.' ') then
        ind=9
        return
      endif
!
      cet=ce(ibx,iby,ibz)
      cst=omg(ibx,iby,ibz)*ce(ibx,iby,ibz)
      w0t=cst/cet

      if(cet .le. 0) goto 9
      dst=sqrt((pe(1)-ps(1))**2 &
           +(pe(2)-ps(2))**2+(pe(3)-ps(3))**2)
      dtau=cet*dst
! exit to next scattering process
      if(taue<=tau+dtau) exit do11

! next box for continued ray tracing
! if icycl.gt.0 then cyclic boundary for lateral sides
      tau=tau+dtau
9     call chksd(icycl,is,ibs,pe,ibe2,pe2,nb,bx,ind)

        if (iabs(ind).ge.4) then
           ind=9
           return
        endif
        ps(1:3)=pe2(1:3)
        ibs(1:3)=ibe2(1:3)

        if (iabs(ind).le.2) then
           if ((iabs(ind).eq.0).or.(icycl.gt.0)) then
              iexit=0
           else
              iexit=1
           endif
         else if (ind.eq.3) then
           iexit=1   !tr added
           return
! hit ground surface
         else if (ind.eq.-3) then
           agg=ag(ibs(1),ibs(2))

!tr downward photon being hitted by ground
!tr test           fdg=fdg+u*abs(ds(3))  !tr for incident horizontal flux at ground
           nng=nng+1
           if(nng>knng) then
             ind=8
             return
           endif
           u1=u1+u*(1-agg)*bgnd

           if(agg <= 0) return

           u=agg*u
           call ggdir(dseed,ds)
           ibs(3)=1

          else
           ind=9
           return
        endif  ! ind branch end

        if (iexit > 0) return
     enddo do11
! scattering point
     w0t=omg(ibx,iby,ibz)

     dst=(taue-tau)/cet
     ps(1:3)=dst*ds(1:3)+ps(1:3)
     nns=nns+1
! thermal radiation to the receiver
     if(nplk1<=0) then
       bp=0
      else
       bp=cplkz(1,ibz)+cplkz(2,ibz)*(bx(ibz+1,3)-ps(3))
     endif
     if(bp<0) then
       ind=9; return
     endif
     u1=u1+(1-w0t)*u*bp

! next scattering direction
     u=w0t*u
     nrnd=2
     call rdmu2(dseed,nrnd,rndv)
     !mm211219 debugged
     !i=rndv(1)*(na4-1)+1
     i=rndv(1)*(na4-1)+1.5
     if(i<1 .or. i>na4) then
       ind=9
       return
     endif
     scat=angdp(i,ibx,iby,ibz)
     sfi=rndv(2)*360.0

     call atrns(scat,sfi,ds,eps)

  enddo do10
! photon quench
  return
end subroutine erdbxb1