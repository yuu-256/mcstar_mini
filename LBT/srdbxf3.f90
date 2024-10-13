!mm220221
subroutine srdbxf3(init,icycl,indg,dseed,drs,pvw,nb,bx,dbx,ag, &
     ce,omg,angdp,nns,nng,u,ds,ps,ibx,iby,fdg,fdg0,fug, &
     nnsg,ug,dsg,psg,ibxg,ibyg,lim1,mlim,ind)
! forward monte carlo radiative transfer for a box system
! one photon input case
!--- history
! 110616  created from jmonte
!tr210904 modified; drop np
!tr211104 rdbxsf3, ground variables with knng
!tr211106 iexit initialization, Ground reflection by random Ag
!tr211212 fis, photon shooting plane
!mm211219 debugged
!---
!
! use uniform random generator  ggubs
!--- input
! init      i          if <>0 then initialization
! icycl     i          if >0 then cyclic boundary condition
! dseed     d          initial value for random generator
! drs      r(3)        solar direction vector from eyes
! pvw      r(3)        solar insolation position  (x, y, z)
! nb       i(3)        number of boxes along each axis
! bx      r(kb1,3)     coordinates of box interfaces
!                      increasing box number should correspond to increasing coor
! ag     r(kx,ky)      ground albedo, give zero for the no-ground problem
! ce    r(kx,ky,kz)    extinction cross section (/length)
! omg   r(kx,ky,kz)    single scattering albedo
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
! u        r          photon intensity (insolation photon=1)
! ug(knng,2)          dn/up diffuse photons hit ground at nng-th interaction between atmosphere and ground
! ibx,iby  i          emerging x-y point of photon at toa
! ibxg, ibyg(knng)    box location of ground hit at nng-th interaction between atmosphere and ground
! fdg                 downward diffuse horizontal flux at ground !tr
! fdg0                downward direct at ground
! fug                 reflected horizontal flux at ground
! lim1                scattering order
! ds(3),ps(3)         direction and position of photon exit
! dsg(knng,3,2), psg(knng,3)    dn/up directions and location in the box of ground hit at nng-th interaction
! nns                 number of scatterings
! nng                 number of ground hit
!--- parameter
! kx,ky,kz            declared numbers of boxes in x/y/z directions
! kb1                 kb1=max(kx, ky, kz) + 1
! na4                 Number of scattering angles and distribution steps
!$endi
  use paras
  implicit none
   character(len=64) :: erc
  real(8), save :: eps,u1min
  integer, save :: mlim1,mlim2
  real(8) ::  dseed
  real(8) ::  drs(3),pvw(3),bx(kb1,3),dbx(kb1,3,2), &
       ce(kx,ky,kz),ag(kx,ky)
  integer :: nb(3),mlim(*)!,ibotm(2),iside(4,3)
  integer :: init
  real(8) :: u,ug(knng,2)
  integer :: lim2,lim1
  real(8) :: tau,dtau,taue
  integer :: ibx,iby,ibz,ibxg(knng),ibyg(knng)
  integer :: ind
  real(8) :: cet,dst,agg,cst,w0t
  integer :: i,is
  integer :: iexit,icycl,indg(kx,ky)
  real(8) :: trne,rag

! working area
  real(8) ::  ps(3),pe(3),ds(3), &
       pe2(3),cpc(3),dsg(knng,3,2),psg(knng,3)
  integer :: ibs(3),ibe2(3)
   real(8) :: angdp(na4,kx,ky,kz)
  real(8) :: omg(kx,ky,kz)
!  integer :: i1,j,k
!tr
  real(8):: fdg,fdg0,fug,rndv(3),scat
  integer:: nrnd,iw99,nx,ny,nz,nns,nng,nnsg(knng)

! initialization

  if (init.ne.0) then
     init=0
     call cpcon(cpc)
     eps=cpc(1)*10.0
     u1min=10.0**(cpc(2)*0.8)
     mlim1=mlim(1)
     mlim2=mlim(2)
  endif

! start
!  np1=np
  ind=0
  iexit=0
  nx=nb(1); ny=nb(2); nz=nb(3)
  u=1
  fdg=0; fdg0=0; fug=0
  nns=0
  nng=0
! initial point
  ps(1:3)=pvw(1:3)
  ds(1:3)=drs(1:3)
! find box, ibs: found position
  call fbox(ps,bx,nb,ibs,erc)

  if (erc.ne.' ') then
     ind=9
     return
  endif

! extinction of the ray
  lim1=0
  do10 : do
     lim1=lim1+1
!     write(11,*) lim1,'-th order scattering'
     if (lim1.gt.mlim1) then
        ind=8
        return
     endif
     call gttau(dseed,trne,taue)
!tr
     if(trne.le.0) then
        cycle do10
!        ind=9
!        return
     endif

! ray tracing along the incident direction
! ray trace in a box (number of intercepting walls)
     tau=0
     lim2=0
     do11 :do
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

!        if (erc(1:1).ne.' ') then
        if (is.eq.0) then
           ind=9
           return
        endif
!tr
        cet=ce(ibx,iby,ibz)
        cst=omg(ibx,iby,ibz)*ce(ibx,iby,ibz)
        w0t=cst/cet

        if (cet.gt.0) then
          dst=sqrt(sum((pe(1:3)-ps(1:3))**2))
          dtau=cet*dst
! exit to next scattering process
          if (taue <= tau+dtau) exit do11 !tr

! next box for continued ray tracing
! if icycl.gt.0 then cyclic boundary for lateral sides
          tau=tau+dtau
        end if
        call chksd(icycl,is,ibs,pe,ibe2,pe2,nb,bx,ind)

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

!tr downward photon being hitted by ground
!tr test           fdg=fdg+u*abs(ds(3))  !tr for incident horizontal flux at ground
           nng=nng+1
           if(nng>knng) then
             ind=7
             return
           endif
           nnsg(nng)=nns

           if(nns==0 .and. nng==1) then
             fdg0=u
             ug(nng,1)=0
            else
             fdg=fdg+u
             ug(nng,1)=u  ! diffuse downward
           endif
           ibxg(nng)=ibx; ibyg(nng)=iby
           dsg(nng,1:3,1)=ds(1:3)
           psg(nng,1:3)=ps(1:3)
           
           if(indg(ibs(1),ibs(2))==0) then !lambertian surface
             agg=ag(ibs(1),ibs(2))
           else
             !agg=
           end if
           if(agg <= 0) return
! random reflection method
           nrnd=1
           call rdmu2(dseed,nrnd,rndv)

           if(rndv(1) > agg) then
             ug(nng,2)=0
             return
           endif

           call ggdir(dseed,ds)
           if(indg(ibs(1),ibs(2))/=0) then !except lambertian surface
             continue
           end if
           ibs(3)=1
!           u=u*agg   ! no need in the random reflection method
           fug=fug+u
           ug(nng,2)=u
           dsg(nng,1:3,2)=ds(1:3)
         else
          ind=9
          return
        endif  ! ind branch end

        if (iexit > 0) return

        if(u.le.u1min) then
           u=0
           ind=0
           return
        endif
     enddo do11
!! scattering point
! decide scattering or absorbing
!tr
     cet=ce(ibx,iby,ibz)
     cst=omg(ibx,iby,ibz)*ce(ibx,iby,ibz)
     w0t=cst/cet

     dst=(taue-tau)/cet
     ps(1:3)=dst*ds(1:3)+ps(1:3)    !tr
     nns=nns+1
     u=w0t*u

     if (u.le.u1min) then
        u=0
        ind=0
        return
     endif
     if(nns+1>knns) then
        ind=7
        return
     end if
! next scattering direction
     nrnd=2
     call rdmu2(dseed,nrnd,rndv)
     !mm231007 interpolation
     !!mm211219 debugged
     !!i=rndv(1)*(na4-1)+1
     !i=rndv(1)*(na4-1)+1.5
     !if(i>=na4) i=na4
     !scat=angdp(i,ibx,iby,ibz)
     call angintp(rndv(1),na4,angdp(1:na4,ibx,iby,ibz),ind,scat)
     if(ind.eq.9) return
     !end interpolation

     call atrns(scat,rndv(2)*360.0,ds,eps)

  enddo do10
! photon quench
  return
end subroutine srdbxf3
