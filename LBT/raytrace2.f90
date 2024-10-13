!mm220221
subroutine raytrace2(icycl,nb,bx,dbx,ds,ps0,ibs0,ps,ibs,&
  ce,omg,taue,tau,dtau,taus,dtaus,eps,mlim2,ind)
!mm220109 created from raytrace
  use paras, only: kx, ky, kz, kb1
  implicit none
  ! I/O
  integer :: icycl, nb(3)
  real(8) :: bx(kb1,3), dbx(kb1,3,2), ds(3)
  real(8) :: ps0(3), ps(3)
  integer :: ibs0(3), ibs(3)
  real(8) :: ce(kx,ky,kz), omg(kx,ky,kz)
  real(8) :: taue, tau, dtau, taus, dtaus
  real(8) :: eps
  integer :: mlim2, ind
  ! work
  real(8) :: cet, dst, sca
  integer :: lim2, ibx, iby, ibz, is
  real(8) :: pe(3), pe2(3)
  integer :: ibe(3), ibe2(3)
!  character (len=64) :: erc
  
  ps(1:3) = ps0(1:3)
  ibs(1:3) = ibs0(1:3)
!  erc = ' '
  
! ray tracing along the incident direction
! ray trace in a box (number of intercepting walls)
  ind=0 !mm220131
  tau=0
  taus=0
  lim2=0
  do11 : do
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
!    if (erc(1:1).ne.' ') then
    if (is.eq.0) then
      ind=9
      return
    endif
!
    cet=ce(ibx,iby,ibz)
    sca=omg(ibx,iby,ibz)*cet

    if(cet .le. 0) goto 9
    dst=sqrt((pe(1)-ps(1))**2 &
         +(pe(2)-ps(2))**2+(pe(3)-ps(3))**2)
    dtau=cet*dst
    dtaus=sca*dst
! exit to next scattering process
    if(taue<=tau+dtau) exit do11

! next box for continued ray tracing
! if icycl.gt.0 then cyclic boundary for lateral sides
9   call chksd(icycl,is,ibs,pe,ibe2,pe2,nb,bx,ind)

    if (iabs(ind).ge.4) then
      ind=9
      return
    endif
    ps(1:3)=pe2(1:3)
    ibs(1:3)=ibe2(1:3)

    if (iabs(ind).le.2) then
      if ((iabs(ind).eq.0).or.(icycl.gt.0)) then
        tau=tau+dtau !mm220109 debugged
        taus=taus+dtaus
        cycle
      end if
    end if
    return
  enddo do11
end subroutine raytrace2