!mm220221
subroutine strnm(icycl,ds,ibs,ps,ibe,pe,nb,bx,dbx,ce, &
     taus,trns,eps,mlim,ind)
! transmission to the boundary of the system.
!--- history
! 88. 9.21  created
!--- input
! icycl      i        if 1 then cyclic boundary condition
! ds      r(3)        direction vector of the ray
! ibs     i(3)        current box location
! ps      r(3)        current location
! nb      i(3)        number of boxes along each axis
! bx     r(kb1,3)     coordinates of box interfaces
! ce   r(kx,ky,kz) extinction cross section (/length) for each airmass
! mlim     i          muximum loop number for checking
! eps      r          convergence criterion
!--- output
! ibe     i(3)        box location at boundary
! pe      r(3)        location at boundary
! taus     r          optical distance to the boundary
! trns     r          transmissivity
! ind      i           -1: hit -x boundary, +1: +x
!                      -2:     -y         ; +2: +y
!                      -3: -z             ; +3: +z
!                      10: too many boxes
!                      11: ill condition for hitting wall.
!--- parameter
! kb1      i          declared size for bx
! kx,ky,kz i          declared numbers of boxes in x/y/z directions
!$endi
  use paras
  implicit none

  integer, intent(in) :: icycl, ibs(3), nb(3), mlim
  real(8), intent(in) :: ds(3), ps(3), bx(kb1,3), dbx(kb1,3,2), ce(kx,ky,kz), eps
  integer, intent(out) :: ibe(3), ind
  real(8), intent(out) :: pe(3), taus, trns
!c working area
!  character(len=64) erc
  real(8) :: ps1(3),pe1(3)
  integer :: ib1(3)
  integer :: i,lim,is

  taus=0
  ps1(1:3) = ps(1:3)
  ib1(1:3) = ibs(1:3)
!c loop for checking boxes
  do lim = 1, mlim
     call boxss(dbx,ds,ib1,ps1,is,pe1,eps)

!     if(erc(1:1).ne.' ') then
     if(is.eq.0) then
        ind=11
        return
     endif
     taus=taus+ce(ib1(1),ib1(2),ib1(3))*sqrt(sum((pe1(1:3)-ps1(1:3))**2))

!c check boundary condition
!c if ibnd.gt.0 then cyclic boundary condition
     call chksd(icycl,is,ib1,pe1,ibe,pe,nb,bx,ind)

     if(ind .ge. 4) then !cc fatal error
        ind=11
        return
     else if(iabs(ind) .eq. 3) then !cc lateral boundaries
!c transmissivity
        trns=exp(-taus)
        return
     else if(iabs(ind) .le. 2 .and. ind .ne. 0 .and. icycl.le.0) then
!c transmissivity
        trns=exp(-taus)
        return
     endif
!cc update starting position for further checking
     ps1(1:3)=pe(1:3)
     ib1(1:3)=ibe(1:3)
  end do
  ind=10
  return
end subroutine strnm