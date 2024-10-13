!mm220221
subroutine chksd(icycl,is,ib,pe,ibu,peu,nb,bx,ind)
! check emission from a boundary of the system
! is, ib and pe are given by boxss routine before this routine
! icycl       i       if >0 then cyclic lateral boundary condition
! is          i       surface of emergence (+- 1, 2, 3 = x, y, z)
!                     >0 then positive x,y,z side relative to the box center
!                     <0 then negative x,y,z side relative to the box center
! ib        i(3)      current box location.
! pe        r(3)      emergent point
! nb        i(3)      number of boxes along each axis
! bx     r(kb1,3)     coordinates of box interfaces
!--- output
! ibu       i(3)      box location + increment
!                     taking the boundary condition into account
! peu       r(3)      updated location
! ind         i       -1: emit from -x boundary, +1: from +x
!                     -2: -y; +2: +y;  -3: -z;  +3: +z
!                     0: inside the system.
!                     9: error.
!--- parameter
! kb1        i       max(kx,ky,kz)+1.
!$endi
!      parameter (kb1   =11)
  use paras
  implicit none
  integer, intent(in) :: icycl, is, ib(3), nb(3)
  real(8), intent(in) :: pe(3), bx(kb1,3)
  integer, intent(out) :: ibu(3), ind
  real(8), intent(out) ::  peu(3)
  integer :: k

! new box location
  peu(1:3)=pe(1:3)
  ibu(1:3)=ib(1:3)
  k=iabs(is)
  if(k.lt.1 .or. k.gt.3) then
     ind=9
     return
  endif
  ibu(k)=ib(k)+is/k
  if (ibu(k) .gt. 0 .and. ibu(k) .le. nb(k)) then
     ind = 0
     return
  else if(k .eq. 3) then ! check top
     if(ibu(3).gt.nb(3)) then ! top
        ind=3
        return
     else if(ibu(3).le.0) then ! bottom
        ibu(3)=1    !tr
        ind=-3
        return
     endif
  else if(k .le. 2) then ! check lateral sides of the system
     if(ibu(k) .le. 0) then
        ind=-k
        if(icycl.gt.0) then
           ibu(k)=nb(k)
           peu(k)=bx(nb(k)+1,k)
        endif
        return
     else if(ibu(k).gt.nb(k)) then
        ind=k
        if(icycl.gt.0) then
           ibu(k)=1
           peu(k)=bx(1,k)
        endif
        return
     endif
  endif
end subroutine chksd