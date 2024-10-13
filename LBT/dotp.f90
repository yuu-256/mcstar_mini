!mm220221
subroutine dotp(w,u,v,ni)
!c w = u*v
!c--- history
!c 88. 9.16  created
!c--- input
!c u       r(ni)      source 1-dim array u.
!c v       r(ni)      source 1-dim array v.
!c ni        i        w(i) = u(i)*v(i),  i=1,ni
!c--- output
!c w       r          w=u*v
!c$endi
  implicit none

  integer :: i,ni
  real(8) :: u(ni),v(ni)
  real(8) :: w

  w=0
  do i=1,ni
     w=w+u(i)*v(i)
  enddo
  return
end subroutine dotp
