!mm220221
subroutine rdmu1(r0,n,r2)
!c random number generator
!c--- history
!c 88.10.14 created
!c--- input
!c r0    d    initial condition (0 to 1)
!c n     i    number of random numbers
!c--- output
!c r    r(n)  uniform random numbers ( 0 to 1)
!c$endi
  use paras
  implicit none

  real(8) ::  r2(n)
  real(8) :: r0,r1
  integer :: i,n
  do i=1,n
     r0=(pi+r0)**5
     r1=int(r0)
     r0=r0-r1
     r2(i)=r0
  enddo
  return
end subroutine rdmu1