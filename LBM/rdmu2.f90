!mm220221
subroutine rdmu2(s0,n,r)
!c$name rdmu2
!c random number generator
!c--- history
!c 88.10.14 created
!c 11.06.16 first value always close to 0, so several
!c          runs before use
!c--- input
!c s0    d    initial condition (ge. 1)
!c n     i    number of random numbers
!c--- output
!c r    r(n)  uniform random numbers ( 0 to 1)
!c$endi

  implicit none

  integer :: i
  integer, intent(in) :: n
  real(8), intent(out) :: r(n)
  real(8), intent(inout) :: s0

  real(8), parameter :: a=7.0d0**5, b=2.0d0**31-1.0d0
  do  i=1,n
     s0=a*s0
     s0=mod(s0,b)
     r(i)=s0/(b+1)
  enddo
  return
end subroutine rdmu2