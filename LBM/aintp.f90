function aintp(xx,n,x,y)
! linear interpolation
!--- history
! 90. 3. 6 created
!--- input
! xx    r     interpolation point
! n     r     nbr of data
! x   r(nn)   independent variable,  x(i) .lt. x(i+1)
! y   r(nn)   dependent   variable
!--- output
! aintp r     output
!--
  implicit none
  real(8),intent(in):: xx
  integer,intent(in):: n
  real(8),intent(in):: x(n)
  real(8),intent(in):: y(n)
  real(8):: aintp

  integer:: i,k1,k2
!--exec
! 1 point
  if(n<=1) then
     aintp=y(1)
     return
  endif
! 2 points or xx.le.x(1)
  if(xx<=x(1) .or. n==2) then
     aintp=y(1)+(y(2)-y(1))*(xx-x(1))/(x(2)-x(1))
     return
  endif
! xx.ge.xx(n)
  if(xx>=x(n)) then
     aintp=y(n-1)+(y(n)-y(n-1))*(xx-x(n-1))/(x(n)-x(n-1))
     return
  endif

! more than 2 points
  do i=1,n
     k1=i-1
     if(xx<=x(i)) exit
  enddo
  k1=max(1,k1)
  k2=k1+1
  if(k2>n) then
     k2=n
     k1=k2-1
  endif
  aintp=y(k1)+(y(k2)-y(k1))*(xx-x(k1))/(x(k2)-x(k1))
  return
end function aintp
